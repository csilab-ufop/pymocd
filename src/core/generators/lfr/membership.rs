//! Who belongs where, and how many of each node's edges stay inside: the stage
//! that decides whether a draw of degrees and community sizes can be wired
//! exactly.
//!
//! Each node `v` of degree `k` gets the internal degree `(1 - mu) k`, rounded
//! stochastically (down, or up with probability equal to the fraction) as the
//! original code does (benchm.cpp, `int_interno`), so the expected internal
//! degree is `(1 - mu) k`. networkx's `round` is round-half-to-even:
//! deterministic, and biased whenever many nodes share a degree. No later step
//! may move internal degrees in a direction set by the state it finds, or the
//! expectation drifts: rounding the odd member of every odd community down,
//! say, takes half an edge end per community off the internal side.
//!
//! A node may join a community larger than the ceiling of `(1 - mu) k`, so
//! that either rounding fits. networkx admits on `round((1 - mu) k)` instead,
//! so when no placement exists on the ceilings, the fewest nodes that make one
//! exist are admitted on their floor, chosen among those whose `(1 - mu) k`
//! networkx rounds down: whenever networkx can place a draw, it is placed here
//! too. One that rounded up and lands where its ceiling does not fit is
//! rounded down, and another node that rounded down is rounded up in its
//! place: a member of the same community when one can be, else a member of
//! another community of the same size, else any node of the communities being
//! placed (see `repay`). That keeps the expected internal total, unless no
//! node can be rounded up; the shortfall is then kept and reported (see
//! `unpaid`). Admitting any node on its floor would also place draws networkx
//! cannot, such as `mu = 0.001` with communities of 3 nodes of degree 3, at an
//! internal total that rounding others up cannot restore.
//!
//! The placement is drawn uniformly from all placements in which every node
//! fits (see `place`). The homeless-node procedure of the LFR paper, which
//! networkx follows, walks towards such a placement by random evictions, and
//! when the seats that fit some nodes are exactly as many as those nodes, it
//! needs a number of steps that grows faster than `n`.
//!
//! A placement is wirable only if each community's internal degrees form a
//! graphical sequence (even sum, Erdős–Gallai) and the external stubs can be
//! paired across communities. networkx checks neither. Here:
//! - parity: when the communities being repaired have an odd internal total,
//!   one node's rounding is flipped, up or down by a fair coin (the original
//!   code's coin, binary_benchm.cpp); the odd communities, now even in number,
//!   are paired, and each pair exchanges a member of odd internal degree for
//!   one of even internal degree, which changes no node's internal degree. A
//!   pair no exchange fixes rounds one member up and one down. A community
//!   that is still odd and has no member whose rounding can change (every
//!   `(1 - mu) k` an integer, `0 < mu < 1`) moves one stub across, which
//!   leaves that node one edge off its target;
//! - a non-graphical community exchanges members with others, each exchange
//!   lowering its Erdős–Gallai violation and leaving the partner graphical;
//! - a community holding more than half of all external stubs exchanges its
//!   members of highest external degree for outsiders of lower external
//!   degree;
//! - what still fails is re-placed together with a few other communities,
//!   round after round, until nothing fails, the rounds stop making progress
//!   (the original code's `mate_trooper`) or the draw's repair work runs out
//!   (`REPAIR_WORK_PER_STUB`); the draw then fails, and the caller draws
//!   again.
//!
//! Which members an exchange moves depends on their internal degrees, so on
//! how they were rounded: a non-graphical community of 10 hands its
//! rounded-up members to whoever takes them, and a community holding hubs
//! takes in members of high internal degree, which rounded up more often
//! than not. Every exchange therefore looks for its partner among the
//! communities of the same size first, which moves rounding between
//! communities but not between sizes. What exchanges and re-placements
//! between sizes move is booked per community and given back at the end by
//! exchanging members of equal internal degree and different degree, which
//! changes no community's internal degrees (see `settle`). Without either,
//! small communities ended with fewer internal edges than `(1 - mu) k` over
//! their members and large ones with more.
//!
//! Rejecting a whole placement instead of repairing it does not scale: the
//! chance that every community is graphical at once falls exponentially with
//! their number.
//! This Source Code Form is subject to the terms of The GNU General Public License v3.0
//! Copyright 2026 - Guilherme Santos. If a copy of the MPL was not distributed with this
//! file, You can obtain one at https://www.gnu.org/licenses/gpl-3.0.html

use std::cmp::Reverse;
use std::collections::hash_map::Entry;

use rand::seq::SliceRandom;
use rand::{Rng, RngExt};
use rustc_hash::FxHashMap;

use super::graphical::{self, Histogram, Scratch, zverovich};
use super::{ASSIGN, LfrError, WIRE};

const UNPLACED: u32 = u32::MAX;
/// Most repair rounds for one draw. Each costs time linear in `n` besides
/// the repairs themselves, which `REPAIR_WORK_PER_STUB` bounds.
const REPAIR_ROUNDS: usize = 64;
/// Repair rounds without fewer failing communities before the draw fails.
const STALL_ROUNDS: usize = 8;
/// Exchange proposals per member for a non-graphical community, per partner
/// window.
const TRIES_PER_MEMBER: usize = 200;
/// The repair work one draw may spend, per edge end of the graph, plus
/// `REPAIR_WORK_BASE`: steps of the Erdős–Gallai evaluations, partner scans
/// and proposals that the graphical and the external-stub repairs make. A
/// draw that needs more fails, and the caller draws again, so a request that
/// no draw can meet costs time linear in the size of the graph per draw.
/// The successful draws measured used far less (see the user guide).
const REPAIR_WORK_PER_STUB: u64 = 5_000;
const REPAIR_WORK_BASE: u64 = 10_000_000;
/// Draws' worth of repair work a request may spend in all (see
/// `repair_pool`). A request whose draws keep failing in the repairs stops
/// after about this many such draws, even short of `max_iters`; draws that
/// fail before the repairs (placement, sizes) cost no repair work.
const POOL_DRAWS: u64 = 32;
/// Odd communities an odd community tries to pair with by an exchange. The
/// first nearly always works; the bound keeps a request where none can (every
/// member of odd internal degree) linear in the number of communities.
const PARTNER_TRIES: usize = 16;
/// Outsiders `balance_external` considers for each member it tries to move,
/// not counting those it has already moved.
const BALANCE_SCAN: usize = 4096;
/// Members in a row that `balance_external` scans for without finding an
/// exchange before it gives up: past the first members that can move, a
/// balance that is going to succeed rarely misses twice in a row, and one
/// that cannot (its hubs fit no other community, say) otherwise scans every
/// member it holds.
const BALANCE_MISSES: usize = 1024;
/// Draws in a row `settle` makes without finding an exchange before it stops.
const SETTLE_TRIES: usize = 16;
/// Relative distance within which `(1 - mu) k` is taken to be an integer.
const SNAP: f64 = 1e-9;
/// Absolute distance within which two `move_stub` costs count as a tie.
const COST_TIE: f64 = 1e-12;

/// A complete placement with internal degrees that can be wired exactly.
pub struct Assignment {
    /// Degree of each node, as drawn.
    pub deg: Vec<u32>,
    /// Internal degree of each node.
    pub kin: Vec<u32>,
    /// Community of each node, an index into `members`.
    pub comm: Vec<u32>,
    /// Nodes of each community.
    pub members: Vec<Vec<u32>>,
    sizes: Vec<u32>,
    /// Index of each node in `members[comm[v]]`.
    pos: Vec<u32>,
    /// Internal degree as first rounded; a re-placed node starts from it again.
    kin0: Vec<u32>,
    floor: Vec<u32>,
    /// `ceil((1 - mu) k)`, the largest internal degree rounding can give.
    cap: Vec<u32>,
    /// The internal degree a community must be larger than to admit `v`:
    /// `cap`, or `floor` for the nodes `make_room` admits on their floor.
    req: Vec<u32>,
    /// Whether networkx's `round((1 - mu) k)`, in its own floating point and
    /// round-half-to-even, is the floor of a fractional `(1 - mu) k`: the
    /// nodes `make_room` may admit on their floor.
    nx_floor: Vec<bool>,
    /// Fractional part of `(1 - mu) k`: the probability of rounding up.
    frac: Vec<f64>,
    /// Communities by increasing size, ties by index: exchanges look for
    /// partners near a community in this order.
    by_size: Vec<u32>,
    /// Rounding residual, `sum (k_in - (1 - mu) k)`, that exchanges and
    /// re-placements between communities of different sizes moved into each
    /// community; `settle` takes it back out in expectation.
    debt: Vec<f64>,
    /// Internal edge ends each node holds above its first rounding because
    /// `repay` rounded it up (one per round-up), less the one `join` took when
    /// it rounded the node down; reset when the node is re-placed. The fair
    /// parity moves are left out, so minus the sum is what `repay` could not
    /// give back (see `unpaid`).
    credit: Vec<i32>,
    /// Repair work left to this draw (see `REPAIR_WORK_PER_STUB`).
    work: u64,
    scratch: Scratch,
}

/// The repair work a whole request may spend over all its draws: as much as
/// `POOL_DRAWS` draws of `expected_stubs` edge ends may each spend.
pub fn repair_pool(expected_stubs: f64) -> u64 {
    let per_draw = REPAIR_WORK_BASE as f64 + REPAIR_WORK_PER_STUB as f64 * expected_stubs;
    (POOL_DRAWS as f64 * per_draw).min(u64::MAX as f64 / 2.0) as u64
}

/// Places every node and makes every community wirable, or fails with the
/// error of the stage that could not. The repairs spend at most the draw's
/// own budget and what is left of `pool`, from which their work is taken.
pub fn assign<R: Rng + ?Sized>(
    deg: Vec<u32>,
    sizes: &[usize],
    mu: f64,
    pool: &mut u64,
    rng: &mut R,
) -> Result<Assignment, LfrError> {
    let mut a = Assignment::new(deg, sizes, mu, rng);
    a.work = a.work.min(*pool);
    let budget = a.work;
    let outcome = a.complete(mu, rng);
    *pool -= budget - a.work;
    outcome.map(|()| a)
}

impl Assignment {
    /// `assign` on a fresh assignment.
    fn complete<R: Rng + ?Sized>(&mut self, mu: f64, rng: &mut R) -> Result<(), LfrError> {
        if !self.make_room(rng) {
            return Err(LfrError::runtime(ASSIGN));
        }
        let n = self.deg.len();
        let q = self.sizes.len();
        let all: Vec<u32> = (0..q as u32).collect();
        self.place((0..n as u32).collect(), &all, rng)?;

        let mut scope = all;
        let mut fewest = usize::MAX;
        let mut stalled = 0;
        let mut unrealizable = false;
        for _ in 0..REPAIR_ROUNDS {
            self.repay(&scope, rng);
            let mut bad = vec![false; q];
            let odd = self.fix_parity(&scope, mu, rng);
            unrealizable = !odd.is_empty();
            for c in odd {
                bad[c as usize] = true;
            }
            // a community's size is fixed, so re-placing does not help one that
            // keeps more than half of the external stubs after the exchanges; it
            // is checked before the costlier graphical repair
            if !self.balance_external(rng) {
                return Err(LfrError::runtime(WIRE));
            }
            // every failing community is known before any repair, so that none
            // is offered as a partner: an exchange must leave the partner
            // graphical, which almost never happens to one that is not
            let mut nongraphical = vec![false; q];
            for &c in &scope {
                let members = &self.members[c as usize];
                let kin = &self.kin;
                nongraphical[c as usize] =
                    graphical::deficit(members.iter().map(|&v| kin[v as usize]), &mut self.scratch)
                        > 0;
            }
            for &c in &scope {
                if bad[c as usize] || !nongraphical[c as usize] {
                    continue;
                }
                if self.make_graphical(c, &nongraphical, rng) {
                    nongraphical[c as usize] = false;
                } else {
                    bad[c as usize] = true;
                    unrealizable = true;
                }
            }
            self.mark_unwirable(&mut bad);
            let failing: Vec<u32> = (0..q as u32).filter(|&c| bad[c as usize]).collect();
            if failing.is_empty() {
                self.settle(rng);
                let mut bad = vec![false; q];
                self.mark_unwirable(&mut bad);
                if bad.contains(&true) || !self.external_graphical() {
                    return Err(LfrError::runtime(WIRE));
                }
                return Ok(());
            }
            if failing.len() < fewest {
                fewest = failing.len();
                stalled = 0;
            } else {
                stalled += 1;
                if stalled >= STALL_ROUNDS {
                    break;
                }
            }
            if self.work == 0 {
                break;
            }
            scope = self.reshuffle(&failing, rng)?;
        }
        // a community that cannot even hold its internal edges calls for larger
        // communities; one that holds them but cannot send its external stubs
        // anywhere calls for smaller ones
        Err(LfrError::runtime(if unrealizable { ASSIGN } else { WIRE }))
    }
}

/// Free seats per community, in a Fenwick tree over an order of the
/// communities, for drawing a seat uniformly within a prefix of that order.
struct Seats {
    tree: Vec<u64>,
}

impl Seats {
    fn new(free: &[u64]) -> Self {
        let mut tree = vec![0; free.len() + 1];
        for (i, &w) in free.iter().enumerate() {
            let i = i + 1;
            tree[i] += w;
            let parent = i + (i & i.wrapping_neg());
            if parent < tree.len() {
                tree[parent] += tree[i];
            }
        }
        Self { tree }
    }

    /// Free seats in the first `len` communities.
    fn prefix(&self, mut len: usize) -> u64 {
        let mut sum = 0;
        while len > 0 {
            sum += self.tree[len];
            len &= len - 1;
        }
        sum
    }

    /// The community holding the `u`-th free seat, counting from zero; `u` is
    /// below the total.
    fn find(&self, mut u: u64) -> usize {
        let mut at = 0;
        let mut step = (self.tree.len() - 1).next_power_of_two();
        while step > 0 {
            let next = at + step;
            if next < self.tree.len() && self.tree[next] <= u {
                at = next;
                u -= self.tree[next];
            }
            step >>= 1;
        }
        at
    }

    fn take(&mut self, i: usize) {
        let mut i = i + 1;
        while i < self.tree.len() {
            self.tree[i] -= 1;
            i += i & i.wrapping_neg();
        }
    }
}

/// What the Zverovich–Zverovich test needs of a community's internal
/// degrees: how many are positive, a lower bound on the smallest positive one
/// and an upper bound on the largest. Removing a degree leaves the bounds
/// valid, so they are kept through exchanges without rescanning.
#[derive(Clone, Copy)]
struct Spread {
    positive: u64,
    low: u32,
    high: u32,
}

impl Spread {
    const EMPTY: Self = Self {
        positive: 0,
        low: u32::MAX,
        high: 0,
    };

    fn add(&mut self, d: u32) {
        if d > 0 {
            self.positive += 1;
            self.low = self.low.min(d);
            self.high = self.high.max(d);
        }
    }

    /// The spread with one degree `out` replaced by `inn`.
    fn swapped(mut self, out: u32, inn: u32) -> Self {
        if out > 0 {
            self.positive -= 1;
        }
        self.add(inn);
        self
    }

    /// Graphical by the sufficient condition, for an even sum.
    fn graphical(self) -> bool {
        zverovich(self.positive, self.low, self.high)
    }
}

impl Assignment {
    fn new<R: Rng + ?Sized>(deg: Vec<u32>, sizes: &[usize], mu: f64, rng: &mut R) -> Self {
        let n = deg.len();
        let stubs: u64 = deg.iter().map(|&k| u64::from(k)).sum();
        let mut kin0 = Vec::with_capacity(n);
        let mut floor = Vec::with_capacity(n);
        let mut cap = Vec::with_capacity(n);
        let mut fracs = Vec::with_capacity(n);
        let mut nx_floor = Vec::with_capacity(n);
        for &k in &deg {
            let x = internal_target(k, mu);
            let f = x.floor();
            let frac = x - f;
            // one draw per node at every mu, so later stages see the same
            // stream position whatever mu is
            let up = rng.random::<f64>() < frac;
            nx_floor.push(frac > 0.0 && (f64::from(k) * (1.0 - mu)).round_ties_even() == f);
            let f = f as u32;
            floor.push(f);
            cap.push(f + u32::from(frac > 0.0));
            kin0.push(f + u32::from(up));
            fracs.push(frac);
        }
        let mut by_size: Vec<u32> = (0..sizes.len() as u32).collect();
        by_size.sort_unstable_by_key(|&c| (sizes[c as usize], c));
        Self {
            kin: kin0.clone(),
            deg,
            comm: vec![UNPLACED; n],
            members: sizes.iter().map(|&s| Vec::with_capacity(s)).collect(),
            sizes: sizes.iter().map(|&s| s as u32).collect(),
            pos: vec![0; n],
            kin0,
            floor,
            req: cap.clone(),
            cap,
            nx_floor,
            frac: fracs,
            by_size,
            debt: vec![0.0; sizes.len()],
            credit: vec![0; n],
            work: REPAIR_WORK_BASE + REPAIR_WORK_PER_STUB * stubs,
            scratch: Scratch::default(),
        }
    }

    /// Takes `steps` from the draw's repair work; false once it is spent.
    fn spend(&mut self, steps: u64) -> bool {
        self.work = self.work.saturating_sub(steps);
        self.work > 0
    }

    /// Whether `v` can move to `c` keeping its internal degree.
    fn fits(&self, v: u32, c: u32) -> bool {
        let v = v as usize;
        self.req[v].max(self.kin[v]) < self.sizes[c as usize]
    }

    /// The range of `by_size` holding the communities of `c`'s size.
    fn size_group(&self, c: u32) -> (usize, usize) {
        let s = self.sizes[c as usize];
        let lo = self
            .by_size
            .partition_point(|&d| self.sizes[d as usize] < s);
        let hi = self
            .by_size
            .partition_point(|&d| self.sizes[d as usize] <= s);
        (lo, hi)
    }

    /// Hall's condition for the placement, admitting as few nodes as possible
    /// on their floor. The communities that admit `v` are those larger than
    /// `req[v]`, nested by size, so a placement exists exactly when, for every
    /// `t`, the nodes with `req >= t` fit in the communities of size above
    /// `t`. Starting from `req = cap`, the levels are checked from the top; a
    /// level short of `d` seats lowers `req` to the floor for `d` nodes with
    /// `cap = t` and `nx_floor`, which changes that level's count and no
    /// other. With every such node lowered, `req` is networkx's `round`, so
    /// this fails only where networkx's placement does.
    ///
    /// The nodes lowered are drawn without regard to how they were rounded,
    /// so the placement stays independent of the rounding. Lowering those
    /// that rounded down first would spare some roundings, but would leave
    /// the ones that rounded up as the nodes that need a ceiling, and only
    /// the largest communities hold those. One lowered after rounding up
    /// that lands where its ceiling does not fit is rounded down by `join`
    /// and repaid by `repay` when some node can be rounded up.
    fn make_room<R: Rng + ?Sized>(&mut self, rng: &mut R) -> bool {
        let max_size = self.sizes.iter().copied().max().unwrap_or(0) as usize;
        let max_cap = self.cap.iter().copied().max().unwrap_or(0) as usize;
        let top = max_size.max(max_cap);
        let mut seats_of_size = vec![0u64; top + 2];
        for &s in &self.sizes {
            seats_of_size[s as usize] += u64::from(s);
        }
        let mut at_level = vec![0u64; top + 1];
        let mut lowerable_at = vec![0u64; top + 1];
        for (&c, &lowerable) in self.cap.iter().zip(&self.nx_floor) {
            at_level[c as usize] += 1;
            lowerable_at[c as usize] += u64::from(lowerable);
        }
        let mut lower = vec![0u64; top + 1];
        let (mut need, mut seats) = (0u64, 0u64);
        for t in (1..=top).rev() {
            seats += seats_of_size[t + 1];
            need += at_level[t];
            if need > seats {
                let short = need - seats;
                if short > lowerable_at[t] {
                    return false;
                }
                lower[t] = short;
                need -= short;
                at_level[t - 1] += short;
            }
        }
        if lower.iter().all(|&d| d == 0) {
            return true;
        }
        let mut lowerable: Vec<Vec<u32>> = vec![Vec::new(); top + 1];
        for v in 0..self.deg.len() {
            let c = self.cap[v] as usize;
            if lower[c] > 0 && self.nx_floor[v] {
                lowerable[c].push(v as u32);
            }
        }
        for t in 1..=top {
            let take = lower[t] as usize;
            let (picked, _) = lowerable[t].partial_shuffle(rng, take);
            for &v in picked.iter() {
                self.req[v as usize] = self.floor[v as usize];
            }
        }
        true
    }

    /// Seats the nodes of `free` in the communities of `comms`, all empty on
    /// entry, uniformly among the placements in which every node sits in a
    /// community larger than its `req` and every community is full.
    ///
    /// Nodes go in order of decreasing `req`, each to a seat drawn uniformly
    /// among the free seats of the communities larger than its `req`. Every
    /// node seated before `v` has a `req` at least `v`'s and sits in a
    /// community larger than that, so the number of seats open to `v` is the
    /// same whatever the earlier draws were: each placement is drawn with the
    /// same probability, and Hall's condition (`make_room`) means no node
    /// finds every seat taken. networkx's homeless-node procedure has no such
    /// guarantee: when the seats that fit some nodes are exactly as many as
    /// those nodes, it needs a number of evictions that grows faster than `n`.
    fn place<R: Rng + ?Sized>(
        &mut self,
        mut free: Vec<u32>,
        comms: &[u32],
        rng: &mut R,
    ) -> Result<(), LfrError> {
        let mut order = comms.to_vec();
        order.sort_unstable_by_key(|&c| (Reverse(self.sizes[c as usize]), c));
        let open: Vec<u64> = order
            .iter()
            .map(|&c| u64::from(self.sizes[c as usize]))
            .collect();
        let mut seats = Seats::new(&open);
        free.sort_unstable_by_key(|&v| (Reverse(self.req[v as usize]), v));
        for v in free {
            let r = self.req[v as usize];
            let fit = order.partition_point(|&c| self.sizes[c as usize] > r);
            let total = seats.prefix(fit);
            // make_room checked Hall's condition, and a re-placed scope held
            // its nodes before; reported rather than panicking across the FFI
            if total == 0 {
                return Err(LfrError::runtime(ASSIGN));
            }
            let i = seats.find(rng.random_range(0..total));
            seats.take(i);
            self.join(v, order[i]);
        }
        Ok(())
    }

    /// Seats `v` in `c` with its first rounding, rounded down if `c` cannot
    /// hold its ceiling (only a node `make_room` admitted on its floor).
    fn join(&mut self, v: u32, c: u32) {
        let members = &mut self.members[c as usize];
        self.pos[v as usize] = members.len() as u32;
        members.push(v);
        self.comm[v as usize] = c;
        let v = v as usize;
        self.kin[v] = self.kin0[v].min(self.sizes[c as usize] - 1);
        // at most one below: `make_room` admits no node below its floor
        self.credit[v] = -((self.kin0[v] - self.kin[v]) as i32);
    }

    /// Exchanges two nodes of different communities; both sizes are kept.
    fn exchange(&mut self, a: u32, b: u32) {
        let (a, b) = (a as usize, b as usize);
        let (ca, cb) = (self.comm[a], self.comm[b]);
        let (pa, pb) = (self.pos[a], self.pos[b]);
        self.members[ca as usize][pa as usize] = b as u32;
        self.members[cb as usize][pb as usize] = a as u32;
        self.comm.swap(a, b);
        self.pos.swap(a, b);
    }

    /// `k_in - (1 - mu) k` of `v`: its rounding, less the expected rounding.
    fn residual(&self, v: u32) -> f64 {
        let v = v as usize;
        f64::from(self.kin[v]) - f64::from(self.floor[v]) - self.frac[v]
    }

    /// Exchanges two members for a repair. Which members a repair picks
    /// depends on their internal degrees, so on their rounding: a community
    /// that needs partners for its hubs takes in nodes of high internal
    /// degree, which rounded up more often than not, and a parity exchange
    /// requires a parity that only one of a node's two roundings has. Between
    /// communities of one size that moves residual within a size; between
    /// sizes the residual moved is booked in `debt`.
    fn trade(&mut self, a: u32, b: u32) {
        let (ca, cb) = (
            self.comm[a as usize] as usize,
            self.comm[b as usize] as usize,
        );
        if self.sizes[ca] != self.sizes[cb] {
            let moved = self.residual(b) - self.residual(a);
            self.debt[ca] += moved;
            self.debt[cb] -= moved;
        }
        self.exchange(a, b);
    }

    /// Takes the residual booked in `debt` back out of each community that
    /// took some in, in expectation, without changing any node's internal
    /// degree: a member `u` is exchanged for a member `w` of a community
    /// that gave residual away (picked with probability proportional to how
    /// much) with the same internal degree and a higher degree. Both
    /// communities keep their internal degree sequences, so their parity and
    /// graphicality, and `(1 - mu) (k_w - k_u)` of residual goes back. The
    /// last exchange that would overshoot is made with the probability that
    /// makes the expected amount the debt; a community for which
    /// `SETTLE_TRIES` draws in a row find no exchange keeps the rest.
    ///
    /// Without it, communities without a partner of their size (at small `n`,
    /// the largest, which hold the hubs) took in rounding: `tau1 = 2`,
    /// `average_degree = 5`, `max_degree = 50`, communities of 10 to 50,
    /// `n = 2000`, `mu = 0.05`: 0.017 more internal degree per node than
    /// `(1 - mu) k` in the communities of 50, `z = 15.7`. Settling by
    /// rounding members the other way instead leaves the sizes neutral but
    /// rounds the hubs down: 0.07 fewer internal edges per node of degree 46.
    fn settle<R: Rng + ?Sized>(&mut self, rng: &mut R) {
        let q = self.sizes.len();
        let mut stubs = vec![0u64; q];
        for v in 0..self.deg.len() {
            stubs[self.comm[v] as usize] += u64::from(self.deg[v] - self.kin[v]);
        }
        let total: u64 = stubs.iter().sum();
        for c in 0..q {
            let mut left = self.debt[c];
            let mut misses = 0;
            while left > 0.0 && misses < SETTLE_TRIES {
                let credit: Vec<f64> = self.debt.iter().map(|&d| (-d).max(0.0)).collect();
                if !credit.iter().any(|&w| w > 0.0) {
                    break;
                }
                let d = pick_weighted(&credit, rng);
                // members of d by internal degree, for the members of c in
                // random order: the first with a partner of its internal
                // degree and a higher degree that keeps c within half of the
                // external stubs
                let mut theirs: Vec<(u32, u32)> = self.members[d]
                    .iter()
                    .filter(|&&w| self.fits(w, c as u32))
                    .map(|&w| (self.kin[w as usize], w))
                    .collect();
                theirs.sort_unstable();
                let mut ours = self.members[c].clone();
                ours.shuffle(rng);
                let found = ours.iter().find_map(|&u| {
                    if !self.fits(u, d as u32) {
                        return None;
                    }
                    let (ku, du) = (self.kin[u as usize], self.deg[u as usize]);
                    let lo = theirs.partition_point(|&(k, _)| k < ku);
                    let hi = theirs.partition_point(|&(k, _)| k <= ku);
                    let options: Vec<u32> = theirs[lo..hi]
                        .iter()
                        .map(|&(_, w)| w)
                        .filter(|&w| {
                            let dw = self.deg[w as usize];
                            dw > du && 2 * (stubs[c] + u64::from(dw - du)) <= total
                        })
                        .collect();
                    (!options.is_empty()).then_some((u, options))
                });
                let Some((u, options)) = found else {
                    misses += 1;
                    continue;
                };
                let w = options[rng.random_range(0..options.len())];
                let back = self.residual(u) - self.residual(w);
                if back > left && rng.random::<f64>() * back >= left {
                    break;
                }
                let more = u64::from(self.deg[w as usize] - self.deg[u as usize]);
                stubs[c] += more;
                stubs[d] -= more;
                self.exchange(u, w);
                self.debt[d] += back;
                left -= back;
            }
        }
        self.debt.fill(0.0);
    }

    fn internal_sum(&self, c: u32) -> u64 {
        self.members[c as usize]
            .iter()
            .map(|&v| u64::from(self.kin[v as usize]))
            .sum()
    }

    /// For each member that `join` rounded down, rounds up one member of the
    /// same community that rounded down by chance, picked with `flip_weight`:
    /// rounding a node down where its ceiling does not fit then costs its
    /// community no expected internal stub, only moves one. What a community
    /// cannot repay itself (no member on a floor it can leave) is repaid by
    /// the communities of `scope` of its size, and what they cannot by all of
    /// `scope`. What `scope` cannot repay either (every node on a floor sits
    /// in a community too small for its ceiling: degree 2 in communities of 2
    /// at `mu = 0.3`, say) stays unpaid. The draw is kept, since the shortfall
    /// comes from small communities that the request lets every draw have,
    /// and `unpaid` reports it.
    fn repay<R: Rng + ?Sized>(&mut self, scope: &[u32], rng: &mut R) {
        let mut unpaid: Vec<(u32, usize)> = Vec::new();
        for &c in scope {
            let owed = self.members[c as usize]
                .iter()
                .filter(|&&v| self.kin[v as usize] < self.kin0[v as usize])
                .count();
            if owed > 0 {
                let left = self.round_up(&[c], owed, rng);
                if left > 0 {
                    unpaid.push((self.sizes[c as usize], left));
                }
            }
        }
        if unpaid.is_empty() {
            return;
        }
        unpaid.sort_unstable();
        let mut by_size = scope.to_vec();
        by_size.sort_unstable_by_key(|&c| (self.sizes[c as usize], c));
        let mut spill = 0;
        for same in unpaid.chunk_by(|x, y| x.0 == y.0) {
            let s = same[0].0;
            let lo = by_size.partition_point(|&c| self.sizes[c as usize] < s);
            let hi = by_size.partition_point(|&c| self.sizes[c as usize] <= s);
            let owed = same.iter().map(|&(_, left)| left).sum();
            spill += self.round_up(&by_size[lo..hi], owed, rng);
        }
        if spill > 0 {
            // what is left here is counted by `unpaid`, from `credit`
            self.round_up(scope, spill, rng);
        }
    }

    /// Internal edge ends that `join` took from nodes rounded down for lack
    /// of room and that `repay` could not give back to other nodes: the
    /// internal total falls short of the first roundings', whose expectation
    /// is `(1 - mu)` times the degree sum, by this much (the fair parity moves
    /// aside). Zero whenever every such node was repaid.
    pub fn unpaid(&self) -> u64 {
        let held: i64 = self.credit.iter().map(|&c| i64::from(c)).sum();
        held.min(0).unsigned_abs()
    }

    /// Rounds up `count` members of `comms`, drawn without replacement with
    /// weights `flip_weight(v, true)`, as `count` calls of `flip` would: the
    /// `count` largest keys `ln(u) / w` (Efraimidis and Spirakis), in one pass
    /// instead of one pass per member. Returns how many it could not round up.
    fn round_up<R: Rng + ?Sized>(&mut self, comms: &[u32], count: usize, rng: &mut R) -> usize {
        let mut keyed: Vec<(f64, u32)> = Vec::new();
        for &c in comms {
            for &v in &self.members[c as usize] {
                let w = self.flip_weight(v, true);
                if w > 0.0 {
                    let u: f64 = 1.0 - rng.random::<f64>();
                    keyed.push((u.ln() / w, v));
                }
            }
        }
        let take = count.min(keyed.len());
        if take < keyed.len() {
            keyed.select_nth_unstable_by(take, |x, y| y.0.total_cmp(&x.0));
        }
        for &(_, v) in &keyed[..take] {
            self.kin[v as usize] += 1;
            self.credit[v as usize] += 1;
        }
        count - take
    }

    /// Makes the internal-degree sum of every community of `scope` even and
    /// returns, sorted, those it could not fix (possible only at `mu = 0`).
    ///
    /// Only fair moves are used, so the expected internal degrees are kept:
    /// one flip in a direction drawn by a fair coin when the scope's total is
    /// odd, exchanges that change no internal degree, and, for a pair of odd
    /// communities no exchange fixes, one flip up and one down. A direction
    /// is forced only when no member can move the other way: near `mu = 0`,
    /// when no fractional member of the scope rounded down, which each of
    /// them does with probability about `mu k`. The odd communities are
    /// paired in order of size, so that the members an exchange moves stay
    /// among communities of one size wherever there are two odd ones.
    fn fix_parity<R: Rng + ?Sized>(&mut self, scope: &[u32], mu: f64, rng: &mut R) -> Vec<u32> {
        let total: u64 = scope.iter().map(|&c| self.internal_sum(c)).sum();
        if total % 2 == 1 {
            let up = rng.random_bool(0.5);
            if !self.flip(scope, up, rng) {
                self.flip(scope, !up, rng);
            }
        }
        let mut odd: Vec<u32> = scope
            .iter()
            .copied()
            .filter(|&c| self.internal_sum(c) % 2 == 1)
            .collect();
        // in order of size, at random within one size
        odd.shuffle(rng);
        odd.sort_by_key(|&c| self.sizes[c as usize]);
        let unpaired = self.pair_by_exchange(&odd, rng);
        let mut left = Vec::new();
        for pair in unpaired.chunks(2) {
            if let [x, y] = *pair {
                let (a, b) = if rng.random_bool(0.5) { (x, y) } else { (y, x) };
                if self.flip_apart(a, b, rng) || self.flip_apart(b, a, rng) {
                    continue;
                }
            }
            left.extend_from_slice(pair);
        }
        if 0.0 < mu && mu < 1.0 {
            for &c in &left {
                let up = rng.random_bool(0.5);
                if !self.flip(&[c], up, rng) && !self.flip(&[c], !up, rng) {
                    self.move_stub(c, mu, rng);
                }
            }
            left.clear();
        }
        left.sort_unstable();
        left
    }

    /// Pairs the odd communities of `odd`, sorted by size, by `swap_parity`:
    /// each tries up to `PARTNER_TRIES` of the unpaired ones after it, the
    /// nearest in size first. Returns those left unpaired, in the same order.
    fn pair_by_exchange<R: Rng + ?Sized>(&mut self, odd: &[u32], rng: &mut R) -> Vec<u32> {
        let mut paired = vec![false; odd.len()];
        let mut unpaired = Vec::new();
        for i in 0..odd.len() {
            if paired[i] {
                continue;
            }
            paired[i] = true;
            let mut tries = 0;
            let mut found = false;
            for j in i + 1..odd.len() {
                if tries == PARTNER_TRIES {
                    break;
                }
                if paired[j] {
                    continue;
                }
                tries += 1;
                if self.swap_parity(odd[i], odd[j], rng) {
                    paired[j] = true;
                    found = true;
                    break;
                }
            }
            if !found {
                unpaired.push(odd[i]);
            }
        }
        unpaired
    }

    /// Rounds one member of `x` up and one of `y` down, or neither when one of
    /// the two moves is impossible; the internal total is kept.
    fn flip_apart<R: Rng + ?Sized>(&mut self, x: u32, y: u32, rng: &mut R) -> bool {
        let can = |c: u32, up: bool| {
            self.members[c as usize]
                .iter()
                .any(|&v| self.flip_weight(v, up) > 0.0)
        };
        if !(can(x, true) && can(y, false)) {
            return false;
        }
        self.flip(&[x], true, rng);
        self.flip(&[y], false, rng);
        true
    }

    /// The weight with which a flip in direction `up` picks `v`: the
    /// probability that its rounding went the other way, its fraction `f` on
    /// the floor (going up) and `1 - f` on the ceiling (going down), and zero
    /// when `v` cannot move that way: an integer `(1 - mu) k`, already on that
    /// side, or a ceiling its community cannot hold.
    ///
    /// A member sits on its floor with probability `1 - f` and on its ceiling
    /// with probability `f`, so its expected weight is `f (1 - f)` in either
    /// direction; with the direction drawn by a fair coin, it is then moved up
    /// as often as down. Picking uniformly instead moves members with a large
    /// fraction down more often than up: they are rarely on their floor.
    fn flip_weight(&self, v: u32, up: bool) -> f64 {
        let v = v as usize;
        if self.cap[v] == self.floor[v] {
            0.0
        } else if up {
            let room = self.cap[v] < self.sizes[self.comm[v] as usize];
            if self.kin[v] == self.floor[v] && room {
                self.frac[v]
            } else {
                0.0
            }
        } else if self.kin[v] == self.cap[v] {
            1.0 - self.frac[v]
        } else {
            0.0
        }
    }

    /// Rounds one member of `comms` the other way in direction `up`, picked
    /// with `flip_weight`; false when no member can move that way.
    fn flip<R: Rng + ?Sized>(&mut self, comms: &[u32], up: bool, rng: &mut R) -> bool {
        let total: f64 = comms
            .iter()
            .flat_map(|&c| &self.members[c as usize])
            .map(|&v| self.flip_weight(v, up))
            .sum();
        if total.is_nan() || total <= 0.0 {
            return false;
        }
        let mut u = rng.random::<f64>() * total;
        // the last member with positive weight absorbs any rounding residue
        let mut pick = None;
        'scan: for &c in comms {
            for &v in &self.members[c as usize] {
                let w = self.flip_weight(v, up);
                if w > 0.0 {
                    pick = Some(v as usize);
                    if u < w {
                        break 'scan;
                    }
                    u -= w;
                }
            }
        }
        let Some(v) = pick else {
            return false;
        };
        if up {
            self.kin[v] += 1;
        } else {
            self.kin[v] -= 1;
        }
        true
    }

    /// Makes two odd communities even by exchanging a member of odd internal
    /// degree for a member of even internal degree. Which of the two gives the
    /// odd one is a fair coin, then each member is uniform among those that
    /// fit. Weighting the direction by the number of such pairs instead lets a
    /// large community take a small one's single odd member about four times
    /// in five, so nodes that rounded up drift into large communities.
    fn swap_parity<R: Rng + ?Sized>(&mut self, x: u32, y: u32, rng: &mut R) -> bool {
        let split = |own: u32, other: u32| {
            let (mut odd, mut even) = (Vec::new(), Vec::new());
            for &v in &self.members[own as usize] {
                if self.fits(v, other) {
                    if self.kin[v as usize] % 2 == 1 {
                        odd.push(v);
                    } else {
                        even.push(v);
                    }
                }
            }
            (odd, even)
        };
        let (x_odd, x_even) = split(x, y);
        let (y_odd, y_even) = split(y, x);
        let x_gives_odd = match (
            !x_odd.is_empty() && !y_even.is_empty(),
            !x_even.is_empty() && !y_odd.is_empty(),
        ) {
            (true, true) => rng.random_bool(0.5),
            (true, false) => true,
            (false, true) => false,
            (false, false) => return false,
        };
        let (from_x, from_y) = if x_gives_odd {
            (&x_odd, &y_even)
        } else {
            (&x_even, &y_odd)
        };
        let a = from_x[rng.random_range(0..from_x.len())];
        let b = from_y[rng.random_range(0..from_y.len())];
        self.trade(a, b);
        true
    }

    /// Moves one stub of one member of `c` between its internal and external
    /// sides, choosing the move that puts that node's own mixing least far
    /// from `mu` (ties at random). Used only when no member of `c` can be
    /// rounded the other way, so the node ends one edge off its target range.
    fn move_stub<R: Rng + ?Sized>(&mut self, c: u32, mu: f64, rng: &mut R) {
        let size = self.sizes[c as usize];
        let mut best: Option<(f64, usize, bool)> = None;
        let mut ties = 0u32;
        for &v in &self.members[c as usize] {
            let v = v as usize;
            let (k, i) = (self.deg[v], self.kin[v]);
            let kf = f64::from(k);
            let now = ((kf - f64::from(i)) / kf - mu).abs();
            for up in [false, true] {
                let ok = if up { i < k.min(size - 1) } else { i > 0 };
                if !ok {
                    continue;
                }
                let next = if up { i + 1 } else { i - 1 };
                let cost = ((kf - f64::from(next)) / kf - mu).abs() - now;
                // Costs equal in exact arithmetic can differ in the last bits
                // (at mu = 0.3 an integral (1 - mu) k makes "up" cheaper by
                // 5e-17), which would turn every tie into the same direction;
                // compare within COST_TIE so the RNG breaks real ties.
                match best {
                    Some((b, _, _)) if cost > b + COST_TIE => {}
                    Some((b, _, _)) if cost >= b - COST_TIE => {
                        ties += 1;
                        if rng.random_range(0..ties) == 0 {
                            best = Some((cost, v, up));
                        }
                    }
                    _ => {
                        best = Some((cost, v, up));
                        ties = 1;
                    }
                }
            }
        }
        // an odd sum has a member with a positive internal degree
        if let Some((_, v, up)) = best {
            if up {
                self.kin[v] += 1;
            } else {
                self.kin[v] -= 1;
            }
        }
    }

    /// Lowers the Erdős–Gallai violation of `c` to zero by exchanging a random
    /// member for a member of another community with the same internal-degree
    /// parity (both sums keep their parity), accepting an exchange only when
    /// it lowers the violation of `c` and leaves the partner graphical.
    ///
    /// Partners come from the communities of `c`'s size first, then from
    /// windows of the size order that widen around them, `TRIES_PER_MEMBER`
    /// proposals per member each; a community flagged in `nongraphical` is
    /// never a partner. False when the last window, every community, did not
    /// get there, or when the draw's repair work is spent.
    fn make_graphical<R: Rng + ?Sized>(
        &mut self,
        c: u32,
        nongraphical: &[bool],
        rng: &mut R,
    ) -> bool {
        let ci = c as usize;
        let mut own = Histogram::new(self.members[ci].iter().map(|&v| self.kin[v as usize]));
        let (mut now, steps) = own.deficit();
        if now == 0 {
            return true;
        }
        if !self.spend(steps) {
            return false;
        }
        // the partners' histograms, built on first use; only looked up, so
        // the hash order never matters
        let mut theirs: FxHashMap<u32, Histogram> = FxHashMap::default();
        let size = self.members[ci].len();
        let q = self.by_size.len();
        let (group_lo, group_hi) = self.size_group(c);
        let mut widen = 0;
        loop {
            let lo = group_lo.saturating_sub(widen);
            let hi = (group_hi + widen).min(q);
            let open = self.by_size[lo..hi]
                .iter()
                .any(|&d| d != c && !nongraphical[d as usize]);
            if open {
                for _ in 0..TRIES_PER_MEMBER * size {
                    if now == 0 || !self.spend(1) {
                        break;
                    }
                    let d = self.by_size[rng.random_range(lo..hi)];
                    if d == c || nongraphical[d as usize] {
                        continue;
                    }
                    let di = d as usize;
                    let a = self.members[ci][rng.random_range(0..size)];
                    let b = self.members[di][rng.random_range(0..self.members[di].len())];
                    let (ka, kb) = (self.kin[a as usize], self.kin[b as usize]);
                    if ka == kb || (ka ^ kb) & 1 == 1 || !self.fits(b, c) || !self.fits(a, d) {
                        continue;
                    }
                    let (next, steps) = own.deficit_if(ka, kb);
                    if !self.spend(steps) {
                        break;
                    }
                    if next >= now {
                        continue;
                    }
                    let (partner, built) = match theirs.entry(d) {
                        Entry::Occupied(slot) => (slot.into_mut(), 0),
                        Entry::Vacant(slot) => {
                            let members = &self.members[di];
                            let kin = &self.kin;
                            let h = Histogram::new(members.iter().map(|&v| kin[v as usize]));
                            (slot.insert(h), members.len() as u64)
                        }
                    };
                    let (after, steps) = partner.deficit_if(kb, ka);
                    if !self.spend(built + steps) {
                        break;
                    }
                    if after > 0 {
                        continue;
                    }
                    partner.replace(kb, ka);
                    own.replace(ka, kb);
                    self.trade(a, b);
                    now = next;
                }
            }
            if now == 0 || self.work == 0 || (lo == 0 && hi == q) {
                break;
            }
            widen = (widen * 4).max(2);
        }
        now == 0
    }

    /// The external stubs pair up across communities only if no community
    /// holds more than half of them. When one does (one community with most
    /// of the nodes, say), its members of highest external degree are
    /// exchanged, one at a time, for outsiders of lower external degree, the
    /// lowest first, while it still holds more than half. An exchange keeps
    /// both internal sums' parity (same parity of internal degree), fits both
    /// nodes, raises neither community's Erdős–Gallai violation, and leaves
    /// the partner at most half. False when the community still holds more
    /// than half, or when the draw's repair work is spent.
    ///
    /// The exchanges needed grow with `n` (a community holding 54% of the
    /// stubs of a graph of 100,000 nodes needs hundreds), so each must cost
    /// little: the Zverovich–Zverovich condition settles graphicality in
    /// constant time for a large sparse community, and only where it fails is
    /// the deficit evaluated, from a histogram. Each member that leaves looks
    /// at the next `BALANCE_SCAN` outsiders not yet moved in, and the balance
    /// stops after `BALANCE_MISSES` members in a row found none.
    fn balance_external<R: Rng + ?Sized>(&mut self, rng: &mut R) -> bool {
        let n = self.deg.len();
        let q = self.sizes.len();
        let mut stubs = vec![0u64; q];
        for v in 0..n {
            stubs[self.comm[v] as usize] += u64::from(self.deg[v] - self.kin[v]);
        }
        let total: u64 = stubs.iter().sum();
        let Some(big) = (0..q).find(|&c| 2 * stubs[c] > total) else {
            return true;
        };
        let ext = |s: &Self, v: u32| s.deg[v as usize] - s.kin[v as usize];
        // no exchange can take it below the least any set of its size holds
        let mut least: Vec<u32> = (0..n as u32).map(|v| ext(self, v)).collect();
        let size = self.members[big].len();
        if size < n {
            least.select_nth_unstable(size);
        }
        if 2 * least[..size].iter().map(|&e| u64::from(e)).sum::<u64>() > total {
            return false;
        }
        // highest external degree first, lowest first; at random within one
        let mut leave = self.members[big].clone();
        leave.shuffle(rng);
        leave.sort_by_key(|&v| Reverse(ext(self, v)));
        let mut enter: Vec<u32> = (0..n as u32)
            .filter(|&v| self.comm[v as usize] as usize != big)
            .collect();
        enter.shuffle(rng);
        enter.sort_by_key(|&v| ext(self, v));
        let mut moved_in = vec![false; enter.len()];
        let mut first = 0;

        let mut spread = vec![Spread::EMPTY; q];
        for v in 0..n {
            spread[self.comm[v] as usize].add(self.kin[v]);
        }
        // exact deficits, from histograms built on first need
        let mut hist: FxHashMap<usize, (Histogram, u64)> = FxHashMap::default();
        // the verdicts for the member being moved out, by partner community
        // and the incoming internal degree (`usize::MAX` for the dominant
        // community itself): they depend on nothing else
        let mut verdict: FxHashMap<(usize, u32), bool> = FxHashMap::default();
        let roomiest = (0..q)
            .filter(|&c| c != big)
            .map(|c| self.sizes[c])
            .max()
            .unwrap_or(0);
        // the most the exchanges could move: every member another community
        // can hold, for the outsider of least external degree
        let fewest = enter.first().map_or(0, |&b| ext(self, b));
        let movable: u64 = leave
            .iter()
            .filter(|&&a| self.req[a as usize].max(self.kin[a as usize]) < roomiest)
            .map(|&a| u64::from(ext(self, a).saturating_sub(fewest)))
            .sum();
        if 2 * movable < 2 * stubs[big] - total {
            return false;
        }
        let mut misses = 0;
        for &a in &leave {
            if 2 * stubs[big] <= total || self.work == 0 {
                break;
            }
            // a member no other community can hold is skipped without a scan
            let need = self.req[a as usize].max(self.kin[a as usize]);
            if need >= roomiest {
                continue;
            }
            if misses == BALANCE_MISSES {
                break;
            }
            misses += 1;
            verdict.clear();
            let ea = ext(self, a);
            // an exchange moves ea - eb stubs, which some partner must have
            // room for: the outsiders below ea - room are skipped
            let room = (0..q)
                .filter(|&d| d != big)
                .map(|d| (total - 2 * stubs[d].min(total / 2)) / 2)
                .max()
                .unwrap_or(0);
            let low = u64::from(ea).saturating_sub(room);
            let from = enter.partition_point(|&b| u64::from(ext(self, b)) < low);
            let mut seen = 0;
            for i in first.max(from)..enter.len() {
                if seen == BALANCE_SCAN || !self.spend(1) {
                    break;
                }
                if moved_in[i] {
                    continue;
                }
                seen += 1;
                let b = enter[i];
                let eb = ext(self, b);
                if eb >= ea {
                    break;
                }
                let d = self.comm[b as usize] as usize;
                let (ka, kb) = (self.kin[a as usize], self.kin[b as usize]);
                let moved = u64::from(ea - eb);
                if d == big
                    || (ka ^ kb) & 1 == 1
                    || !self.fits(a, d as u32)
                    || !self.fits(b, big as u32)
                    || 2 * (stubs[d] + moved) > total
                {
                    continue;
                }
                let big_after = spread[big].swapped(ka, kb);
                let d_after = spread[d].swapped(kb, ka);
                let big_ok = big_after.graphical()
                    || *verdict
                        .entry((usize::MAX, kb))
                        .or_insert_with(|| self.no_worse(&mut hist, big, ka, kb));
                if !big_ok {
                    continue;
                }
                let d_ok = d_after.graphical()
                    || *verdict
                        .entry((d, kb))
                        .or_insert_with(|| self.no_worse(&mut hist, d, kb, ka));
                if !d_ok {
                    continue;
                }
                for (c, out, inn) in [(big, ka, kb), (d, kb, ka)] {
                    if let Some((h, now)) = hist.get_mut(&c) {
                        h.replace(out, inn);
                        *now = h.deficit().0;
                    }
                }
                spread[big] = big_after;
                spread[d] = d_after;
                self.trade(a, b);
                stubs[big] -= moved;
                stubs[d] += moved;
                moved_in[i] = true;
                while first < enter.len() && moved_in[first] {
                    first += 1;
                }
                misses = 0;
                break;
            }
        }
        2 * stubs[big] <= total
    }

    /// Whether replacing one internal degree `out` of community `c` by `inn`
    /// leaves its Erdős–Gallai deficit no larger, from a histogram of `c`
    /// kept in `hist` with its current deficit.
    fn no_worse(
        &mut self,
        hist: &mut FxHashMap<usize, (Histogram, u64)>,
        c: usize,
        out: u32,
        inn: u32,
    ) -> bool {
        let ((h, now), built) = match hist.entry(c) {
            Entry::Occupied(slot) => (slot.into_mut(), 0),
            Entry::Vacant(slot) => {
                let members = &self.members[c];
                let kin = &self.kin;
                let mut h = Histogram::new(members.iter().map(|&v| kin[v as usize]));
                let (now, steps) = h.deficit();
                (slot.insert((h, now)), members.len() as u64 + steps)
            }
        };
        let now = *now;
        let (after, steps) = h.deficit_if(out, inn);
        self.spend(built + steps) && after <= now
    }

    /// Flags the communities whose external stubs cannot all be paired across
    /// communities, by two necessary conditions: no community holds more than
    /// half of all external stubs, and no node needs more external neighbours
    /// than there are nodes with external stubs outside its community.
    fn mark_unwirable(&self, bad: &mut [bool]) {
        let q = self.sizes.len();
        let mut stubs = vec![0u64; q];
        let mut senders = vec![0u64; q];
        let (mut all_stubs, mut all_senders) = (0u64, 0u64);
        for (v, (&k, &i)) in self.deg.iter().zip(&self.kin).enumerate() {
            let out = u64::from(k - i);
            if out > 0 {
                let c = self.comm[v] as usize;
                stubs[c] += out;
                senders[c] += 1;
                all_stubs += out;
                all_senders += 1;
            }
        }
        for c in 0..q {
            if 2 * stubs[c] > all_stubs {
                bad[c] = true;
            }
        }
        for (v, (&k, &i)) in self.deg.iter().zip(&self.kin).enumerate() {
            let c = self.comm[v] as usize;
            if u64::from(k - i) > all_senders - senders[c] {
                bad[c] = true;
            }
        }
    }

    /// Whether the external degrees, taken as one sequence over all nodes,
    /// pass the Erdős–Gallai test: necessary for a simple graph between
    /// communities, which adds the constraint that no edge lies inside one.
    /// Hubs of a heavy-tailed draw can need more distinct external
    /// neighbours than the other nodes' external stubs allow; no exchange of
    /// members changes that, and the configuration model would spend its
    /// whole budget finding it out.
    fn external_graphical(&mut self) -> bool {
        let ext = self.deg.iter().zip(&self.kin).map(|(&k, &i)| k - i);
        graphical::deficit(ext, &mut self.scratch) == 0
    }

    /// Re-places each failing community together with two others, drawn from
    /// its size group when it has two not yet chosen, else from the nearest
    /// sizes, every node starting again from its first rounding. Returns the
    /// communities re-placed, sorted.
    fn reshuffle<R: Rng + ?Sized>(
        &mut self,
        failing: &[u32],
        rng: &mut R,
    ) -> Result<Vec<u32>, LfrError> {
        let q = self.sizes.len();
        let mut chosen = vec![false; q];
        for &c in failing {
            chosen[c as usize] = true;
        }
        let mut scope = failing.to_vec();
        for &c in failing {
            let (group_lo, group_hi) = self.size_group(c);
            let mut widen = 0;
            let mut pool = loop {
                let lo = group_lo.saturating_sub(widen);
                let hi = (group_hi + widen).min(q);
                let pool: Vec<u32> = self.by_size[lo..hi]
                    .iter()
                    .copied()
                    .filter(|&d| !chosen[d as usize])
                    .collect();
                if pool.len() >= 2 || (lo == 0 && hi == q) {
                    break pool;
                }
                widen = (widen * 2).max(1);
            };
            let take = pool.len().min(2);
            let (picked, _) = pool.partial_shuffle(rng, take);
            for &d in picked.iter() {
                chosen[d as usize] = true;
                scope.push(d);
            }
        }
        scope.sort_unstable();
        let sizes_differ = scope
            .iter()
            .any(|&c| self.sizes[c as usize] != self.sizes[scope[0] as usize]);
        if sizes_differ {
            for &c in &scope {
                let held: f64 = self.members[c as usize]
                    .iter()
                    .map(|&v| self.residual(v))
                    .sum();
                self.debt[c as usize] -= held;
            }
        }

        let mut free = Vec::new();
        for &c in &scope {
            free.append(&mut self.members[c as usize]);
        }
        for &v in &free {
            self.kin[v as usize] = self.kin0[v as usize];
            self.comm[v as usize] = UNPLACED;
        }
        self.place(free, &scope, rng)?;
        if sizes_differ {
            // at the first rounding: what `join` rounds down, `repay` restores
            for &c in &scope {
                let held: f64 = self.members[c as usize]
                    .iter()
                    .map(|&v| {
                        let v = v as usize;
                        f64::from(self.kin0[v]) - f64::from(self.floor[v]) - self.frac[v]
                    })
                    .sum();
                self.debt[c as usize] += held;
            }
        }
        Ok(scope)
    }
}

/// An index drawn with probability proportional to `weights` (not all zero).
fn pick_weighted<R: Rng + ?Sized>(weights: &[f64], rng: &mut R) -> usize {
    let total: f64 = weights.iter().sum();
    let mut u = rng.random::<f64>() * total;
    // the last positive weight absorbs any rounding residue
    let mut pick = 0;
    for (i, &w) in weights.iter().enumerate() {
        if w > 0.0 {
            pick = i;
            if u < w {
                break;
            }
            u -= w;
        }
    }
    pick
}

/// `(1 - mu) k`, the product networkx rounds, snapped to the nearest integer
/// when within a relative `SNAP` of it. For a decimal `mu` the exact product is
/// often an integer that the float lands a rounding error away from (`k = 10`,
/// `mu = 0.7` gives `3.0000000000000004`); unsnapped, that integer target would
/// get a fractional part of `4e-16`, become a flip candidate, and end a whole
/// edge off. Snapping moves a genuine target by at most `SNAP` times itself.
fn internal_target(k: u32, mu: f64) -> f64 {
    let x = f64::from(k) * (1.0 - mu);
    let r = x.round();
    if (x - r).abs() <= SNAP * r.max(1.0) {
        r
    } else {
        x
    }
}

/// The smallest internal degree a community must exceed to admit a node of
/// degree `k`: the floor of `(1 - mu) k` when networkx rounds it down, else
/// its ceiling (see `make_room`). It does not decrease with `k`.
pub fn lowest_admission(k: u32, mu: f64) -> u32 {
    let x = internal_target(k, mu);
    let f = x.floor();
    if x > f && (f64::from(k) * (1.0 - mu)).round_ties_even() == f {
        f as u32
    } else {
        x.ceil() as u32
    }
}

#[cfg(test)]
mod tests {
    use super::super::powerlaw::{PowerLaw, community_sizes, degree_sequence};
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    /// A repair pool no test runs out of.
    const fn unlimited() -> u64 {
        u64::MAX
    }

    fn sample(n: usize, mu: f64, seed: u64) -> Result<Assignment, LfrError> {
        let mut rng = StdRng::seed_from_u64(seed);
        let deg = degree_sequence(&PowerLaw::new(2.5, 12, 50), n, 500, &mut rng)?;
        let sizes = community_sizes(1.5, 20, 100, n, 500, &mut rng)?;
        assign(deg, &sizes, mu, &mut unlimited(), &mut rng)
    }

    /// Checks a finished assignment at `mu = tenths / 10`, against the exact
    /// targets `(10 - tenths) k / 10` in integers.
    fn check(a: &Assignment, tenths: u32) {
        let mut s = Scratch::default();
        for (c, members) in a.members.iter().enumerate() {
            assert_eq!(members.len(), a.sizes[c] as usize, "sizes are kept");
            for (j, &v) in members.iter().enumerate() {
                assert_eq!(a.comm[v as usize] as usize, c);
                assert_eq!(a.pos[v as usize] as usize, j);
            }
            assert_eq!(a.internal_sum(c as u32) % 2, 0, "even internal sum");
            let seq = members.iter().map(|&v| a.kin[v as usize]);
            assert_eq!(graphical::deficit(seq, &mut s), 0, "graphical");
        }
        let mut off = 0;
        for v in 0..a.deg.len() {
            let num = (10 - tenths) * a.deg[v];
            let (lo, hi) = (num / 10, num.div_ceil(10));
            let i = a.kin[v];
            assert_eq!((a.floor[v], a.cap[v]), (lo, hi), "node {v}: snapped target");
            if lo < hi {
                assert!(lo <= i && i <= hi, "node {v}: {i} vs [{lo}, {hi}]");
            } else if i != lo {
                assert_eq!(i.abs_diff(lo), 1, "node {v}: a stub move is one edge");
                off += 1;
            }
            let size = a.sizes[a.comm[v] as usize];
            assert!(a.kin[v] < size && a.req[v] < size, "fits");
        }
        assert!(off <= a.sizes.len(), "at most one stub move per community");
    }

    /// Pooled z-score of the internal stubs against their expectation
    /// `sum (1 - mu) k` over many assignments: the rounding variance
    /// `sum f (1 - f)` plus one for the coin flip of each graph.
    fn internal_total_z(n: usize, mu: f64, sizes: (usize, usize), seeds: u64) -> f64 {
        let (mut dev, mut var) = (0.0, 0.0);
        for seed in 0..seeds {
            let mut rng = StdRng::seed_from_u64(seed);
            // degrees below the smallest community: every node fits anywhere
            let deg = degree_sequence(&PowerLaw::new(2.5, 3, 9), n, 500, &mut rng).unwrap();
            let sizes = community_sizes(1.5, sizes.0, sizes.1, n, 500, &mut rng).unwrap();
            let a = assign(deg, &sizes, mu, &mut unlimited(), &mut rng).unwrap();
            for v in 0..n {
                dev += f64::from(a.kin[v]) - f64::from(a.floor[v]) - a.frac[v];
                var += a.frac[v] * (1.0 - a.frac[v]);
            }
            var += 1.0;
        }
        dev / var.sqrt()
    }

    #[test]
    fn every_community_is_even_and_graphical() {
        for tenths in [0, 1, 3, 5, 7, 8, 9, 10] {
            for seed in 0..3 {
                let a = sample(1500, f64::from(tenths) / 10.0, seed).unwrap();
                check(&a, tenths);
            }
        }
    }

    #[test]
    fn float_noise_is_snapped() {
        // decimal products that are integers but not in floating point
        assert_eq!(10.0 * (1.0 - 0.7), 3.0000000000000004);
        assert_eq!(internal_target(10, 0.7), 3.0);
        assert_eq!(internal_target(15, 0.8), 3.0);
        assert_eq!(internal_target(30, 0.7), 9.0);
        // genuine fractions are kept
        assert_eq!(internal_target(7, 0.3), 7.0 * (1.0 - 0.3));
        assert_eq!(internal_target(3, 0.7), 3.0 * (1.0 - 0.7));
        // admission: 2.1 rounds down in networkx, 5.6 up, 3 is exact
        assert_eq!(lowest_admission(3, 0.3), 2);
        assert_eq!(lowest_admission(8, 0.3), 6);
        assert_eq!(lowest_admission(10, 0.7), 3);
    }

    #[test]
    fn flips_keep_the_rounding() {
        // at mu = 0.3 the 0.7 k are fractional unless 10 divides k, so parity
        // is fixed by flips and every internal degree stays on its floor or
        // ceiling; at mu = 0.8 integer targets (5 divides k) must not move
        for (tenths, seed) in [(3, 4), (8, 5), (7, 6)] {
            let a = sample(2000, f64::from(tenths) / 10.0, seed).unwrap();
            check(&a, tenths);
            for v in 0..a.deg.len() {
                let num = (10 - tenths) * a.deg[v];
                let i = a.kin[v];
                assert!(
                    i == num / 10 || i == num.div_ceil(10),
                    "node {v}: {i} vs {num}/10"
                );
            }
        }
    }

    #[test]
    fn flips_pick_by_the_other_rounding() {
        let mut rng = StdRng::seed_from_u64(3);
        // mu = 0.7: k = 3 has fraction 0.9, k = 7 has fraction 0.1, k = 10 none
        let fresh = |rng: &mut StdRng| {
            let mut a = Assignment::new(vec![3, 7, 10, 3], &[4], 0.7, rng);
            for v in 0..4 {
                a.join(v, 0);
            }
            a.kin = vec![a.floor[0], a.floor[1] + 1, 3, a.floor[3]];
            a
        };
        let a = fresh(&mut rng);
        // on the floor the weight is f going up, on the ceiling 1 - f going
        // down, zero the other way; never the integer
        assert!((a.flip_weight(0, true) - 0.9).abs() < 1e-12);
        assert_eq!(a.flip_weight(0, false), 0.0);
        assert!((a.flip_weight(1, false) - 0.9).abs() < 1e-12);
        assert_eq!(a.flip_weight(1, true), 0.0);
        assert_eq!(
            (a.flip_weight(2, true), a.flip_weight(2, false)),
            (0.0, 0.0)
        );
        // no room for a ceiling: a community of 1 cannot hold k_in = 1
        let mut tight = Assignment::new(vec![3], &[1], 0.7, &mut rng);
        tight.join(0, 0);
        assert_eq!((tight.kin[0], tight.flip_weight(0, true)), (0, 0.0));
        // two floor nodes of equal weight: each goes up half the time, by one
        // flip or by round_up
        let mut ups = [[0u32; 2]; 2];
        for _ in 0..20_000 {
            for (way, count) in ups.iter_mut().enumerate() {
                let mut b = fresh(&mut rng);
                if way == 0 {
                    assert!(b.flip(&[0], true, &mut rng));
                } else {
                    assert_eq!(b.round_up(&[0], 1, &mut rng), 0);
                }
                assert_eq!((b.kin[1], b.kin[2]), (a.kin[1], 3));
                count[usize::from(b.kin[3] > b.floor[3])] += 1;
            }
        }
        // within 5 sd (sd = 71)
        for count in ups {
            assert!(count[0].abs_diff(10_000) < 355, "{ups:?}");
        }
        // round_up: two of the two can go up, the third cannot
        let mut b = fresh(&mut rng);
        assert_eq!(b.round_up(&[0], 3, &mut rng), 1);
        assert_eq!(b.kin, vec![b.cap[0], b.cap[1], 3, b.cap[3]]);
    }

    #[test]
    fn parity_repair_keeps_the_expected_internal_total() {
        // near mu = 0 nearly every node sits on its ceiling and near mu = 1 on
        // its floor, so a repair that moves whichever member can move drifts
        // by about half a stub per odd community: z of -20 or less here
        for mu in [0.001, 0.01, 0.5, 0.99, 0.999] {
            for sizes in [(10, 30), (20, 100)] {
                let z = internal_total_z(500, mu, sizes, 100);
                assert!(z.abs() < 5.0, "mu {mu}, sizes {sizes:?}: z = {z}");
            }
        }
    }

    #[test]
    fn swaps_repair_odd_pairs_at_mu_zero() {
        // at mu = 0 nothing is fractional: parity must come from exchanges
        let a = sample(1500, 0.0, 6).unwrap();
        check(&a, 0);
        assert!(a.kin.iter().zip(&a.deg).all(|(i, k)| i == k));
    }

    #[test]
    fn graphical_repair_exchanges_members() {
        let mut rng = StdRng::seed_from_u64(8);
        // at mu = 0 the internal degree is the degree. Community 0 holds
        // 2, 2, 0: two nodes cannot both have degree 2 with one isolated
        // neighbour. Community 1 holds 2, 2, 2, 2, 0. Exchanging the 0 for a 2
        // leaves a triangle in each.
        let mut a = Assignment::new(vec![2, 2, 0, 2, 2, 2, 2, 0], &[3, 5], 0.0, &mut rng);
        for v in 0..8 {
            a.join(v, u32::from(v >= 3));
        }
        let mut s = Scratch::default();
        assert_eq!(graphical::deficit([2, 2, 0], &mut s), 2);
        assert!(a.make_graphical(0, &[true, false], &mut rng));
        check(&a, 0);
        assert_eq!(
            a.members[0].iter().map(|&v| a.kin[v as usize]).sum::<u32>(),
            6
        );
    }

    /// An assignment at `mu` with node `v` in community `at[v]`.
    fn seated(deg: Vec<u32>, sizes: &[usize], mu: f64, at: &[u32], seed: u64) -> Assignment {
        let mut rng = StdRng::seed_from_u64(seed);
        let mut a = Assignment::new(deg, sizes, mu, &mut rng);
        for (v, &c) in at.iter().enumerate() {
            a.join(v as u32, c);
        }
        a
    }

    #[test]
    fn graphical_repair_prefers_partners_of_the_same_size() {
        // mu = 0, so internal degree = degree. Community 0 (3 nodes) is
        // 2, 2, 0, not graphical; it is fixed by taking a 2 for its 0.
        // Community 1 (3 nodes, 2, 1, 1) and community 2 (5 nodes,
        // 2, 2, 2, 2, 0) can both give one and stay graphical; community 1,
        // of the same size, is always the one used.
        let deg = vec![2, 2, 0, 2, 1, 1, 2, 2, 2, 2, 0];
        let at = [0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 2];
        for seed in 0..20 {
            let mut a = seated(deg.clone(), &[3, 3, 5], 0.0, &at, seed);
            let other_size = a.members[2].clone();
            let mut rng = StdRng::seed_from_u64(seed + 100);
            assert!(a.make_graphical(0, &[true, false, false], &mut rng));
            check(&a, 0);
            assert_eq!(a.members[2], other_size, "seed {seed}");
        }
        // when community 1 is 2, 2, 2, giving a 2 would leave it 2, 2, 0:
        // the window widens to community 2
        let deg = vec![2, 2, 0, 2, 2, 2, 2, 2, 2, 2, 0];
        for seed in 0..20 {
            let mut a = seated(deg.clone(), &[3, 3, 5], 0.0, &at, seed);
            let same_size = a.members[1].clone();
            let mut rng = StdRng::seed_from_u64(seed + 100);
            assert!(a.make_graphical(0, &[true, false, false], &mut rng));
            check(&a, 0);
            assert_eq!(a.members[1], same_size, "seed {seed}");
        }
    }

    #[test]
    fn placement_is_uniform_over_placements_that_fit() {
        // mu = 0: nodes 0 and 1 (degree 2) need a community above 2, nodes 2
        // and 3 (degree 0) fit anywhere. Of communities of 3 and 1, the 1
        // holds node 2 or node 3: two placements, each half the time.
        let mut counts = [0u32; 2];
        let mut rng = StdRng::seed_from_u64(9);
        for _ in 0..20_000 {
            let mut a = Assignment::new(vec![2, 2, 0, 0], &[3, 1], 0.0, &mut rng);
            a.place(vec![0, 1, 2, 3], &[0, 1], &mut rng).unwrap();
            let single = a.members[1][0];
            assert!(single == 2 || single == 3);
            counts[single as usize - 2] += 1;
        }
        // within 5 sd (sd = 71)
        assert!(counts[0].abs_diff(10_000) < 355, "{counts:?}");
        // mu = 0.5: degree 1 (0.5) and degree 9 (4.5), half and half, with
        // 10,000 communities of 1 and 200 of 50. The singletons hold only
        // nodes admitted on their floor of 0, and make_room admits exactly
        // as many as there are singletons: the tight case in which the
        // homeless-node walk needs many evictions per node
        let n = 20_000;
        let deg: Vec<u32> = (0..n).map(|v| if v % 2 == 0 { 1 } else { 9 }).collect();
        let mut sizes = vec![1usize; n / 2];
        sizes.extend(std::iter::repeat_n(50, n / 100));
        let mut a = Assignment::new(deg, &sizes, 0.5, &mut rng);
        assert!(a.make_room(&mut rng));
        let all: Vec<u32> = (0..sizes.len() as u32).collect();
        a.place((0..n as u32).collect(), &all, &mut rng).unwrap();
        for (c, m) in a.members.iter().enumerate() {
            assert_eq!(m.len(), sizes[c]);
            assert!(m.iter().all(|&v| a.req[v as usize] < sizes[c] as u32));
        }
    }

    #[test]
    fn dominant_community_is_balanced() {
        // mu = 0.5. Community 0 (6 nodes) holds three nodes of degree 6 and
        // three of degree 2: 3 * 3 + 3 * 1 = 12 of the 20 external stubs.
        // Exchanging one degree-6 node for a degree-2 node of community 1
        // leaves it 10, half.
        let deg = vec![6, 6, 6, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2];
        let at = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2];
        let external = |a: &Assignment, c: usize| -> u32 {
            a.members[c]
                .iter()
                .map(|&v| a.deg[v as usize] - a.kin[v as usize])
                .sum()
        };
        for seed in 0..20 {
            let mut a = seated(deg.clone(), &[6, 4, 4], 0.5, &at, seed);
            assert_eq!((external(&a, 0), external(&a, 1)), (12, 4));
            let mut rng = StdRng::seed_from_u64(seed + 100);
            assert!(a.balance_external(&mut rng));
            assert_eq!(external(&a, 0), 10);
            check(&a, 5);
            let mut bad = vec![false; 3];
            a.mark_unwirable(&mut bad);
            assert_eq!(bad, vec![false; 3]);
        }
        // one community holds every stub, and no exchange changes that
        let mut one = seated(vec![2, 2, 2, 2], &[4], 0.5, &[0, 0, 0, 0], 1);
        assert!(!one.balance_external(&mut StdRng::seed_from_u64(1)));
    }

    #[test]
    fn settling_gives_rounding_back_without_changing_internal_degrees() {
        // mu = 0.05: degree 3 (2.85) rounded up and degree 4 (3.8) rounded
        // down both have internal degree 3, residuals +0.15 and -0.8.
        // Community 0 holds six of the first, community 1 six of the second,
        // and community 0 owes 1.9: two exchanges of 0.95 each pay it
        let deg = vec![3, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4];
        let at = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1];
        for seed in 0..20 {
            let mut a = seated(deg.clone(), &[6, 6], 0.05, &at, seed);
            a.kin = vec![3; 12];
            a.debt = vec![1.9, -1.9];
            a.settle(&mut StdRng::seed_from_u64(seed));
            assert_eq!(a.kin, vec![3; 12], "no internal degree changes");
            let fours = a.members[0]
                .iter()
                .filter(|&&v| a.deg[v as usize] == 4)
                .count();
            assert_eq!(fours, 2, "seed {seed}");
            assert_eq!(a.debt, vec![0.0, 0.0]);
        }
        // trades book residual only between communities of different sizes
        let mut a = seated(
            vec![3, 4, 3, 4, 4, 3],
            &[2, 2, 2],
            0.05,
            &[0, 0, 1, 1, 2, 2],
            1,
        );
        a.kin = vec![3; 6];
        a.trade(0, 3);
        assert_eq!(a.debt, vec![0.0; 3]);
        let mut a = seated(vec![3, 4, 3, 4, 4], &[2, 3], 0.05, &[0, 0, 1, 1, 1], 1);
        a.kin = vec![3; 5];
        a.trade(0, 3);
        assert!((a.debt[0] + 0.95).abs() < 1e-12 && (a.debt[1] - 0.95).abs() < 1e-12);
    }

    #[test]
    fn hall_condition() {
        let mut rng = StdRng::seed_from_u64(1);
        // three nodes need a community above 10; one of 11 seats them
        let mut a = Assignment::new(vec![10, 10, 10, 1, 1, 1, 1, 1], &[11, 4], 0.0, &mut rng);
        assert!(a.make_room(&mut rng));
        assert_eq!(a.req, a.cap, "nobody admitted on the floor");
        // the same nodes with only communities of 10: integer targets
        let mut a = Assignment::new(vec![10, 10, 1, 1], &[10, 10], 0.0, &mut rng);
        assert!(!a.make_room(&mut rng));
        // twelve nodes need a community above 6, which only the 11 offers
        let mut a = Assignment::new(vec![6; 12], &[11, 6, 6], 0.0, &mut rng);
        assert!(!a.make_room(&mut rng));
        // mu = 0.3, degree 8: 5.6, which networkx rounds up to 6, so no
        // community of 6 admits it, there or here
        let mut a = Assignment::new(vec![8; 12], &[6, 6], 0.3, &mut rng);
        assert!(!a.make_room(&mut rng));
        // degree 3: 2.1, which networkx rounds down to 2: a community of 3
        // holds it on its floor
        let mut a = Assignment::new(vec![3; 6], &[3, 3], 0.3, &mut rng);
        assert!(a.make_room(&mut rng));
        assert_eq!(a.req, vec![2; 6]);
    }

    #[test]
    fn floors_admit_what_ceilings_cannot() {
        let mut rng = StdRng::seed_from_u64(2);
        // n = 2, mu = 0.5, degree 1, two communities of one node: 0.5 rounds
        // to 1 and fits nowhere, 0 fits; networkx places it on round(0.5) = 0
        let a = assign(vec![1, 1], &[1, 1], 0.5, &mut unlimited(), &mut rng).unwrap();
        assert_eq!((a.kin.clone(), a.req.clone()), (vec![0, 0], vec![0, 0]));
        // mu = 0.5: seven nodes of degree 9 (4.5, ceiling 5) and four of
        // degree 1 (0.5), communities of 5 and 6. Only the 6 holds a 5, so one
        // node of degree 9 must be admitted on its floor of 4, whichever way
        // it rounded: each of the seven about as often
        let deg = vec![9, 9, 9, 9, 9, 9, 9, 1, 1, 1, 1];
        let mut chosen = [0u32; 7];
        for _ in 0..7000 {
            let mut a = Assignment::new(deg.clone(), &[5, 6], 0.5, &mut rng);
            a.kin0 = vec![5, 5, 5, 4, 5, 5, 5, 0, 0, 0, 0];
            a.kin.clone_from(&a.kin0);
            assert!(a.make_room(&mut rng));
            let floored: Vec<usize> = (0..11).filter(|&v| a.req[v] < a.cap[v]).collect();
            assert_eq!(floored.len(), 1);
            chosen[floored[0]] += 1;
        }
        // 1000 each, sd 30
        assert!(chosen.iter().all(|&c| c.abs_diff(1000) < 150), "{chosen:?}");
        // all seven rounded up: one of them is admitted on its floor, and if
        // it lands in the 5 it is rounded down there and repaid by a member of
        // the same community going up, so each community's internal total is
        // what was drawn
        for seed in 0..50 {
            let mut rng = StdRng::seed_from_u64(seed);
            let mut a = Assignment::new(deg.clone(), &[5, 6], 0.5, &mut rng);
            a.kin0 = vec![5, 5, 5, 5, 5, 5, 5, 0, 0, 0, 0];
            a.kin.clone_from(&a.kin0);
            assert!(a.make_room(&mut rng));
            assert_eq!((0..11).filter(|&v| a.req[v] < a.cap[v]).count(), 1);
            a.place((0..11).collect(), &[0, 1], &mut rng).unwrap();
            let drawn: Vec<u32> = (0..2)
                .map(|c| a.members[c].iter().map(|&v| a.kin0[v as usize]).sum())
                .collect();
            a.repay(&[0, 1], &mut rng);
            for (c, &want) in drawn.iter().enumerate() {
                let now: u32 = a.members[c].iter().map(|&v| a.kin[v as usize]).sum();
                assert_eq!(now, want, "seed {seed}, community {c}");
            }
            assert!((0..11).all(|v| a.kin[v] < a.sizes[a.comm[v] as usize]));
            assert_eq!(a.unpaid(), 0);
        }
    }

    #[test]
    fn unpaid_counts_what_no_node_can_take_back() {
        // degree 2 at mu = 0.3 (1.4, which networkx rounds to 1) in
        // communities of 2, which hold one internal edge per node: every node
        // is admitted on its floor, and none can be rounded up for another
        let mut rng = StdRng::seed_from_u64(4);
        let mut a = Assignment::new(vec![2; 6], &[2, 2, 2], 0.3, &mut rng);
        a.kin0 = vec![2, 1, 2, 2, 1, 1];
        a.kin.clone_from(&a.kin0);
        assert!(a.make_room(&mut rng));
        a.place((0..6).collect(), &[0, 1, 2], &mut rng).unwrap();
        a.repay(&[0, 1, 2], &mut rng);
        assert_eq!(a.kin, vec![1; 6]);
        assert_eq!(a.unpaid(), 3);
        // a re-placement starts the count of the nodes it moves again
        let scope = a.reshuffle(&[0, 1, 2], &mut rng).unwrap();
        a.repay(&scope, &mut rng);
        assert_eq!(a.unpaid(), 3);
        // a community of 3 holds the ceiling of 2, so four of seven nodes are
        // admitted on their floor and seated in the communities of 2. If
        // those four rounded up and the three others down, the three are
        // rounded up for four nodes: one edge end stays unpaid
        let mut a = Assignment::new(vec![2; 7], &[2, 2, 3], 0.3, &mut rng);
        assert!(a.make_room(&mut rng));
        a.kin0 = (0..7)
            .map(|v| if a.req[v] < a.cap[v] { 2 } else { 1 })
            .collect();
        assert_eq!(a.kin0.iter().filter(|&&k| k == 2).count(), 4);
        a.kin.clone_from(&a.kin0);
        a.place((0..7).collect(), &[0, 1, 2], &mut rng).unwrap();
        a.repay(&[0, 1, 2], &mut rng);
        assert_eq!(
            a.members[2].iter().map(|&v| a.kin[v as usize]).sum::<u32>(),
            6
        );
        assert_eq!(a.unpaid(), 1);
    }
}
