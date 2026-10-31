//! LFR benchmark graphs (Lancichinetti, Fortunato & Radicchi, Phys. Rev. E 78,
//! 046110, 2008) with crisp communities: undirected, unweighted, no overlap.
//!
//! The API, the defaults and the parameter checks are those of networkx 3.4.2's
//! `LFR_benchmark_graph`, and so are the degree law for a `min_degree` (a
//! discrete power law truncated to `[min_degree, max_degree]`, redrawn until
//! its sum is even) and the community-size law (a truncated power law
//! conditioned on summing to `n`). For an `average_degree` the degree law is
//! the power law with less mass at its lowest degree, so that its mean is the
//! request exactly (see `powerlaw::law_for_average`).
//! The edges are wired the way the original algorithm wires them, not the way
//! networkx does. networkx's wiring loop compares a node's whole degree, edges
//! other nodes already attached included, against its internal target, and
//! draws each new neighbour without checking whether that neighbour already
//! has its full degree, so its graphs overshoot both the degrees and the mixing
//! (networkx issue 6809). Here every node ends with exactly its drawn degree
//! `k`, of which `floor((1 - mu) k)` or `ceil((1 - mu) k)` edges stay inside its
//! community, with `(1 - mu) k` as the expected value of their total, so the
//! realized mixing is `mu` up to rounding noise. The one exception to the
//! rounding is a community whose internal degrees sum to an odd number with no
//! member whose rounding can change: one node ends one edge off (see
//! `membership`). The one exception to the expected total is a request whose
//! communities are too small for their members' ceilings, so that nodes that
//! rounded up are rounded down and no other node can be rounded up in their
//! place; `mixing_shortfall` warns when that moves the mixing by more than
//! `MIXING_SHIFT`.
//!
//! Not every draw of degrees and sizes can be wired that way: a hub may need a
//! community larger than any drawn, or one community may hold most of the
//! nodes and so most of the external stubs. Such a draw is replaced by a new
//! one, up to `max_iters` draws, the error of the last one raised after that.
//! The repairs that try to save a draw share one work budget per request
//! (`membership::repair_pool`), so a request whose draws all fail there stops
//! in time linear in the size of the graph. The graphs returned are
//! conditioned on their draw having been completed: where many draws fail,
//! their degrees and sizes lean away from the requested laws (fewer hubs,
//! for instance), and `shortfall` warns when the mean degree misses
//! `average_degree` by far; a `min_degree` request has no such warning. The
//! repairs are heuristics, so a rejected draw
//! may have had an exact wiring that they did not find. The original code
//! stops at the first failure instead (benchm.cpp, `benchmark`,
//! `return -1`).
//!
//! Pipeline, one stage per file:
//! - `params`: validation, in networkx's order and wording;
//! - `powerlaw`: the degree law for an `average_degree`, the degree sequence
//!   and the community sizes;
//! - `membership`: internal degrees, community assignment, and the repairs that
//!   make each community's internal degree sequence realizable;
//! - `graphical`: the Erdős–Gallai test used by those repairs;
//! - `wiring`: Havel–Hakimi plus degree-preserving swaps inside communities, a
//!   configuration model plus repair (or a greedy construction where that
//!   stalls) and swaps between them.
//!
//! Each draw has its own `StdRng`, seeded from the caller's seed, the degree
//! and size laws and the draw's number (see `stream_seed`), and drives every
//! stage of that draw in a fixed order; no hash container is ever iterated. A
//! seed therefore fixes the output for a given build of the crate, and calls
//! that differ only in `mu` (which no law depends on) draw the same degrees
//! and community sizes in their `i`-th draw.
//! This Source Code Form is subject to the terms of The GNU General Public License v3.0
//! Copyright 2026 - Guilherme Santos. If a copy of the MPL was not distributed with this
//! file, You can obtain one at https://www.gnu.org/licenses/gpl-3.0.html

mod graphical;
mod membership;
mod params;
mod powerlaw;
mod wiring;

use std::fmt;

use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};

use params::DegreeTarget;
pub use params::LfrParams;

/// networkx's message when no minimum degree gives the requested average.
pub const MATCH_AVERAGE: &str = "Could not match average_degree";
/// networkx's message when no degree or size sequence satisfies its condition.
pub const POWER_LAW: &str = "Could not create power law sequence";
/// networkx's message when the nodes cannot be fitted into the communities.
pub const ASSIGN: &str = "Could not assign communities; try increasing min_community";
/// No networkx counterpart: networkx either loops forever when a node finds no
/// new neighbour outside its community, or returns a graph whose degrees
/// overshoot.
pub const WIRE: &str = "Could not wire the inter-community edges; try decreasing max_community";

/// Why a request produced no graph. The binding maps `Value` to `ValueError`
/// (networkx's `NetworkXError`) and `Runtime` to `RuntimeError` (networkx's
/// `ExceededMaxIterations`).
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum LfrError {
    /// The parameters are invalid on their own.
    Value(String),
    /// The parameters are valid but no graph was found for them.
    Runtime(String),
}

impl LfrError {
    pub(crate) fn value(msg: &str) -> Self {
        Self::Value(msg.to_owned())
    }

    pub(crate) fn runtime(msg: &str) -> Self {
        Self::Runtime(msg.to_owned())
    }

    pub fn message(&self) -> &str {
        match self {
            Self::Value(m) | Self::Runtime(m) => m,
        }
    }
}

impl fmt::Display for LfrError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.message())
    }
}

impl std::error::Error for LfrError {}

/// A generated graph on nodes `0..n`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct LfrGraph {
    /// Every undirected edge once, as `(u, v)` with `u < v`, sorted.
    pub edges: Vec<(u32, u32)>,
    /// Community of each node, numbered `0..k` by first appearance in node order.
    pub communities: Vec<u32>,
    /// Draws made, the last one the graph's; the others were rejected.
    pub draws: usize,
    /// At most one message per cause, for the caller to raise as warnings:
    /// the mixing moved by internal edge ends that no node could take back
    /// (see `mixing_shortfall`), and a mean degree that misses
    /// `average_degree` by more than `SHORTFALL_SE` standard errors of the
    /// degree law and by more than `SHORTFALL_SHARE` of it (see `shortfall`).
    pub warnings: Vec<String>,
}

/// Standard errors of the mean degree, `sd / sqrt(n)` of the degree law, and
/// fraction of `average_degree`, by both of which a graph's mean degree must
/// miss the request before `generate` warns. Graphs from a single draw
/// crossed both in 12 of the 141,600 graphs measured for the user guide,
/// about one in 11,800; the share keeps the warning to misses large enough
/// to matter at large `n`.
const SHORTFALL_SE: f64 = 4.0;
const SHORTFALL_SHARE: f64 = 0.02;

/// The shift of the fraction of edges between communities, from internal
/// edge ends that no node could take back, above which `generate` warns; a
/// smaller shift is not reported.
const MIXING_SHIFT: f64 = 1e-3;

/// SplitMix64's output function: a bijective mix of a 64-bit word.
const fn splitmix64(mut x: u64) -> u64 {
    x = x.wrapping_add(0x9E37_79B9_7F4A_7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^ (x >> 31)
}

/// The seed of a request's first draw: the caller's seed mixed with every
/// parameter that shapes the degree and size laws, `mu` excluded. The degree
/// law is `(kmin, alpha)` (see `powerlaw::law_for_average`; `alpha = 0` for a
/// `min_degree`).
///
/// Each degree draw consumes one uniform whatever the law, so with the bare
/// seed two requests that differ only in `tau1` or the degree law reach the
/// size stage at the same point of the stream, and their community sizes come
/// out correlated. Mixing the law into the seed makes such requests draw
/// independently, while requests that differ only in `mu` (which is used after
/// both laws are sampled) keep the same degrees and community sizes.
fn stream_seed(seed: u64, req: &params::Request, kmin: usize, alpha: f64) -> u64 {
    let unset = u64::MAX;
    let words = [
        req.n as u64,
        req.tau1.to_bits(),
        kmin as u64,
        alpha.to_bits(),
        req.kmax as u64,
        req.tau2.to_bits(),
        req.min_community.map_or(unset, |c| c as u64),
        req.max_community.map_or(unset, |c| c as u64),
    ];
    words
        .iter()
        .fold(splitmix64(seed), |h, &w| splitmix64(h ^ splitmix64(w)))
}

/// The seed of draw `i`. A draw's stream does not depend on how much of its
/// own stream an earlier, failed draw used, so the `i`-th draws of two
/// requests that differ only in `mu` share their degrees and sizes.
const fn draw_seed(first: u64, i: usize) -> u64 {
    if i == 0 {
        first
    } else {
        splitmix64(first ^ splitmix64(i as u64))
    }
}

/// The error for `min_degree = n`: no node of a simple graph has `n`
/// neighbours. networkx 3.4.2 gets past its degree stage only for an even `n`
/// (the sum `n * n` is odd otherwise) and past its size stage as we do. It then
/// gives every node the internal degree `round((1 - mu) n)` and seats a node
/// only in a community larger than that, so a community no larger leaves seats
/// unfilled and it raises "Could not assign communities". Otherwise every node
/// is seated and its wiring loop never ends; the degree law has no value a
/// simple graph can realize, which is networkx's power-law error.
fn degree_n_error(req: &params::Request, sizes: &[usize]) -> LfrError {
    // networkx's own product and round-half-to-even
    let internal = (req.n as f64 * (1.0 - req.mu)).round_ties_even();
    if sizes.iter().any(|&s| s as f64 <= internal) {
        LfrError::runtime(ASSIGN)
    } else {
        LfrError::runtime(POWER_LAW)
    }
}

/// Generates one LFR graph. `seed = None` seeds from OS entropy.
///
/// Draws degrees, community sizes, an assignment and the edges, and draws
/// all of them again when the draw cannot be completed: up to `max_iters`
/// draws, then the error of the stage that failed in the last one. Failures
/// that every draw would repeat are raised at once: a degree law whose every
/// sequence has an odd sum, no size sequence when the size law does not
/// depend on the draw, `min_degree = n`, and a lowest degree that fits no
/// community any draw can have.
pub fn generate(params: &LfrParams, seed: Option<u64>) -> Result<LfrGraph, LfrError> {
    let req = params::validate(params)?;
    let (kmin, alpha) = match req.degree {
        DegreeTarget::Min(k) => (k, 0.0),
        DegreeTarget::Average(avg) => powerlaw::law_for_average(req.tau1, avg, req.kmax, req.tol)?,
    };
    let seed = seed.unwrap_or_else(|| rand::make_rng::<StdRng>().random());
    let first = stream_seed(seed, &req, kmin, alpha);
    // min_degree = n leaves kmin above kmax; its law is then the point [n, n]
    let degree_law = powerlaw::PowerLaw::mixture(req.tau1, kmin, req.kmax.max(kmin), alpha);
    if powerlaw::always_odd(&degree_law, req.n) {
        return Err(LfrError::runtime(POWER_LAW));
    }
    // the size law depends on the draw only through a default bound
    let fixed_sizes = req.min_community.is_some() && req.max_community.is_some();
    // the largest community any draw can have: the default max_community is
    // the largest drawn degree, at most kmax
    let largest = req.max_community.unwrap_or(req.kmax).min(req.n);
    let lowest = membership::lowest_admission(kmin as u32, req.mu) as usize;
    // the repairs of all draws together, so that a request whose draws keep
    // failing there stops in time linear in the size of the graph
    let mut pool = membership::repair_pool(req.n as f64 * degree_law.moments().0);

    let mut last = LfrError::runtime(POWER_LAW);
    for i in 0..req.max_iters.max(1) {
        if pool == 0 {
            break;
        }
        let mut rng = StdRng::seed_from_u64(draw_seed(first, i));
        // max_iters attempts at an even sum within the draw, as networkx makes
        let deg = match powerlaw::degree_sequence(&degree_law, req.n, req.max_iters, &mut rng) {
            Ok(deg) => deg,
            Err(e) => {
                last = e;
                continue;
            }
        };
        // networkx's defaults: the community sizes span the drawn degrees
        let lo = req
            .min_community
            .unwrap_or_else(|| deg.iter().copied().min().unwrap_or(1) as usize);
        let hi = req
            .max_community
            .unwrap_or_else(|| deg.iter().copied().max().unwrap_or(1) as usize);
        let sizes =
            match powerlaw::community_sizes(req.tau2, lo, hi, req.n, req.max_iters, &mut rng) {
                Ok(sizes) => sizes,
                Err(e) if fixed_sizes || req.min_community.is_some_and(|c| c > req.n) => {
                    return Err(e);
                }
                Err(e) => {
                    last = e;
                    continue;
                }
            };
        if kmin > req.kmax {
            return Err(degree_n_error(&req, &sizes));
        }
        // every node needs a community larger than `lowest`
        if lowest >= largest {
            return Err(LfrError::runtime(ASSIGN));
        }
        match build(deg, &sizes, req.mu, &mut pool, &mut rng) {
            Ok(mut graph) => {
                graph.draws = i + 1;
                if let DegreeTarget::Average(avg) = req.degree {
                    graph.warnings.extend(shortfall(&graph, avg, &degree_law));
                }
                return Ok(graph);
            }
            Err(e) => last = e,
        }
    }
    Err(last)
}

/// The warning for a graph whose mean degree misses `average_degree` by
/// more than `SHORTFALL_SE` standard errors of the degree law and more than
/// `SHORTFALL_SHARE` of it. A single draw did that in 12 of the 141,600
/// graphs measured for the user guide; a request most of whose draws are
/// rejected does it every time, since the draws that can be placed and wired
/// are not a fair sample of the degree law (fewer large hubs, typically).
///
/// A request that `tol` accepted outside the reachable range uses the law at
/// the nearer end of it, so the graph is compared with that law's mean: an
/// `average_degree` of `-1e299` with `tol = 1e300` is not "missed by" 100%.
fn shortfall(graph: &LfrGraph, avg: f64, law: &powerlaw::PowerLaw) -> Option<String> {
    let n = graph.communities.len() as f64;
    let mean = 2.0 * graph.edges.len() as f64 / n;
    let (law_mean, sd) = law.moments();
    // the law for a request inside the range has the request as its mean, up
    // to floating-point error
    let inside = (law_mean - avg).abs() <= 1e-6 * law_mean;
    let target = if inside { avg } else { law_mean };
    let miss = mean - target;
    let se = sd / n.sqrt();
    if miss.abs() <= SHORTFALL_SE * se || miss.abs() <= SHORTFALL_SHARE * target {
        return None;
    }
    let rejected = graph.draws - 1;
    let why = if rejected == 0 {
        "No draw was rejected, so this is an unusual draw of the degree law itself.".to_owned()
    } else {
        format!(
            "{rejected} of the {} draws of degrees and community sizes were rejected because \
             they could not be placed or wired, and the draws that can be are not a fair \
             sample of the degree law{}. A lower max_degree makes fewer draws fail, if \
             average_degree stays reachable.",
            graph.draws,
            if miss < 0.0 {
                ": draws with more high-degree nodes fail more often"
            } else {
                ""
            }
        )
    };
    let against = if inside {
        format!("average_degree={avg}")
    } else {
        format!(
            "{law_mean:.3}, the mean of the degree law used for average_degree={avg} (the end of \
             the reachable range within tol of it)"
        )
    };
    Some(format!(
        "lfr_graph: the mean degree of the returned graph is {mean:.3}, {:.1}% {} {against} \
         ({:.1} standard errors of the degree law). {why}",
        100.0 * miss.abs() / target,
        if miss < 0.0 { "below" } else { "above" },
        miss.abs() / se,
    ))
}

/// The warning for a graph whose fraction of edges between communities is
/// raised by more than `MIXING_SHIFT` by the `unpaid` internal edge ends of
/// its assignment (see `membership::Assignment::unpaid`): nodes that rounded
/// `(1 - mu) k` up sit in communities too small for that ceiling, and no node
/// elsewhere could be rounded up in their place. The draw is kept, since the
/// cause is the community sizes the request allows, so the message gives the
/// mixing of the graph returned and the remedy.
fn mixing_shortfall(
    unpaid: u64,
    stubs: u64,
    mu: f64,
    edges: &[(u32, u32)],
    comm: &[u32],
) -> Option<String> {
    if stubs == 0 || edges.is_empty() {
        return None;
    }
    let shift = unpaid as f64 / stubs as f64;
    if shift <= MIXING_SHIFT {
        return None;
    }
    let between = edges
        .iter()
        .filter(|&&(u, v)| comm[u as usize] != comm[v as usize])
        .count();
    let realized = between as f64 / edges.len() as f64;
    Some(format!(
        "lfr_graph: the fraction of edges between communities in the returned graph is \
         {realized:.4}, against mu={mu}. Some communities are too small to hold their members' \
         (1 - mu) * k internal edges: nodes that rounded (1 - mu) * k up were rounded down to \
         fit, and for {unpaid} of them no other node could be rounded up in their place, which \
         adds {shift:.4} to the expected fraction. A larger min_community (with max_community \
         at least as large) avoids this."
    ))
}

/// One draw's assignment and edges, with communities relabelled by first
/// appearance in node order.
fn build(
    deg: Vec<u32>,
    sizes: &[usize],
    mu: f64,
    pool: &mut u64,
    rng: &mut StdRng,
) -> Result<LfrGraph, LfrError> {
    let assignment = membership::assign(deg, sizes, mu, pool, rng)?;
    let mut edges = wiring::wire(&assignment, rng)?;
    edges.sort_unstable();
    let stubs = assignment.deg.iter().map(|&k| u64::from(k)).sum();
    let warnings = mixing_shortfall(assignment.unpaid(), stubs, mu, &edges, &assignment.comm)
        .into_iter()
        .collect();

    let mut relabel = vec![u32::MAX; sizes.len()];
    let mut next = 0;
    let communities = assignment
        .comm
        .iter()
        .map(|&c| {
            let slot = &mut relabel[c as usize];
            if *slot == u32::MAX {
                *slot = next;
                next += 1;
            }
            *slot
        })
        .collect();
    Ok(LfrGraph {
        edges,
        communities,
        draws: 1,
        warnings,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn base(n: i64, mu: f64) -> LfrParams {
        LfrParams {
            average_degree: Some(20.0),
            max_degree: Some(50.0),
            min_community: Some(20.0),
            max_community: Some(100.0),
            ..LfrParams::new(n, 2.5, 1.5, mu)
        }
    }

    /// networkx's docstring example: about half of its draws cannot be
    /// placed, a hub needing a community larger than any drawn.
    fn docstring(mu: f64) -> LfrParams {
        LfrParams {
            average_degree: Some(5.0),
            min_community: Some(20.0),
            ..LfrParams::new(250, 3.0, 1.5, mu)
        }
    }

    /// Realized degree, internal degree and community sizes of a graph.
    fn realized(g: &LfrGraph) -> (Vec<u32>, Vec<u32>, Vec<usize>) {
        let n = g.communities.len();
        let mut deg = vec![0; n];
        let mut kin = vec![0; n];
        for &(u, v) in &g.edges {
            deg[u as usize] += 1;
            deg[v as usize] += 1;
            if g.communities[u as usize] == g.communities[v as usize] {
                kin[u as usize] += 1;
                kin[v as usize] += 1;
            }
        }
        let k = g.communities.iter().max().map_or(0, |&c| c as usize + 1);
        let mut sizes = vec![0; k];
        for &c in &g.communities {
            sizes[c as usize] += 1;
        }
        (deg, kin, sizes)
    }

    #[test]
    fn simple_sorted_and_labelled() {
        let g = generate(&base(1000, 0.3), Some(7)).unwrap();
        assert!(
            g.edges.windows(2).all(|w| w[0] < w[1]),
            "sorted, no duplicates"
        );
        assert!(g.edges.iter().all(|&(u, v)| u < v && v < 1000));
        let mut first = Vec::new();
        for &c in &g.communities {
            if !first.contains(&c) {
                assert_eq!(c as usize, first.len(), "labels by first appearance");
                first.push(c);
            }
        }
    }

    #[test]
    fn degrees_mixing_and_sizes_hit_the_request() {
        // mu = tenths / 10, so (1 - mu) k = (10 - tenths) k / 10 exactly
        for tenths in [0u32, 1, 3, 5, 7, 8, 9, 10] {
            let mu = f64::from(tenths) / 10.0;
            let g = generate(&base(1000, mu), Some(3)).unwrap();
            let (deg, kin, sizes) = realized(&g);
            assert!(
                deg.iter().all(|&d| (11..=50).contains(&d)),
                "the law for avg 20 starts at 11"
            );
            let mut off = 0;
            for (v, (&d, &i)) in deg.iter().zip(&kin).enumerate() {
                let num = (10 - tenths) * d;
                let (lo, hi) = (num / 10, num.div_ceil(10));
                let c = g.communities[v] as usize;
                assert!(
                    (i as usize) < sizes[c],
                    "internal degree fits its community"
                );
                if lo < hi {
                    assert!(lo <= i && i <= hi, "node {v}: k {d}, k_in {i}, mu {mu}");
                } else if i != lo {
                    // a stub move, in a community whose targets are all integers
                    assert_eq!(i.abs_diff(lo), 1, "node {v}: k {d}, k_in {i}, mu {mu}");
                    off += 1;
                }
            }
            assert!(off <= sizes.len(), "at most one stub move per community");
            assert!(sizes.iter().all(|&s| (20..=100).contains(&s)));
            assert_eq!(sizes.iter().sum::<usize>(), 1000);
            let out: u32 = deg.iter().zip(&kin).map(|(d, i)| d - i).sum();
            let total: u32 = deg.iter().sum();
            let mixing = f64::from(out) / f64::from(total);
            assert!((mixing - mu).abs() < 0.01, "mu {mu}: realized {mixing}");
            if mu == 0.0 {
                assert_eq!(out, 0, "no edge leaves a community at mu = 0");
            }
            if mu == 1.0 {
                assert!(kin.iter().all(|&i| i == 0), "no edge inside at mu = 1");
            }
        }
    }

    #[test]
    fn average_degree_is_met_in_expectation() {
        // the law for avg 20 has standard deviation 8.84, so the mean of
        // 20 * 1000 degrees has standard error 0.063; 0.3 is about 5 of them
        let mut total = 0u64;
        for seed in 0..20 {
            let g = generate(&base(1000, 0.3), Some(seed)).unwrap();
            total += 2 * g.edges.len() as u64;
        }
        let mean = total as f64 / 20_000.0;
        assert!((mean - 20.0).abs() < 0.3, "{mean}");
    }

    #[test]
    fn same_seed_same_graph() {
        let p = base(500, 0.4);
        assert_eq!(generate(&p, Some(11)), generate(&p, Some(11)));
        assert_ne!(
            generate(&p, Some(11)).unwrap().edges,
            generate(&p, Some(12)).unwrap().edges
        );
        // with redraws too
        for seed in 0..5 {
            assert_eq!(
                generate(&docstring(0.1), Some(seed)),
                generate(&docstring(0.1), Some(seed))
            );
        }
    }

    #[test]
    fn mu_alone_keeps_degrees_and_sizes() {
        let profile = |p: &LfrParams| {
            let (deg, _, mut sizes) = realized(&generate(p, Some(5)).unwrap());
            sizes.sort_unstable();
            (deg, sizes)
        };
        let at_03 = profile(&base(1000, 0.3));
        assert_eq!(at_03, profile(&base(1000, 0.6)));
        // a different degree law draws an independent size sequence
        let other = LfrParams {
            tau1: 3.0,
            ..base(1000, 0.3)
        };
        assert_ne!(at_03.1, profile(&other).1);
    }

    #[test]
    fn failed_draws_are_drawn_again() {
        // networkx 3.4.2 returns a graph for 26 of seeds 0 to 49 here; the
        // first draw fails on about half of them
        for seed in 0..50 {
            let g = generate(&docstring(0.1), Some(seed)).unwrap();
            let (deg, kin, sizes) = realized(&g);
            assert!(sizes.iter().all(|&s| s >= 20));
            for (v, (&k, &i)) in deg.iter().zip(&kin).enumerate() {
                assert!((i as usize) < sizes[g.communities[v] as usize]);
                let x = f64::from(k) * 0.9;
                assert!(f64::from(i) >= x.floor() - 1.0 && f64::from(i) <= x.ceil() + 1.0);
            }
        }
        // every draw fails: degree 10 at mu = 0 needs a community above 10,
        // and only sizes of 10 sum to 100. The error is the assignment's,
        // after max_iters draws
        let never = LfrParams {
            min_degree: Some(10.0),
            max_degree: Some(10.0),
            min_community: Some(10.0),
            max_community: Some(11.0),
            max_iters: 3,
            ..LfrParams::new(100, 2.0, 1.5, 0.0)
        };
        assert_eq!(generate(&never, Some(1)), Err(LfrError::runtime(ASSIGN)));
    }

    /// Every degree `d`, communities of `s`, one draw.
    fn regular(n: i64, s: f64, d: f64, mu: f64) -> LfrParams {
        LfrParams {
            min_degree: Some(d),
            max_degree: Some(d),
            min_community: Some(s),
            max_community: Some(s),
            max_iters: 1,
            ..LfrParams::new(n, 3.0, 1.5, mu)
        }
    }

    #[test]
    fn dense_external_graphs_are_wired() {
        // two communities of 50 at mu = 1: a d-regular bipartite graph exists
        // for every d <= 50; at mu = 0.5 the external degree is d / 2; three
        // communities of 33 with degree 66 at mu = 1 are K(33, 33, 33)
        let cases = [
            (100, 50.0, 35.0, 1.0),
            (100, 50.0, 45.0, 1.0),
            (100, 50.0, 49.0, 1.0),
            (100, 50.0, 50.0, 1.0),
            (100, 50.0, 60.0, 0.5),
            (100, 50.0, 96.0, 0.5),
            (99, 33.0, 66.0, 1.0),
            (100, 25.0, 75.0, 1.0),
        ];
        for (n, s, d, mu) in cases {
            for seed in 0..3 {
                let g = generate(&regular(n, s, d, mu), Some(seed))
                    .unwrap_or_else(|e| panic!("n {n}, size {s}, degree {d}, mu {mu}: {e}"));
                let (deg, kin, _) = realized(&g);
                assert!(deg.iter().all(|&k| f64::from(k) == d));
                let want = (d * (1.0 - mu)) as u32;
                assert!(kin.iter().all(|&i| i == want));
                assert_eq!(g.draws, 1);
            }
        }
    }

    /// 41 nodes of degree 5 or 6 (the law for `average_degree = 5.98` on
    /// `[5, 6]` puts 0.98 on 6) at `mu = 0`, in communities of 6 and 7 nodes:
    /// the only sizes summing to 41 are one 6 and five 7s, so a draw can be
    /// placed only with at least six nodes of degree 5, about one draw in
    /// 7,000, and the graphs returned have a mean of 5.85 or less.
    fn starved(max_iters: i64) -> LfrParams {
        LfrParams {
            average_degree: Some(5.98),
            max_degree: Some(6.0),
            min_community: Some(6.0),
            max_community: Some(7.0),
            max_iters,
            ..LfrParams::new(41, 2.0, 1.5, 0.0)
        }
    }

    #[test]
    fn a_shortfall_from_rejected_draws_is_reported() {
        for seed in 0..10 {
            let g = generate(&starved(100_000), Some(seed)).unwrap();
            let mean = 2.0 * g.edges.len() as f64 / 41.0;
            assert!(mean <= 5.86 && g.draws > 1, "{mean} {}", g.draws);
            // mu = 0 rounds nothing, so the mixing has no shortfall to report
            let [m] = g.warnings.as_slice() else {
                panic!("one warning, 2% and 5.8 standard errors below the request")
            };
            assert!(m.contains(&format!("is {mean:.3}, ")), "{m}");
            assert!(m.contains("below average_degree=5.98"), "{m}");
            assert!(
                m.contains(&format!("{} of the {} draws", g.draws - 1, g.draws)),
                "{m}"
            );
            assert!(m.contains("max_degree"), "{m}");
        }
        // the docstring example misses by about 2% on average, within four
        // standard errors (0.37 at n = 250); the grid request by far less
        for seed in 0..20 {
            let quiet = |p: &LfrParams| generate(p, Some(seed)).unwrap().warnings.is_empty();
            assert!(quiet(&docstring(0.1)));
            assert!(quiet(&base(1000, 0.3)));
        }
        // a min_degree request has no average to miss
        let fixed = LfrParams {
            average_degree: None,
            min_degree: Some(5.0),
            ..starved(100_000)
        };
        assert_eq!(generate(&fixed, Some(0)).map(|g| g.warnings), Ok(vec![]));
        // a request that tol accepts far outside the range uses the law on
        // [1, max_degree], and the graph is compared with that law's mean
        let outside = LfrParams {
            average_degree: Some(-1e299),
            tol: 1e300,
            ..base(1000, 0.3)
        };
        for seed in 0..5 {
            assert!(generate(&outside, Some(seed)).unwrap().warnings.is_empty());
        }
    }

    #[test]
    fn shortfall_needs_both_margins() {
        let (j, alpha) = powerlaw::law_for_average(2.0, 5.98, 6, 1e-7).unwrap();
        let law = powerlaw::PowerLaw::mixture(2.0, j, 6, alpha);
        let (_, sd) = law.moments();
        let graph = |n: usize, m: usize, draws: usize| LfrGraph {
            edges: vec![(0, 1); m],
            communities: vec![0; n],
            draws,
            warnings: Vec::new(),
        };
        // sd 0.140: at n = 41 the standard error is 0.0219, four of them
        // 0.088, and 2% of 5.98 is 0.120
        assert!((sd - 0.1403).abs() < 1e-3, "{sd}");
        // 5.854 (0.126 below, 5.8 standard errors): both margins exceeded
        assert!(shortfall(&graph(41, 120, 9), 5.98, &law).is_some());
        // 5.878 at n = 82 (0.102 below, 6.6 standard errors): within 2%
        assert!(shortfall(&graph(82, 241, 9), 5.98, &law).is_none());
        // 5.833 at n = 12 (0.147 below, 3.6 standard errors): within four
        assert!(shortfall(&graph(12, 35, 9), 5.98, &law).is_none());
        // above the request, and without a rejected draw
        let m = shortfall(&graph(41, 130, 1), 5.98, &law).unwrap();
        assert!(
            m.contains("above average_degree=5.98") && m.contains("No draw was rejected"),
            "{m}"
        );
    }

    /// Fraction of the edges of `g` that join two communities.
    fn edge_mixing(g: &LfrGraph) -> f64 {
        let c = &g.communities;
        let between = g
            .edges
            .iter()
            .filter(|&&(u, v)| c[u as usize] != c[v as usize])
            .count();
        between as f64 / g.edges.len() as f64
    }

    #[test]
    fn a_mixing_shortfall_is_reported() {
        // every degree 2 at mu = 0.3 (1.4, which networkx rounds to 1) with
        // the default community bounds, so communities of 2: no node has room
        // for a second internal edge, so the 40% that rounded up are rounded
        // down and nothing takes their edge back, a shift of 0.2
        let pairs = LfrParams {
            min_degree: Some(2.0),
            max_degree: Some(2.0),
            ..LfrParams::new(10_000, 2.5, 1.5, 0.3)
        };
        // communities of 10 nodes of degree 10 at mu = 0.06 (9.4): every
        // node is held to 9, a shift of 0.04
        let tens = LfrParams {
            min_degree: Some(10.0),
            max_degree: Some(10.0),
            min_community: Some(10.0),
            max_community: Some(10.0),
            ..LfrParams::new(2000, 2.0, 1.5, 0.06)
        };
        for (p, mixing, shift) in [(&pairs, 0.5, 0.2), (&tens, 0.1, 0.04)] {
            for seed in 0..3 {
                let g = generate(p, Some(seed)).unwrap();
                let realized = edge_mixing(&g);
                assert!((realized - mixing).abs() < 1e-12, "{realized}");
                let [m] = g.warnings.as_slice() else {
                    panic!("one warning: {:?}", g.warnings)
                };
                assert!(
                    m.contains(&format!("graph is {realized:.4}, against mu={}", p.mu)),
                    "{m}"
                );
                assert!(m.contains("min_community"), "{m}");
                // the shift stated is the unpaid ends over the degree sum
                let at = m.find("which adds ").unwrap() + "which adds ".len();
                let stated: f64 = m[at..at + 6].parse().unwrap();
                assert!((stated - shift).abs() < 0.02, "{m}");
            }
        }
        // the grid's requests at every mu, and a request whose nodes held to
        // their floor are repaid by other communities, are not reported
        for tenths in 0..=10 {
            let g = generate(&base(1000, f64::from(tenths) / 10.0), Some(1)).unwrap();
            assert!(g.warnings.is_empty(), "{:?}", g.warnings);
        }
        let repaid = LfrParams {
            min_degree: Some(5.0),
            ..LfrParams::new(10_000, 3.5, 3.0, 0.1)
        };
        let g = generate(&repaid, Some(0)).unwrap();
        assert!(g.warnings.is_empty(), "{:?}", g.warnings);
        assert!((edge_mixing(&g) - 0.1).abs() < 0.003);
    }

    #[test]
    fn mixing_shortfall_needs_more_than_a_thousandth() {
        let edges = [(0, 1), (1, 2)];
        let comm = [0, 0, 1];
        assert!(mixing_shortfall(10, 10_000, 0.3, &edges, &comm).is_none());
        let m = mixing_shortfall(11, 10_000, 0.3, &edges, &comm).unwrap();
        assert!(m.contains("is 0.5000, against mu=0.3"), "{m}");
        assert!(
            m.contains("for 11 of them") && m.contains("adds 0.0011"),
            "{m}"
        );
    }

    #[test]
    fn an_odd_degree_sum_is_drawn_again() {
        // degrees 3 or 4 on 101 nodes: with max_iters = 2 a draw's two
        // attempts at an even sum both fail one time in four. Raised at once,
        // that error would end at least 500 of 2,000 calls (sd 19); drawn
        // again, it ends a call only if the second draw also fails there:
        // 0.16 of calls here, 320 (sd 16)
        let p = LfrParams {
            min_degree: Some(3.0),
            max_degree: Some(4.0),
            min_community: Some(20.0),
            max_community: Some(101.0),
            max_iters: 2,
            ..LfrParams::new(101, 2.5, 1.5, 0.3)
        };
        let (mut power_law, mut second) = (0, 0);
        for seed in 0..2000 {
            match generate(&p, Some(seed)) {
                Err(e) if e == LfrError::runtime(POWER_LAW) => power_law += 1,
                Ok(g) if g.draws == 2 => second += 1,
                _ => {}
            }
        }
        assert!(power_law < 410, "{power_law}");
        assert!(second > 0);
    }

    #[test]
    fn infeasible_requests_are_runtime_errors() {
        // every degree is 3 and n is odd: the degree sum is always odd
        let odd = LfrParams {
            min_degree: Some(3.0),
            max_degree: Some(3.0),
            min_community: Some(11.0),
            max_community: Some(11.0),
            ..LfrParams::new(11, 3.0, 1.5, 0.1)
        };
        assert_eq!(generate(&odd, Some(1)), Err(LfrError::runtime(POWER_LAW)));
        // no size sequence of 30s sums to 100
        let sizes = LfrParams {
            min_degree: Some(5.0),
            max_degree: Some(10.0),
            min_community: Some(30.0),
            max_community: Some(30.0),
            ..LfrParams::new(100, 3.0, 1.5, 0.1)
        };
        assert_eq!(generate(&sizes, Some(1)), Err(LfrError::runtime(POWER_LAW)));
        // internal degree 25..30 cannot fit communities of 20
        let fit = LfrParams {
            min_degree: Some(25.0),
            max_degree: Some(30.0),
            min_community: Some(20.0),
            max_community: Some(20.0),
            ..LfrParams::new(100, 3.0, 1.5, 0.0)
        };
        assert_eq!(generate(&fit, Some(1)), Err(LfrError::runtime(ASSIGN)));
        // min_degree = n: networkx's errors where it stops (an odd degree sum;
        // a community not above round((1 - mu) n)), a power-law error where it
        // hangs (every community above it)
        let degree_n = |n: i64, mu: f64, lo: Option<f64>, hi: Option<f64>| LfrParams {
            min_degree: Some(n as f64),
            min_community: lo,
            max_community: hi,
            ..LfrParams::new(n, 2.5, 1.5, mu)
        };
        let power_law = Err(LfrError::runtime(POWER_LAW));
        let assign = Err(LfrError::runtime(ASSIGN));
        assert_eq!(generate(&degree_n(1, 0.3, None, None), Some(1)), power_law);
        assert_eq!(generate(&degree_n(51, 0.3, None, None), Some(1)), power_law);
        assert_eq!(generate(&degree_n(50, 0.0, None, None), Some(1)), assign);
        assert_eq!(
            generate(&degree_n(50, 0.3, Some(10.0), Some(25.0)), Some(1)),
            assign
        );
        assert_eq!(generate(&degree_n(50, 0.3, None, None), Some(1)), power_law);
        // every community size drawn is above n
        let big = LfrParams {
            min_degree: Some(3.0),
            max_degree: Some(20.0),
            min_community: Some(101.0),
            max_community: Some(120.0),
            ..LfrParams::new(100, 2.5, 1.5, 0.3)
        };
        assert_eq!(generate(&big, Some(1)), power_law);
        // a default max_community below min_community: networkx hangs
        let crossed = LfrParams {
            max_community: None,
            ..big
        };
        assert_eq!(generate(&crossed, Some(1)), power_law);
        // a single community cannot host edges between communities
        let one = LfrParams {
            min_degree: Some(5.0),
            max_degree: Some(10.0),
            min_community: Some(60.0),
            max_community: Some(60.0),
            ..LfrParams::new(60, 3.0, 1.5, 0.5)
        };
        assert_eq!(generate(&one, Some(1)), Err(LfrError::runtime(WIRE)));
    }
}
