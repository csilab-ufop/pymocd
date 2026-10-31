//! Edges from an `Assignment`: every node gets exactly its internal degree
//! inside its community and the rest of its degree across communities, with no
//! self-loop and no repeated edge.
//!
//! Inside a community the internal degrees are graphical (`membership` made
//! them so), and Havel–Hakimi builds a simple graph on them; it is then
//! randomized by degree-preserving double-edge swaps, ten per edge, the
//! original code's ten sweeps of switches (benchm.cpp, `build_subgraph`).
//! Between communities a configuration model pairs the external stubs, and each
//! pair that is a loop, a duplicate, or inside one community is rewired against
//! another pair, the move of Lancichinetti & Fortunato, Phys. Rev. E 80,
//! 016118 (2009), Fig. 2; ten swaps per edge then randomize the result.
//! Every move preserves every node's internal and external degree, which is
//! what makes the mixing exact.
//!
//! When the external graph is dense and spread over few communities, the
//! rewiring walk stalls: with two communities a pair inside one can only be
//! fixed against a pair inside the other, and near a complete bipartite graph
//! almost every rewiring repeats an edge. The external edges are then built
//! greedily instead (`greedy_external`), and the same swaps randomize them.
//! This Source Code Form is subject to the terms of The GNU General Public License v3.0
//! Copyright 2026 - Guilherme Santos. If a copy of the MPL was not distributed with this
//! file, You can obtain one at https://www.gnu.org/licenses/gpl-3.0.html

use std::collections::{BTreeSet, BinaryHeap};

use rand::seq::SliceRandom;
use rand::{Rng, RngExt};
use rustc_hash::{FxBuildHasher, FxHashMap, FxHashSet};

use super::membership::Assignment;
use super::{ASSIGN, LfrError, WIRE};

/// Swap attempts per edge when randomizing.
const SWAPS_PER_EDGE: usize = 10;
/// Fresh stub shuffles the external wiring may start from.
const RESHUFFLES: usize = 3;

type Edge = (u32, u32);

#[inline]
const fn ordered(u: u32, v: u32) -> Edge {
    if u < v { (u, v) } else { (v, u) }
}

#[inline]
fn key((u, v): Edge) -> u64 {
    let (a, b) = ordered(u, v);
    (u64::from(a) << 32) | u64::from(b)
}

/// All edges of the graph, unsorted.
///
/// An edge inside a community never equals an edge between communities, nor
/// one inside another community, so each community's swaps only check its own
/// edges and the external swaps only the external ones: the same decisions as
/// one set of every edge, from sets small enough to stay in cache.
pub fn wire<R: Rng + ?Sized>(a: &Assignment, rng: &mut R) -> Result<Vec<Edge>, LfrError> {
    let m = a.deg.iter().map(|&k| k as usize).sum::<usize>() / 2;
    let mut present: FxHashSet<u64> = FxHashSet::default();
    let mut edges = Vec::with_capacity(m);
    let mut local = Vec::new();
    for members in &a.members {
        local.clear();
        present.clear();
        // membership guarantees graphical sequences; failing here is a bug,
        // reported rather than panicking across the FFI boundary
        if !havel_hakimi(members, &a.kin, &mut local, rng) {
            return Err(LfrError::runtime(ASSIGN));
        }
        present.extend(local.iter().map(|&e| key(e)));
        randomize(&mut local, &mut present, |_, _| true, rng);
        edges.extend_from_slice(&local);
    }
    let mut external = external_edges(a, rng)?;
    let mut present: FxHashSet<u64> =
        FxHashSet::with_capacity_and_hasher(external.len(), FxBuildHasher);
    present.extend(external.iter().map(|&e| key(e)));
    let comm = &a.comm;
    randomize(
        &mut external,
        &mut present,
        |u, v| comm[u as usize] != comm[v as usize],
        rng,
    );
    edges.append(&mut external);
    Ok(edges)
}

/// Havel–Hakimi on `nodes` with degrees `degree[v]`: the node of largest
/// remaining degree `d` is joined to the `d` others of largest remaining
/// degree. On a graphical sequence this never runs short. Ties are broken by
/// a random key per node, so equal degrees are not wired in node order.
fn havel_hakimi<R: Rng + ?Sized>(
    nodes: &[u32],
    degree: &[u32],
    out: &mut Vec<Edge>,
    rng: &mut R,
) -> bool {
    let mut left: BTreeSet<(u32, u32, u32)> = BTreeSet::new();
    for &v in nodes {
        let d = degree[v as usize];
        if d > 0 {
            left.insert((d, rng.random::<u32>(), v));
        }
    }
    let mut partners = Vec::new();
    while let Some((d, _, v)) = left.pop_last() {
        if d as usize > left.len() {
            return false;
        }
        partners.clear();
        for _ in 0..d {
            partners.push(left.pop_last().expect("d <= len"));
        }
        for &(r, tie, u) in &partners {
            out.push(ordered(u, v));
            if r > 1 {
                left.insert((r - 1, tie, u));
            }
        }
    }
    true
}

/// Degree-preserving double-edge swaps, `SWAPS_PER_EDGE` attempts per edge:
/// `(a, b), (c, d) -> (a, d), (c, b)` with the second edge's orientation at
/// random, kept only if neither new edge is a loop, exists already, or is
/// refused by `allowed`.
fn randomize<R, F>(edges: &mut [Edge], present: &mut FxHashSet<u64>, allowed: F, rng: &mut R)
where
    R: Rng + ?Sized,
    F: Fn(u32, u32) -> bool,
{
    let m = edges.len();
    if m < 2 {
        return;
    }
    for _ in 0..SWAPS_PER_EDGE * m {
        let i = rng.random_range(0..m);
        let j = rng.random_range(0..m);
        if i == j {
            continue;
        }
        let (a, b) = edges[i];
        let (mut c, mut d) = edges[j];
        if rng.random_bool(0.5) {
            std::mem::swap(&mut c, &mut d);
        }
        if a == d || c == b || !allowed(a, d) || !allowed(c, b) {
            continue;
        }
        let (e1, e2) = (ordered(a, d), ordered(c, b));
        if present.contains(&key(e1)) || present.contains(&key(e2)) {
            continue;
        }
        present.remove(&key(edges[i]));
        present.remove(&key(edges[j]));
        present.insert(key(e1));
        present.insert(key(e2));
        edges[i] = e1;
        edges[j] = e2;
    }
}

/// Pairs the external stubs into a simple graph with no edge inside a
/// community: a shuffled stub list, then rewiring of each bad pair, within
/// `50 m + 1000` attempts per shuffle, and `greedy_external` when no shuffle
/// gets there.
///
/// A bad pair inside a community can only be rewired against a pair whose
/// two ends both lie outside it, and a draw from all `m` pairs rarely finds
/// one when that community holds about half of the external stubs: then
/// nearly every edge of the finished graph has one end in it, and the pairs
/// that avoid it are few. Half of the partners of a pair inside the
/// community with the most external stubs are therefore drawn from the pairs
/// that avoid it, kept in `Avoiding`; half of the partners of any other bad
/// pair from the pairs still waiting for repair, which with two communities
/// are the pairs inside the other one. The rest come from all pairs.
fn external_edges<R: Rng + ?Sized>(a: &Assignment, rng: &mut R) -> Result<Vec<Edge>, LfrError> {
    let q = a.members.len();
    let mut stubs = Vec::new();
    let mut held = vec![0u64; q];
    for (v, (&k, &i)) in a.deg.iter().zip(&a.kin).enumerate() {
        stubs.extend(std::iter::repeat_n(v as u32, (k - i) as usize));
        held[a.comm[v] as usize] += u64::from(k - i);
    }
    let m = stubs.len() / 2;
    if m == 0 {
        return Ok(Vec::new());
    }
    let comm = &a.comm;
    let same = |u: u32, v: u32| comm[u as usize] == comm[v as usize];
    // the first of the communities holding the most external stubs
    let heavy = (0..q).fold(0, |best, c| if held[c] > held[best] { c } else { best }) as u32;
    let inside_heavy = |(u, v): Edge| comm[u as usize] == heavy && comm[v as usize] == heavy;
    let budget = 50 * m + 1000;
    for _ in 0..RESHUFFLES {
        stubs.shuffle(rng);
        let mut edges: Vec<Edge> = stubs.chunks_exact(2).map(|p| ordered(p[0], p[1])).collect();
        let mut count: FxHashMap<u64, u32> = FxHashMap::with_capacity_and_hasher(m, FxBuildHasher);
        for &e in &edges {
            *count.entry(key(e)).or_insert(0) += 1;
        }
        let bad = |e: Edge, count: &FxHashMap<u64, u32>| {
            e.0 == e.1 || same(e.0, e.1) || count[&key(e)] > 1
        };
        let mut avoiding = Avoiding::new(&edges, |e| !touches(e, comm, heavy));
        let mut todo: Vec<usize> = (0..m).filter(|&i| bad(edges[i], &count)).collect();
        let mut attempts = 0;
        let mut done = true;
        'repair: while let Some(i) = todo.pop() {
            // an earlier rewiring may have fixed it (the other copy of a duplicate)
            if !bad(edges[i], &count) {
                continue;
            }
            loop {
                attempts += 1;
                if attempts > budget {
                    done = false;
                    break 'repair;
                }
                let targeted = rng.random_bool(0.5);
                let j = if targeted && inside_heavy(edges[i]) && !avoiding.is_empty() {
                    avoiding.draw(rng)
                } else if targeted && !inside_heavy(edges[i]) && !todo.is_empty() {
                    todo[rng.random_range(0..todo.len())]
                } else {
                    rng.random_range(0..m)
                };
                if j == i {
                    continue;
                }
                let (x, y) = edges[i];
                let (mut u, mut w) = edges[j];
                if rng.random_bool(0.5) {
                    std::mem::swap(&mut u, &mut w);
                }
                // (x, y), (u, w) -> (x, u), (y, w)
                if x == u || y == w || same(x, u) || same(y, w) {
                    continue;
                }
                let (e1, e2) = (ordered(x, u), ordered(y, w));
                if e1 == e2 || count.contains_key(&key(e1)) || count.contains_key(&key(e2)) {
                    continue;
                }
                for old in [edges[i], edges[j]] {
                    let k = key(old);
                    let c = count.get_mut(&k).expect("listed edges are counted");
                    *c -= 1;
                    if *c == 0 {
                        count.remove(&k);
                    }
                }
                edges[i] = e1;
                edges[j] = e2;
                count.insert(key(e1), 1);
                count.insert(key(e2), 1);
                avoiding.set(i, !touches(e1, comm, heavy));
                avoiding.set(j, !touches(e2, comm, heavy));
                break;
            }
        }
        if done {
            return Ok(edges);
        }
    }
    greedy_external(a, rng).ok_or_else(|| LfrError::runtime(WIRE))
}

/// Whether either end of `e` lies in community `c`.
fn touches((u, v): Edge, comm: &[u32], c: u32) -> bool {
    comm[u as usize] == c || comm[v as usize] == c
}

/// A set of pair indices with constant-time insertion, removal and uniform
/// draws: the pairs with no end in one community.
struct Avoiding {
    list: Vec<usize>,
    /// Position of each pair in `list`, or `usize::MAX`.
    slot: Vec<usize>,
}

impl Avoiding {
    fn new(edges: &[Edge], keep: impl Fn(Edge) -> bool) -> Self {
        let mut list = Vec::new();
        let mut slot = vec![usize::MAX; edges.len()];
        for (i, &e) in edges.iter().enumerate() {
            if keep(e) {
                slot[i] = list.len();
                list.push(i);
            }
        }
        Self { list, slot }
    }

    fn is_empty(&self) -> bool {
        self.list.is_empty()
    }

    fn draw<R: Rng + ?Sized>(&self, rng: &mut R) -> usize {
        self.list[rng.random_range(0..self.list.len())]
    }

    fn set(&mut self, i: usize, member: bool) {
        let at = self.slot[i];
        if member && at == usize::MAX {
            self.slot[i] = self.list.len();
            self.list.push(i);
        } else if !member && at != usize::MAX {
            let last = *self.list.last().expect("i is listed");
            self.list.swap_remove(at);
            if last != i {
                self.slot[last] = at;
            }
            self.slot[i] = usize::MAX;
        }
    }
}

/// A node's remaining external degree, a random tie-break, and the node.
type Slot = (u32, u32, u32);

/// Builds the external edges by a Havel–Hakimi rule for graphs with no edge
/// inside a community: the node with the most external stubs left, in the
/// community with the most stubs left, is joined to the nodes with the most
/// stubs left in other communities. None when a node runs out of such
/// partners. Ties are broken at random.
///
/// With two communities this is the greedy construction of a bipartite graph
/// with given degrees (Gale, Ryser), which succeeds whenever such a graph
/// exists, whichever node is joined first: if one exists in which `v` misses
/// some `y` of larger degree than a neighbour `y'`, then `y` has a neighbour
/// `z` that `y'` lacks, and the swap `(v, y'), (z, y) -> (v, y), (z, y')`
/// keeps every degree. With more communities `z` may share `y'`'s community,
/// so the rule can fail where a graph exists, and the draw is then rejected.
/// Serving the community with the most stubs first keeps one that holds
/// about half of them from being left with no partners.
///
/// Each community keeps its nodes ordered by stubs left, and a node's
/// partners are merged from the communities' tops, so a node costs about
/// `d log n` for `d` partners, not a scan past the larger nodes of its own
/// community.
fn greedy_external<R: Rng + ?Sized>(a: &Assignment, rng: &mut R) -> Option<Vec<Edge>> {
    let q = a.members.len();
    let mut parts: Vec<BTreeSet<Slot>> = vec![BTreeSet::new(); q];
    let mut left = vec![0u64; q];
    for (v, (&k, &i)) in a.deg.iter().zip(&a.kin).enumerate() {
        if k > i {
            let c = a.comm[v] as usize;
            parts[c].insert((k - i, rng.random(), v as u32));
            left[c] += u64::from(k - i);
        }
    }
    let m = left.iter().sum::<u64>() as usize / 2;
    // each community's largest slot, and the communities by stubs left
    let mut tops: BTreeSet<(Slot, u32)> = BTreeSet::new();
    let mut by_left: BTreeSet<(u64, u32)> = BTreeSet::new();
    for (c, part) in parts.iter().enumerate() {
        if let Some(&top) = part.last() {
            tops.insert((top, c as u32));
            by_left.insert((left[c], c as u32));
        }
    }
    let mut edges = Vec::with_capacity(m);
    let mut picked: Vec<(Slot, u32)> = Vec::new();
    let mut heap: BinaryHeap<(Slot, u32)> = BinaryHeap::new();
    let mut touched: Vec<u32> = Vec::new();
    let mut is_touched = vec![false; q];
    while let Some(&(_, p)) = by_left.last() {
        let pi = p as usize;
        let source = *parts[pi]
            .last()
            .expect("a community with stubs left has a node");
        let (d, _, v) = source;
        picked.clear();
        heap.clear();
        {
            let mut unopened = tops.iter().rev().filter(|&&(_, c)| c != p).peekable();
            while picked.len() < d as usize {
                let from_heap = match (heap.peek(), unopened.peek()) {
                    (Some(h), Some(o)) => h.0 > o.0,
                    (Some(_), None) => true,
                    (None, Some(_)) => false,
                    (None, None) => return None,
                };
                if from_heap {
                    let (slot, c) = heap.pop().expect("peeked");
                    picked.push((slot, c));
                    if let Some(&below) = parts[c as usize].range(..slot).next_back() {
                        heap.push((below, c));
                    }
                } else {
                    let &(slot, c) = unopened.next().expect("peeked");
                    heap.push((slot, c));
                }
            }
        }
        // the source leaves its community, every partner loses one stub
        let mut retire =
            |c: u32, tops: &mut BTreeSet<(Slot, u32)>, by_left: &mut BTreeSet<(u64, u32)>| {
                let ci = c as usize;
                if !is_touched[ci] {
                    is_touched[ci] = true;
                    touched.push(c);
                    let top = *parts[ci].last().expect("a touched community has a node");
                    tops.remove(&(top, c));
                    by_left.remove(&(left[ci], c));
                }
            };
        retire(p, &mut tops, &mut by_left);
        for &(_, c) in &picked {
            retire(c, &mut tops, &mut by_left);
        }
        parts[pi].remove(&source);
        left[pi] -= u64::from(d);
        for &(slot, c) in &picked {
            let ci = c as usize;
            parts[ci].remove(&slot);
            left[ci] -= 1;
            let (stubs, tie, u) = slot;
            if stubs > 1 {
                parts[ci].insert((stubs - 1, tie, u));
            }
            edges.push(ordered(u, v));
        }
        for c in touched.drain(..) {
            let ci = c as usize;
            is_touched[ci] = false;
            if let Some(&top) = parts[ci].last() {
                tops.insert((top, c));
                by_left.insert((left[ci], c));
            }
        }
    }
    Some(edges)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    #[test]
    fn havel_hakimi_realizes_the_sequence() {
        let mut rng = StdRng::seed_from_u64(2);
        let nodes: Vec<u32> = (10..16).collect();
        let mut degree = vec![0u32; 16];
        for (v, d) in nodes.iter().zip([5, 3, 3, 2, 2, 1]) {
            degree[*v as usize] = d;
        }
        let mut out = Vec::new();
        assert!(havel_hakimi(&nodes, &degree, &mut out, &mut rng));
        let mut got = vec![0u32; 16];
        for &(u, v) in &out {
            assert!(u < v);
            got[u as usize] += 1;
            got[v as usize] += 1;
        }
        assert_eq!(got, degree);
        let mut keys: Vec<u64> = out.iter().map(|&e| key(e)).collect();
        keys.sort_unstable();
        keys.dedup();
        assert_eq!(keys.len(), out.len(), "no repeated edge");
        // 3, 3, 1, 1 is not graphical
        let mut degree = vec![0u32; 16];
        for (v, d) in nodes.iter().zip([3, 3, 1, 1, 0, 0]) {
            degree[*v as usize] = d;
        }
        assert!(!havel_hakimi(&nodes, &degree, &mut Vec::new(), &mut rng));
    }

    #[test]
    fn swaps_preserve_degrees_and_the_constraint() {
        let mut rng = StdRng::seed_from_u64(4);
        // a 20-cycle between the even and the odd nodes' "communities"
        let comm: Vec<u32> = (0..20).map(|v| v % 2).collect();
        let mut edges: Vec<Edge> = (0..20).map(|v| ordered(v, (v + 1) % 20)).collect();
        let mut present: FxHashSet<u64> = edges.iter().map(|&e| key(e)).collect();
        randomize(
            &mut edges,
            &mut present,
            |u, v| comm[u as usize] != comm[v as usize],
            &mut rng,
        );
        let mut deg = [0u32; 20];
        for &(u, v) in &edges {
            assert_ne!(u, v);
            assert_ne!(comm[u as usize], comm[v as usize]);
            deg[u as usize] += 1;
            deg[v as usize] += 1;
        }
        assert!(deg.iter().all(|&d| d == 2));
        assert_eq!(present.len(), 20, "set and list agree, no duplicate");
        assert!(edges.iter().all(|&e| present.contains(&key(e))));
    }
}
