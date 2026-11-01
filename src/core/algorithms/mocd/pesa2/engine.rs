//! The PESA-II generational loop.
//! This Source Code Form is subject to the terms of The GNU General Public License v3.0
//! Copyright 2025 - Guilherme Santos. If a copy of the MPL was not distributed with this
//! file, You can obtain one at https://www.gnu.org/licenses/gpl-3.0.html

use rand::RngExt;

use crate::core::algorithms::mocd::locus::{Genome, NodeIndex, random_genome};
use crate::core::algorithms::mocd::operators::{mutate, mutation_size, uniform_crossover};
use crate::core::graph::Graph;

use super::archive::{insert_nondominated, truncate};
use super::grid::{GRID_DIVISIONS, assign_cells, select_member};
use super::solution::{Member, Solution, evaluate};

/// Runs PESA-II and returns the final external archive as the Pareto front.
/// `pop_size` is PESA-II's `ipsize`, `ep_size` its `epsize`.
pub fn evolutionary_phase(
    graph: &Graph,
    debug_level: i8,
    num_gens: usize,
    pop_size: usize,
    ep_size: usize,
    cross_rate: f64,
    mut_rate: f64,
) -> Vec<Solution> {
    if graph.nodes.is_empty() || graph.edges.is_empty() {
        return Vec::new();
    }

    for p in [cross_rate, mut_rate] {
        assert!((0.0..=1.0).contains(&p), "p={p:?} is outside range [0.0, 1.0]");
    }

    let idx = NodeIndex::build(graph);
    let degrees: Vec<usize> = idx.index_to_node.iter().map(|n| graph.degree(n)).collect();
    let ipsize = pop_size.max(1);
    let epsize = ep_size.max(1);
    let genes = mutation_size(idx.n());

    let mut rng = rand::rng();

    let mut ep = seed_archive(graph, &idx, &degrees, ipsize, epsize, &mut rng);

    let mut occ: Vec<usize> = Vec::new();
    for generation in 0..num_gens {
        if ep.is_empty() {
            break;
        }

        assign_cells(&mut ep, GRID_DIVISIONS, &mut occ);
        let occupied = occupied_cells(&occ);

        let new_ip = breed(
            &ep, &occ, &occupied, &idx, ipsize, cross_rate, mut_rate, genes, &mut rng,
        );

        let evaluated: Vec<Member> = new_ip
            .into_iter()
            .map(|g| evaluate(graph, &idx, &degrees, g))
            .collect();

        for m in evaluated {
            insert_nondominated(&mut ep, m);
        }

        if ep.len() > epsize {
            truncate(&mut ep, epsize, &mut rng);
        }

        if debug_level >= 1 {
            crate::debug!(debug, "gen {} | EP size: {}", generation, ep.len());
        }
    }

    ep.into_iter().map(|m| m.solution).collect()
}

/// Generation 0 (Corne et al. 2001, §3): a random internal population, of
/// which the non-dominated members become the archive.
fn seed_archive(
    graph: &Graph,
    idx: &NodeIndex,
    degrees: &[usize],
    ipsize: usize,
    epsize: usize,
    rng: &mut impl rand::Rng,
) -> Vec<Member> {
    let initial_ip: Vec<Member> = (0..ipsize)
        .map(|_| evaluate(graph, idx, degrees, random_genome(idx, rng)))
        .collect();

    let mut ep: Vec<Member> = Vec::new();
    for m in initial_ip {
        insert_nondominated(&mut ep, m);
    }
    if ep.len() > epsize {
        truncate(&mut ep, epsize, rng);
    }
    ep
}

fn occupied_cells(occ: &[usize]) -> Vec<usize> {
    occ.iter()
        .enumerate()
        .filter(|&(_, &c)| c > 0)
        .map(|(cell, _)| cell)
        .collect()
}

/// One generation of offspring, as the authors' `GA::runGA` breeds it:
/// `ipsize` clones of archive members, each taken from a uniformly random
/// occupied cell and mutated; then `ceil(ipsize / 2)` variation steps, each of
/// which mutates one random clone with probability `mut_rate` or else, with
/// probability `cross_rate`, crosses two distinct random clones in place.
#[allow(clippy::too_many_arguments)]
fn breed(
    ep: &[Member],
    occ: &[usize],
    occupied: &[usize],
    idx: &NodeIndex,
    ipsize: usize,
    cross_rate: f64,
    mut_rate: f64,
    genes: usize,
    rng: &mut impl rand::Rng,
) -> Vec<Genome> {
    let mut ip: Vec<Genome> = (0..ipsize)
        .map(|_| {
            let mut g = select_member(ep, occ, occupied, rng).genome.clone();
            mutate(&mut g, idx, genes, rng);
            g
        })
        .collect();
    for _ in (0..ipsize).step_by(2) {
        let r: f64 = rng.random();
        if r < mut_rate {
            let i = rng.random_range(0..ipsize);
            mutate(&mut ip[i], idx, genes, rng);
        } else if r < cross_rate + mut_rate && ipsize >= 2 {
            let i = rng.random_range(0..ipsize);
            let mut j = rng.random_range(0..ipsize - 1);
            if j >= i {
                j += 1;
            }
            let (head, tail) = ip.split_at_mut(i.max(j));
            uniform_crossover(&mut head[i.min(j)], &mut tail[0], rng);
        }
    }
    ip
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::algorithms::mocd::EPSIZE_CAP;
    use crate::core::graph::Graph;
    use std::collections::HashSet;

    fn two_triangles() -> Graph {
        let mut g = Graph::new();
        for (a, b) in [(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5), (2, 3)] {
            g.add_edge(a, b);
        }
        g.finalize();
        g
    }

    #[test]
    fn archive_is_pareto_nondominated() {
        let g = two_triangles();
        let front = evolutionary_phase(&g, 0, 30, 30, 30, 0.6, 0.4);
        assert!(!front.is_empty());
        for (i, a) in front.iter().enumerate() {
            for (j, b) in front.iter().enumerate() {
                if i != j {
                    assert!(!a.dominates(b), "front member {i} dominates {j}");
                }
            }
        }
    }

    #[test]
    fn breed_returns_ipsize_legal_genomes() {
        let g = two_triangles();
        let idx = NodeIndex::build(&g);
        let degrees: Vec<usize> = idx.index_to_node.iter().map(|n| g.degree(n)).collect();
        let mut rng = rand::rng();
        let mut ep: Vec<Member> = (0..10)
            .map(|_| evaluate(&g, &idx, &degrees, random_genome(&idx, &mut rng)))
            .collect();
        let mut occ = Vec::new();
        assign_cells(&mut ep, GRID_DIVISIONS, &mut occ);
        let occupied = occupied_cells(&occ);
        for ipsize in [1, 2, 7, 30] {
            let ip = breed(&ep, &occ, &occupied, &idx, ipsize, 0.6, 0.4, 3, &mut rng);
            assert_eq!(ip.len(), ipsize);
            for genome in &ip {
                for (i, &a) in genome.iter().enumerate() {
                    assert!(idx.neighbor_candidates[i].contains(&a));
                }
            }
        }
    }

    #[test]
    #[should_panic(expected = "outside range")]
    fn rate_outside_unit_interval_panics() {
        evolutionary_phase(&two_triangles(), 0, 5, 10, 10, 1.5, 0.4);
    }

    #[test]
    fn epsize_cap_respected() {
        let g = two_triangles();
        let front = evolutionary_phase(&g, 0, 5, 200, EPSIZE_CAP, 0.6, 0.4);
        assert!(front.len() <= EPSIZE_CAP);
    }

    #[test]
    fn shi_mocd_max_q_is_two_community_split() {
        let g = two_triangles();
        let archive = evolutionary_phase(&g, 0, 100, 100, 100, 0.6, 0.4);
        assert!(!archive.is_empty(), "empty PESA-II archive");
        // MOCD-Q (Shi Eq. 3.8): argmin(intra + inter) = argmax Q, since
        // Q = 1 - intra - inter (Eq. 3.7) and objectives = [intra, inter].
        let best = archive
            .iter()
            .min_by(|a, b| {
                (a.objectives[0] + a.objectives[1])
                    .partial_cmp(&(b.objectives[0] + b.objectives[1]))
                    .unwrap()
            })
            .unwrap();
        let q = 1.0 - best.objectives[0] - best.objectives[1];
        assert!(q > 0.0, "Q = {q}");
        // Dense labels follow g.nodes_vec() (sorted), so position == node id here.
        assert_ne!(best.labels[0], best.labels[3], "triangles not split");
        let comms: HashSet<i32> = best.labels.iter().copied().collect();
        assert_eq!(comms.len(), 2, "communities = {comms:?}");
    }
}
