//! The variation operators over the locus genome.
//! This Source Code Form is subject to the terms of The GNU General Public License v3.0
//! Copyright 2025 - Guilherme Santos. If a copy of the MPL was not distributed with this
//! file, You can obtain one at https://www.gnu.org/licenses/gpl-3.0.html
//!
//! The paper says only that "some genes" are mutated and that `pm` is "the ratio of
//! mutation". The operators below follow the authors' C++ implementation
//! (`MOGAForClustering`, `GA.cpp`), which the first author's publications page links as
//! the code of the paper.

use rand::RngExt;

use super::locus::{Genome, NodeIndex};

/// Number of genes one mutation redraws, as a function of the node count `n`
/// (`GA::setMaxCluster`): a 34-node graph gets 11, a 128-node graph 2, a
/// 4941-node graph 24.
pub const fn mutation_size(n: usize) -> usize {
    if n < 20 {
        n / 2
    } else if n < 50 {
        n / 3
    } else if n < 100 {
        n / 5
    } else if n < 1000 {
        n / 50
    } else if n < 10_000 {
        n / 200
    } else {
        n / 1000
    }
}

/// Shi's "uniform two-point crossover", applied in place to two genomes
/// (`GA::crossover2`): every position is swapped between them with
/// probability 1/2. Each allele stays at its own position, so two legal
/// genomes remain legal.
pub fn uniform_crossover(a: &mut Genome, b: &mut Genome, rng: &mut impl rand::Rng) {
    for (ga, gb) in a.iter_mut().zip(b.iter_mut()) {
        if rng.random_bool(0.5) {
            std::mem::swap(ga, gb);
        }
    }
}

/// `genes` times, a uniformly random position is redrawn uniformly from its
/// node's neighbors (`GA::mutate`). Positions may repeat, and a redraw may
/// return the current allele.
pub fn mutate(genome: &mut Genome, idx: &NodeIndex, genes: usize, rng: &mut impl rand::Rng) {
    if genome.is_empty() {
        return;
    }
    for _ in 0..genes {
        let i = rng.random_range(0..genome.len());
        let cands = &idx.neighbor_candidates[i];
        genome[i] = cands[rng.random_range(0..cands.len())];
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::core::graph::Graph;

    fn path(n: usize) -> Graph {
        let mut g = Graph::new();
        for i in 0..n - 1 {
            g.add_edge(i as i32, i as i32 + 1);
        }
        g.finalize();
        g
    }

    #[test]
    fn mutation_size_follows_the_authors_table() {
        let cases = [
            (2, 1),
            (19, 9),
            (34, 11),
            (49, 16),
            (62, 12),
            (99, 19),
            (100, 2),
            (128, 2),
            (999, 19),
            (1000, 5),
            (4941, 24),
            (9999, 49),
            (10_000, 10),
        ];
        for (n, k) in cases {
            assert_eq!(mutation_size(n), k, "n = {n}");
        }
    }

    #[test]
    fn mutation_touches_at_most_genes_positions_and_stays_legal() {
        let g = path(60);
        let idx = NodeIndex::build(&g);
        let mut rng = rand::rng();
        for _ in 0..200 {
            let before: Genome = (0..idx.n()).map(|i| idx.neighbor_candidates[i][0]).collect();
            let mut after = before.clone();
            mutate(&mut after, &idx, 3, &mut rng);
            let changed = before.iter().zip(&after).filter(|(a, b)| a != b).count();
            assert!(changed <= 3, "{changed} genes changed");
            for (i, &a) in after.iter().enumerate() {
                assert!(idx.neighbor_candidates[i].contains(&a), "gene {i} = {a}");
            }
        }
    }

    #[test]
    fn crossover_swaps_alleles_position_by_position() {
        let mut rng = rand::rng();
        let a0: Genome = (0..50).collect();
        let b0: Genome = (100..150).collect();
        let (mut a, mut b) = (a0.clone(), b0.clone());
        uniform_crossover(&mut a, &mut b, &mut rng);
        for i in 0..50 {
            let pair = (a[i], b[i]);
            assert!(pair == (a0[i], b0[i]) || pair == (b0[i], a0[i]), "position {i}");
        }
    }
}
