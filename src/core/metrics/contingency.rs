//! Contingency table between two label vectors + marginal entropy.
//! This Source Code Form is subject to the terms of The GNU General Public License v3.0
//! Copyright 2026 - Guilherme Santos. If a copy of the MPL was not distributed with this
//! file, You can obtain one at https://www.gnu.org/licenses/gpl-3.0.html

use rustc_hash::FxHashMap;

use super::sum::Neumaier;

/// Sparse joint count table of two labelings. Counts are stored as `f64`
/// because every consumer divides by `n` or takes a log; zero cells are absent.
pub struct Contingency {
    /// `(label_a, label_b)` -> number of nodes with that pair.
    pub cells: FxHashMap<(i64, i64), f64>,
    /// Marginal counts of the first labeling.
    pub rows: FxHashMap<i64, f64>,
    /// Marginal counts of the second labeling.
    pub cols: FxHashMap<i64, f64>,
    /// Number of nodes scored.
    pub n: f64,
}

/// Joint count table of two label vectors, which must be the same length and
/// aligned node by node.
pub fn contingency(a: &[i64], b: &[i64]) -> Contingency {
    let mut cells: FxHashMap<(i64, i64), f64> = FxHashMap::default();
    let mut rows: FxHashMap<i64, f64> = FxHashMap::default();
    let mut cols: FxHashMap<i64, f64> = FxHashMap::default();
    for (&x, &y) in a.iter().zip(b) {
        *cells.entry((x, y)).or_insert(0.0) += 1.0;
        *rows.entry(x).or_insert(0.0) += 1.0;
        *cols.entry(y).or_insert(0.0) += 1.0;
    }
    Contingency {
        cells,
        rows,
        cols,
        n: a.len() as f64,
    }
}

impl Contingency {
    /// True when the two labelings are the same partition up to renaming the
    /// communities: every row and every column holds exactly one non-empty
    /// cell. Also true for empty input.
    pub fn same_partition(&self) -> bool {
        self.cells.len() == self.rows.len() && self.cells.len() == self.cols.len()
    }
}

/// Shannon entropy of a marginal, in NATS (natural log). NMI is a ratio of
/// quantities from this function and `mutual_info`, so the base cancels — but
/// only as long as the two use the same one.
pub fn entropy(marg: &FxHashMap<i64, f64>, n: f64) -> f64 {
    let mut sum = Neumaier::default();
    for &c in marg.values().filter(|&&c| c > 0.0) {
        sum.add(c / n * ln_ratio(n as u128, c as u128));
    }
    sum.value()
}

/// `ln(p / q)` for positive integers `p`, `q` below 2^106 (counts or products
/// of two counts), to a few ulp. For `p/q ≥ 1/2` it is `ln_1p((p − q)/q)` with
/// `p − q` exact, so a ratio close to 1, such as a community holding almost
/// every node, keeps its low bits; `ln(p as f64 / q as f64)` would keep only
/// the rounding error of the quotient there. Below 1/2 the plain quotient is
/// accurate and `ln_1p` near −1 would not be.
pub(super) fn ln_ratio(p: u128, q: u128) -> f64 {
    if 2 * p >= q {
        ((p as i128 - q as i128) as f64 / q as f64).ln_1p()
    } else {
        (p as f64 / q as f64).ln()
    }
}

/// Same-community node pairs as exact integers: `(T, R, C)` with
/// `T = Σ C(n_ij, 2)` over cells, `R = Σ C(a_i, 2)` over rows and
/// `C = Σ C(b_j, 2)` over columns. ARI and F1 are ratios of these; summing
/// them in `f64` loses exactness once `n²/2` passes 2^53.
pub fn pair_counts(ct: &Contingency) -> (u128, u128, u128) {
    (
        pairs(ct.cells.values()),
        pairs(ct.rows.values()),
        pairs(ct.cols.values()),
    )
}

fn pairs<'a>(counts: impl Iterator<Item = &'a f64>) -> u128 {
    counts
        .map(|&c| {
            let c = c as u128;
            c * c.saturating_sub(1) / 2
        })
        .sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_and_marginals() {
        // table: [[2,1,0],[0,1,2]]
        let a = vec![0, 0, 0, 1, 1, 1];
        let b = vec![0, 0, 1, 1, 2, 2];
        let ct = contingency(&a, &b);
        assert_eq!(ct.n, 6.0);
        assert_eq!(ct.cells[&(0, 0)], 2.0);
        assert_eq!(ct.cells[&(0, 1)], 1.0);
        assert_eq!(ct.cells[&(1, 1)], 1.0);
        assert_eq!(ct.cells[&(1, 2)], 2.0);
        assert_eq!(ct.cells.len(), 4); // zero cells absent
        assert_eq!(ct.rows[&0], 3.0);
        assert_eq!(ct.rows[&1], 3.0);
        assert_eq!(ct.cols[&0], 2.0);
        assert_eq!(ct.cols[&1], 2.0);
        assert_eq!(ct.cols[&2], 2.0);
    }

    #[test]
    fn entropy_of_k_equal_clusters_is_ln_k() {
        let a = vec![0, 0, 0, 1, 1, 1];
        let b = vec![0, 0, 1, 1, 2, 2];
        let ct = contingency(&a, &b);
        assert!((entropy(&ct.rows, ct.n) - 2.0_f64.ln()).abs() < 1e-12);
        assert!((entropy(&ct.cols, ct.n) - 3.0_f64.ln()).abs() < 1e-12);
    }

    #[test]
    fn same_partition_ignores_label_names() {
        assert!(contingency(&[0, 0, 1, 2], &[5, 5, 9, 7]).same_partition());
        assert!(contingency(&[], &[]).same_partition());
        assert!(!contingency(&[0, 0, 1, 1], &[0, 0, 0, 1]).same_partition());
    }

    #[test]
    fn pair_counts_of_a_small_table() {
        // table [[2,1,0],[0,1,2]]: T = 1+0+0+1, R = 3+3, C = 1+1+1
        let ct = contingency(&[0, 0, 0, 1, 1, 1], &[0, 0, 1, 1, 2, 2]);
        assert_eq!(pair_counts(&ct), (2, 6, 3));
    }

    #[test]
    fn ln_ratio_keeps_the_low_bits_near_one() {
        // ln(1 + 1/n) = 1/n - 1/(2n²) + ...; the quotient (n+1)/n rounds to
        // within 1.1e-16 of 1 + 1/n, a relative error of 1e-10 at n = 1e6
        let n: u128 = 1_000_000;
        let exact = 1.0 / 1e6 - 0.5 / 1e12 + 1.0 / 3e18;
        assert!((ln_ratio(n + 1, n) - exact).abs() < 1e-15 * exact);
        assert!((ln_ratio(n, n + 1) + exact).abs() < 1e-15 * exact);
        assert_eq!(ln_ratio(7, 7), 0.0);
        assert!((ln_ratio(1, 1_000_000) + 1e6_f64.ln()).abs() < 1e-14);
        assert!((ln_ratio(1_000_000, 1) - 1e6_f64.ln()).abs() < 1e-14);
    }

    #[test]
    fn entropy_single_cluster_is_zero() {
        let a = vec![0, 0, 0, 0, 0, 0];
        let b = vec![0, 1, 2, 3, 4, 5];
        let ct = contingency(&a, &b);
        assert_eq!(entropy(&ct.rows, ct.n), 0.0);
        assert!((entropy(&ct.cols, ct.n) - 6.0_f64.ln()).abs() < 1e-12);
    }
}
