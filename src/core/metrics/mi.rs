//! Mutual information of a contingency table (shared by NMI and AMI).
//! This Source Code Form is subject to the terms of The GNU General Public License v3.0
//! Copyright 2026 - Guilherme Santos. If a copy of the MPL was not distributed with this
//! file, You can obtain one at https://www.gnu.org/licenses/gpl-3.0.html

use super::Contingency;
use super::contingency::ln_ratio;
use super::sum::Neumaier;

/// Mutual information in NATS, matching [`super::entropy`]'s base.
///
/// As in scikit-learn's `mutual_info_score`, it is exactly 0 when either side
/// has a single community, and a sum that rounds below 0 is clipped to 0. The
/// log of `n·n_ij / (a_i·b_j)` goes through [`ln_ratio`]: when both labelings
/// are nearly one community, the big cell's ratio is within ~1/n of 1.
pub fn mutual_info(ct: &Contingency) -> f64 {
    if ct.rows.len() <= 1 || ct.cols.len() <= 1 {
        return 0.0;
    }
    let n = ct.n;
    let mut sum = Neumaier::default();
    for (&(i, j), &nij) in &ct.cells {
        let (a, b) = (ct.rows[&i] as u128, ct.cols[&j] as u128);
        sum.add((nij / n) * ln_ratio(n as u128 * nij as u128, a * b));
    }
    sum.value().max(0.0)
}

#[cfg(test)]
mod tests {
    use super::super::contingency::contingency;
    use super::*;

    fn mi(a: &[i64], b: &[i64]) -> f64 {
        mutual_info(&contingency(a, b))
    }

    #[test]
    fn identical_k_cluster_labels_give_ln_k() {
        let a = [0, 0, 1, 1];
        assert!((mi(&a, &a) - 2.0_f64.ln()).abs() < 1e-12);
    }

    #[test]
    fn permuted_labels_give_full_entropy() {
        let a = [0, 0, 1, 1, 2, 2, 3, 3];
        let b = [1, 1, 0, 0, 3, 3, 2, 2];
        assert!((mi(&a, &b) - 4.0_f64.ln()).abs() < 1e-12);
    }

    #[test]
    fn single_cluster_gives_exactly_zero() {
        assert_eq!(mi(&[0, 0, 0, 0], &[0, 1, 1, 2]), 0.0);
        assert_eq!(mi(&[], &[]), 0.0);
    }

    #[test]
    fn nearly_constant_labelings_keep_full_precision() {
        // n = 1e6; nodes 0..3 labelled 1 in a, nodes 2..4 labelled 7 in b
        // (node 2 shared). Reference from mpmath at 60 digits. Before
        // ln_ratio, NMI of this pair was 6.5e-12 off, relative.
        let n = 1_000_000;
        let a: Vec<i64> = (0..n).map(|v| i64::from(v < 3)).collect();
        let b: Vec<i64> = (0..n).map(|v| if (2..4).contains(&v) { 7 } else { 0 }).collect();
        let exact = 1.1519675191964778e-05;
        assert!((mi(&a, &b) - exact).abs() < 1e-14 * exact, "{}", mi(&a, &b));
    }

    #[test]
    fn independent_labelings_give_zero() {
        let a = [0, 1, 0, 1, 0, 1, 0, 1];
        let b = [0, 0, 1, 1, 2, 2, 3, 3];
        assert!(mi(&a, &b).abs() < 1e-12);
    }
}
