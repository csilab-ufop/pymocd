//! Normalised mutual information (arithmetic mean normalisation).
//! This Source Code Form is subject to the terms of The GNU General Public License v3.0
//! Copyright 2026 - Guilherme Santos. If a copy of the MPL was not distributed with this
//! file, You can obtain one at https://www.gnu.org/licenses/gpl-3.0.html

use super::Contingency;

/// Normalised mutual information, `mi / mean(hu, hv)`.
///
/// The ARITHMETIC mean is deliberate (scikit-learn's default); the geometric
/// mean is the other common convention and returns different numbers.
/// Special cases follow scikit-learn's `normalized_mutual_info_score`: 1 when
/// the two labelings are the same partition (which includes both constant and
/// empty input; sklearn reaches 1 there up to rounding), 0 when `mi` is 0.
pub fn nmi(ct: &Contingency, mi: f64, hu: f64, hv: f64) -> f64 {
    if ct.same_partition() {
        return 1.0;
    }
    if mi == 0.0 {
        return 0.0;
    }
    mi / (0.5 * (hu + hv))
}

#[cfg(test)]
mod tests {
    use super::super::{contingency, entropy, mutual_info};
    use super::nmi;

    fn score(a: &[i64], b: &[i64]) -> f64 {
        let ct = contingency(a, b);
        nmi(
            &ct,
            mutual_info(&ct),
            entropy(&ct.rows, ct.n),
            entropy(&ct.cols, ct.n),
        )
    }

    #[test]
    fn permuted_labels_are_perfect() {
        assert!((score(&[0, 0, 1, 1, 2, 2, 3, 3], &[1, 1, 0, 0, 3, 3, 2, 2]) - 1.0).abs() < 1e-12);
    }

    #[test]
    fn single_cluster_vs_singletons_is_zero() {
        // convention: one side constant, other side splits -> 0
        assert!(score(&[0, 0, 0, 0, 0, 0], &[0, 1, 2, 3, 4, 5]).abs() < 1e-12);
    }

    #[test]
    fn identical_labelings_are_exactly_one() {
        // mi and mean entropy round differently; this pair gave 1.0000000000000002
        let a: Vec<i64> = (0..8).map(|v| v % 3).collect();
        assert_eq!(score(&a, &a), 1.0);
    }

    #[test]
    fn both_constant_is_one() {
        assert!((score(&[0, 0, 0], &[5, 5, 5]) - 1.0).abs() < 1e-12);
    }
}
