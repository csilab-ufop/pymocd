//! Adjusted Rand index (Hubert & Arabie 1985).
//! This Source Code Form is subject to the terms of The GNU General Public License v3.0
//! Copyright 2026 - Guilherme Santos. If a copy of the MPL was not distributed with this
//! file, You can obtain one at https://www.gnu.org/licenses/gpl-3.0.html

use super::Contingency;
use super::contingency::pair_counts;

/// Adjusted Rand index: the pair-counting agreement, corrected for chance.
///
/// With `T`, `R`, `C` from [`pair_counts`] and `P = C(n, 2)`, the textbook
/// `(T − RC/P) / ((R + C)/2 − RC/P)` equals `2(TP − RC) / ((R + C)P − 2RC)`.
/// Numerator and denominator are computed exactly in 128-bit integers (no
/// overflow below n = 2^32) and rounded once each, as scikit-learn's
/// `adjusted_rand_score` does with Python integers; the result is the same
/// double. In `f64` the `RC/P` term rounds and the subtraction then cancels
/// when both labelings are nearly one community.
///
/// 1.0 when the two labelings are the same partition as far as pairs go
/// (`T = R = C`, sklearn's `fn == fp == 0`), which includes n < 2.
pub fn ari(ct: &Contingency) -> f64 {
    let (t, r, c) = pair_counts(ct);
    if t == r && t == c {
        return 1.0;
    }
    let n = ct.n as u128;
    let p = n * (n - 1) / 2;
    let num = (t * p) as i128 - (r * c) as i128;
    let den = r * (p - c) + c * (p - r);
    2.0 * num as f64 / den as f64
}

#[cfg(test)]
mod tests {
    use super::super::contingency;
    use super::ari;

    fn score(a: &[i64], b: &[i64]) -> f64 {
        ari(&contingency(a, b))
    }

    #[test]
    fn hand_computed_value() {
        // index=2, row_pairs=6, col_pairs=3, expected=6*3/15=1.2, max=4.5
        // -> (2 - 1.2) / (4.5 - 1.2) = 8/33
        assert!((score(&[0, 0, 0, 1, 1, 1], &[0, 0, 1, 1, 2, 2]) - 8.0 / 33.0).abs() < 1e-12);
    }

    #[test]
    fn worse_than_chance_is_negative() {
        assert!(
            score(
                &[0, 1, 2, 0, 1, 2, 0, 1, 2, 0],
                &[0, 0, 0, 1, 1, 1, 2, 2, 2, 2]
            ) < 0.0
        );
        assert!(score(&[0, 1, 0, 1, 0, 1, 0, 1], &[0, 0, 1, 1, 2, 2, 3, 3]) < 0.0);
    }

    #[test]
    fn single_cluster_vs_singletons_is_zero() {
        assert!(score(&[0, 0, 0, 0, 0, 0], &[0, 1, 2, 3, 4, 5]).abs() < 1e-12);
    }

    #[test]
    fn fewer_than_two_nodes_is_one() {
        // sklearn 1.7.2: adjusted_rand_score([], []) = adjusted_rand_score([7], [0]) = 1.0
        assert_eq!(score(&[], &[]), 1.0);
        assert_eq!(score(&[0], &[7]), 1.0);
    }

    #[test]
    fn one_community_against_one_moved_node_is_exactly_zero() {
        // f64 pair arithmetic gave -2.07e-12 here; sklearn 1.7.2 and the exact value are 0
        let a = vec![0; 14375];
        let mut b = a.clone();
        b[0] = 1;
        assert_eq!(score(&a, &b), 0.0);
    }

    #[test]
    fn permuted_labels_are_perfect() {
        assert!((score(&[0, 0, 1, 1, 2, 2, 3, 3], &[1, 1, 0, 0, 3, 3, 2, 2]) - 1.0).abs() < 1e-12);
    }
}
