//! Adjusted mutual information (Vinh, Epps & Bailey 2010).
//! This Source Code Form is subject to the terms of The GNU General Public License v3.0
//! Copyright 2026 - Guilherme Santos. If a copy of the MPL was not distributed with this
//! file, You can obtain one at https://www.gnu.org/licenses/gpl-3.0.html

use rustc_hash::FxHashMap;

use super::Contingency;
use super::contingency::ln_ratio;
use super::sum::Neumaier;

/// Hypergeometric weights below this fraction of the mode's weight are not
/// visited. The pmf is log-concave, so past the cutoff the weights fall at
/// least geometrically, and the skipped mass is below `1e-40 · (n + 2)` of
/// the total.
const TAIL: f64 = 1e-40;

/// `(size, number of communities of that size)`, ascending by size.
fn size_counts(marg: &FxHashMap<i64, f64>) -> Vec<(f64, f64)> {
    let mut sizes: Vec<f64> = marg.values().copied().collect();
    sizes.sort_by(f64::total_cmp);
    let mut out: Vec<(f64, f64)> = Vec::new();
    for s in sizes {
        match out.last_mut() {
            Some((size, count)) if *size == s => *count += 1.0,
            _ => out.push((s, 1.0)),
        }
    }
    out
}

/// `k ln(k / sqrt(a b))` for an overlap `k` of communities of sizes `a`, `b`:
/// the summand of both sums in [`ami`]. It is `≤ 0` since `k ≤ min(a, b)`,
/// and 0 for `k = 0` and for `k = a = b`.
fn overlap_log(k: f64, a: f64, b: f64) -> f64 {
    if k == 0.0 {
        return 0.0;
    }
    0.5 * k * ln_ratio(k as u128 * k as u128, a as u128 * b as u128)
}

/// `E[overlap_log(k, a, b)]` where `k ~ Hypergeometric(n, a, b)` is the
/// overlap of a community of size `a` with a random one of size `b`.
///
/// The weights come from the exact ratio of consecutive hypergeometric
/// probabilities, walked outwards from the mode and normalised by their sum.
/// The log-factorial form `ln a! + ln b! + … − ln n!` cancels terms of size
/// `n ln n` and leaves an error of about `ulp(ln n!)` in every log-weight,
/// ~2e-9 at n = 1e6; the ratio walk has none of that cancellation.
fn expected_overlap_log(a: f64, b: f64, n: f64) -> f64 {
    let lo = (a + b - n).max(0.0);
    let hi = a.min(b);
    let mode = ((a + 1.0) * (b + 1.0) / (n + 2.0)).floor().clamp(lo, hi);
    let rest = n - a - b;
    let mut mass = Neumaier::default();
    let mut sum = Neumaier::default();
    let mut visit = |k: f64, w: f64| {
        mass.add(w);
        sum.add(w * overlap_log(k, a, b));
    };
    visit(mode, 1.0);
    let (mut k, mut w) = (mode, 1.0);
    while k < hi {
        // P(k+1) / P(k)
        w *= (a - k) * (b - k) / ((k + 1.0) * (rest + k + 1.0));
        k += 1.0;
        if w < TAIL {
            break;
        }
        visit(k, w);
    }
    let (mut k, mut w) = (mode, 1.0);
    while k > lo {
        // P(k-1) / P(k)
        w *= k * (rest + k) / ((a - k + 1.0) * (b - k + 1.0));
        k -= 1.0;
        if w < TAIL {
            break;
        }
        visit(k, w);
    }
    sum.value() / mass.value()
}

/// `Σ_cells n_ij ln(n_ij / sqrt(a_i b_j))`, which is `n (I − H̄)` with `H̄` the
/// mean entropy.
fn observed(ct: &Contingency) -> f64 {
    let mut sum = Neumaier::default();
    for (&(i, j), &nij) in &ct.cells {
        sum.add(overlap_log(nij, ct.rows[&i], ct.cols[&j]));
    }
    sum.value()
}

/// The expectation of [`observed`] under the permutation model, which is
/// `n (E[I] − H̄)` (Vinh, Epps & Bailey 2010, eq. 24a, shifted by `H̄`).
///
/// The sum over every (row, column) pair depends on the two communities only
/// through their sizes, so it runs over pairs of DISTINCT sizes weighted by
/// how many communities have each. A labeling of n nodes has fewer than
/// `sqrt(2n)` distinct sizes, which turns the `R·C` pairs of the textbook
/// loop (n² for all-singleton labelings) into at most `2n`.
fn expected(ct: &Contingency) -> f64 {
    let rows = size_counts(&ct.rows);
    let cols = size_counts(&ct.cols);
    // a fixed outer side, so that swapping the arguments adds the same terms
    // in the same order
    let (outer, inner) = if rows <= cols {
        (rows, cols)
    } else {
        (cols, rows)
    };
    let mut sum = Neumaier::default();
    for &(a, count_a) in &outer {
        for &(b, count_b) in &inner {
            sum.add(count_a * count_b * expected_overlap_log(a, b, ct.n));
        }
    }
    sum.value()
}

/// Adjusted mutual information, normalised by the ARITHMETIC mean of the two
/// marginal entropies (scikit-learn's default). Two labelings that agree only
/// by chance score 0; agreement worse than chance is negative.
///
/// The textbook `(I − E[I]) / (H̄ − E[I])` subtracts quantities of size
/// `ln n` and is ill-conditioned when both labelings are nearly all
/// singletons, where `H̄ − E[I]` is O(1/n). Since the expected overlap of row
/// `i` with all columns is `a_i` (and likewise for columns), `H̄`, `I` and
/// `E[I]` share the term `ln n` and their differences reduce to
/// `n (I − H̄) = Y = Σ_cells n_ij ln(n_ij / sqrt(a_i b_j))` and
/// `n (E[I] − H̄) = Z = E[Y]`. Then `AMI = (Z − Y) / Z`. Both `Y` and `Z`
/// are sums of terms `≤ 0`, so neither sum cancels, and `Z < 0` outside the
/// special cases below.
///
/// Special cases, which sklearn's `adjusted_mutual_info_score` shares except
/// where noted:
/// - same partition (including both constant, and empty input): 1. When both
///   sides are all singletons, MI, EMI and the mean entropy are all `ln n`,
///   so the formula is 0/0; sklearn reaches 1 there through bit-identical
///   rounding, the other cases up to rounding.
/// - exactly one side constant: 0.
/// - exactly one side all singletons: 0, because then `Y = Z` exactly (MI and
///   EMI are both the entropy of the other side). Not special-cased in
///   sklearn, whose computed numerator is rounding noise here.
///
/// sklearn also moves its numerator and denominator at least machine epsilon
/// away from zero. Here the denominator `−Z` is a sum of positive terms and
/// cannot vanish, and doing so to the numerator would add an error of up to
/// `eps / (H̄ − E[I])`: 3.2e-10 for one pair on each side of 2e6 nodes.
pub fn ami(ct: &Contingency) -> f64 {
    if ct.same_partition() {
        return 1.0;
    }
    if ct.rows.len() == 1 || ct.cols.len() == 1 {
        return 0.0;
    }
    if ct.rows.len() as f64 == ct.n || ct.cols.len() as f64 == ct.n {
        return 0.0;
    }
    let (y, z) = (observed(ct), expected(ct));
    (z - y) / z
}

#[cfg(test)]
mod tests {
    use super::super::{contingency, entropy};
    use super::{ami, expected};

    fn score(a: &[i64], b: &[i64]) -> f64 {
        ami(&contingency(a, b))
    }

    #[test]
    fn permuted_labels_are_perfect() {
        assert!((score(&[0, 0, 1, 1, 2, 2, 3, 3], &[1, 1, 0, 0, 3, 3, 2, 2]) - 1.0).abs() < 1e-12);
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
    fn identical_singletons_are_perfect() {
        // the 0/0 case: sklearn 1.7.2 returns 1.0 for every n; the old code
        // returned 0.0 at n = 2 and 4 and 1.0071 at n = 309
        for n in (2..=50).chain([309, 1000]) {
            let a: Vec<i64> = (0..n).collect();
            let b: Vec<i64> = (0..n).rev().collect();
            assert_eq!(score(&a, &b), 1.0, "n = {n}");
        }
    }

    #[test]
    fn one_side_all_singletons_is_exactly_zero() {
        // EMI equals MI exactly here; sklearn 1.7.2 returns -5.8e-15 for this pair
        assert_eq!(score(&[0, 1, 2, 3, 4, 5], &[0, 1, 2, 3, 2, 4]), 0.0);
        assert_eq!(score(&[0, 1, 1, 2, 2, 2], &[5, 4, 3, 2, 1, 0]), 0.0);
    }

    #[test]
    fn nearly_all_singletons_keeps_full_precision() {
        // n = 1e5 singletons, plus pairs: {0, 1} in a; {0, 2} in b, then
        // {0, 1} and {5, 6} in c. H̄ − E[I] is ~1e-5 here, and computing
        // (I − E[I]) / (H̄ − E[I]) as written was 5.6e-11 and 5.3e-11 off.
        // References: mpmath, 60 digits.
        let n = 100_000;
        let a: Vec<i64> = (0..n).map(|v| if v == 1 { 0 } else { v }).collect();
        let b: Vec<i64> = (0..n).map(|v| if v == 2 { 0 } else { v }).collect();
        let c: Vec<i64> = (0..n)
            .map(|v| match v {
                1 => 0,
                6 => 5,
                _ => v,
            })
            .collect();
        for (x, y, exact) in [(&a, &b, -2.00002000060001e-10), (&a, &c, 0.6666666665777768)] {
            let got = score(x, y);
            assert!((got - exact).abs() < 1e-14, "{got} vs {exact}");
        }
    }

    /// Vinh et al. eq. 24a as the triple loop over every (row, column) pair
    /// and every overlap `k`, with log-factorials from a running sum of logs.
    fn textbook_emi(a: &[i64], b: &[i64]) -> f64 {
        let ct = contingency(a, b);
        let n = ct.n as usize;
        let mut lf = vec![0.0_f64; n + 1];
        for k in 1..=n {
            lf[k] = lf[k - 1] + (k as f64).ln();
        }
        let mut emi = 0.0;
        for &ai in ct.rows.values() {
            for &bj in ct.cols.values() {
                let (ai, bj) = (ai as usize, bj as usize);
                for k in (ai + bj).saturating_sub(n).max(1)..=ai.min(bj) {
                    let ln_w = lf[ai] + lf[bj] + lf[n - ai] + lf[n - bj]
                        - lf[n]
                        - lf[k]
                        - lf[ai - k]
                        - lf[bj - k]
                        - lf[n + k - ai - bj];
                    let (k, nf) = (k as f64, n as f64);
                    emi += k / nf * (nf * k / (ai as f64 * bj as f64)).ln() * ln_w.exp();
                }
            }
        }
        emi
    }

    #[test]
    fn grouped_expectation_matches_the_textbook_emi() {
        // expected() = n (E[I] − H̄), from distinct sizes and shifted by H̄
        let a: Vec<i64> = (0..200).map(|v| (v * 7 + v / 13) % 23).collect();
        let b: Vec<i64> = (0..200).map(|v| if v < 120 { v % 4 } else { v }).collect();
        let c: Vec<i64> = (0..200).map(|v| v / 50).collect();
        for (x, y) in [(&a, &b), (&b, &c), (&a, &c)] {
            let ct = contingency(x, y);
            let mean_h = 0.5 * (entropy(&ct.rows, ct.n) + entropy(&ct.cols, ct.n));
            let grouped = expected(&ct) / ct.n + mean_h;
            let textbook = textbook_emi(x, y);
            assert!(
                (grouped - textbook).abs() < 1e-12,
                "{grouped} vs {textbook}"
            );
        }
    }
}
