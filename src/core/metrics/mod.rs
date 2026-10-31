//! Ground-truth agreement metrics between two label vectors: NMI, AMI, ARI, F1-Score.
//! This Source Code Form is subject to the terms of The GNU General Public License v3.0
//! Copyright 2026 - Guilherme Santos. If a copy of the MPL was not distributed with this
//! file, You can obtain one at https://www.gnu.org/licenses/gpl-3.0.html

pub mod ami;
pub mod ari;
pub mod contingency;
pub mod f1;
pub mod mi;
pub mod modularity;
pub mod nmi;
mod sum;

pub use contingency::{Contingency, contingency, entropy};
pub use mi::mutual_info;

/// (NMI, AMI, ARI, F1) between two label vectors of equal length.
///
/// NMI and AMI both normalise by the ARITHMETIC mean of the two marginal
/// entropies, matching scikit-learn's default. Geometric normalisation would
/// return different numbers, and AMI is not NMI — every table in this repo
/// depends on which of the four columns it reads. Empty input scores
/// `(1, 1, 1, 1)`, like each metric on its own and like scikit-learn.
pub fn gt_metrics(a: &[i64], b: &[i64]) -> (f64, f64, f64, f64) {
    assert_eq!(a.len(), b.len());
    let ct = contingency(a, b);
    let hu = entropy(&ct.rows, ct.n);
    let hv = entropy(&ct.cols, ct.n);
    let mi = mutual_info(&ct);
    (
        nmi::nmi(&ct, mi, hu, hv),
        ami::ami(&ct),
        ari::ari(&ct),
        f1::f1(&ct),
    )
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identical_labelings_are_perfect() {
        let a = vec![0, 0, 1, 1, 2, 2];
        let (nmi, ami, ari, f1) = gt_metrics(&a, &a);
        assert!((nmi - 1.0).abs() < 1e-12);
        assert!((ami - 1.0).abs() < 1e-12);
        assert!((ari - 1.0).abs() < 1e-12);
        assert!((f1 - 1.0).abs() < 1e-12);
    }

    #[test]
    fn permuted_labels_are_perfect() {
        let a = vec![0, 0, 1, 1];
        let b = vec![7, 7, 3, 3];
        let (nmi, ami, ari, f1) = gt_metrics(&a, &b);
        assert!((nmi - 1.0).abs() < 1e-12);
        assert!((ami - 1.0).abs() < 1e-12);
        assert!((ari - 1.0).abs() < 1e-12);
        assert!((f1 - 1.0).abs() < 1e-12);
    }

    /// `(a, b, nmi, ami, ari)`: NMI, AMI and ARI from scikit-learn 1.7.2,
    /// `f(labels_true=b, labels_pred=a)` with the default arithmetic mean.
    #[allow(clippy::type_complexity)]
    const SKLEARN: [(&[i64], &[i64], f64, f64, f64); 5] = [
        (
            &[0, 0, 0, 1, 1, 1],
            &[0, 0, 1, 1, 2, 2],
            0.5158037429793889,
            0.2987924581708903,
            0.24242424242424243,
        ),
        (
            &[0, 1, 2, 0, 1, 2, 0, 1, 2, 0],
            &[0, 0, 0, 1, 1, 1, 2, 2, 2, 2],
            0.012714017928048455,
            -0.3511857745978148,
            -0.25,
        ),
        (
            &[7, 2, 3, 2, 7, 5, 6, 3, 1],
            &[7, 2, 0, 0, 1, 5, 4, 3, 7],
            0.7874992165205942,
            -0.07142857142858043, // exact: -1/14 = -0.07142857142857142
            -0.07142857142857142,
        ),
        (
            &[
                4, 0, 0, 3, 1, 2, 0, 1, 3, 1, 4, 3, 3, 4, 3, 0, 4, 3, 0, 1, 0, 4, 3, 4, 1, 3, 3, 3,
                0, 2, 3, 4, 3, 2, 2, 1, 0, 1, 0, 1,
            ],
            &[
                5, 3, 0, 3, 1, 4, 5, 0, 3, 3, 6, 2, 1, 1, 3, 4, 2, 3, 4, 2, 6, 1, 6, 6, 6, 5, 3, 4,
                1, 2, 3, 6, 5, 3, 4, 3, 2, 6, 0, 2,
            ],
            0.2178285303531986,
            0.0031137910003660387,
            0.034653465346534656,
        ),
        (
            &[
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 0, 1, 2, 3,
            ],
            &[
                0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                0, 4, 4, 5, 6,
            ],
            0.8056743592562237,
            0.780211503578995,
            0.824945295404814,
        ),
    ];

    #[test]
    fn matches_scikit_learn() {
        for (a, b, nmi, ami, ari) in SKLEARN {
            let got = gt_metrics(a, b);
            assert!((got.0 - nmi).abs() < 1e-14, "nmi {} vs {nmi}", got.0);
            assert!((got.1 - ami).abs() < 1e-12, "ami {} vs {ami}", got.1);
            assert_eq!(got.2, ari, "ari is the same double as sklearn's");
            // symmetric in the two arguments
            let swapped = gt_metrics(b, a);
            assert!((swapped.0 - got.0).abs() < 1e-15);
            assert!((swapped.1 - got.1).abs() < 1e-15);
            assert_eq!(swapped.2, got.2);
            assert_eq!(swapped.3, got.3);
        }
    }

    #[test]
    fn empty_and_single_node_score_one_like_each_metric() {
        // sklearn: nmi = ami = ari = 1.0 on empty and on one-sample labelings
        assert_eq!(gt_metrics(&[], &[]), (1.0, 1.0, 1.0, 1.0));
        assert_eq!(gt_metrics(&[3], &[8]), (1.0, 1.0, 1.0, 1.0));
    }

    #[test]
    fn karate_ground_truth_against_itself() {
        let club: Vec<i64> = crate::core::graph::karate::KARATE_CLUB
            .iter()
            .map(|&c| i64::from(c))
            .collect();
        let (nmi, ami, ari, f1) = gt_metrics(&club, &club);
        assert!((nmi - 1.0).abs() < 1e-12);
        assert!((ami - 1.0).abs() < 1e-12);
        assert!((ari - 1.0).abs() < 1e-12);
        assert!((f1 - 1.0).abs() < 1e-12);
    }
}
