//! Parameter validation for the LFR generator: networkx's checks first, in its
//! order and with its messages, then the checks networkx lacks. Each of those
//! guards an input on which networkx 3.4.2 loops forever, raises an error that
//! does not name the parameter (`min() arg is an empty sequence`, `xmin < 1`),
//! or builds its graph from a value no power law is defined by: an infinite
//! exponent, a non-finite average or tolerance, a fractional degree or size
//! bound. The bound on `n` is this crate's (node ids are `u32`). Inputs on
//! which networkx stops with `ExceededMaxIterations` are left to the stages,
//! which raise the same error.
//! This Source Code Form is subject to the terms of The GNU General Public License v3.0
//! Copyright 2026 - Guilherme Santos. If a copy of the MPL was not distributed with this
//! file, You can obtain one at https://www.gnu.org/licenses/gpl-3.0.html

use super::LfrError;

pub const TAU1: &str = "tau1 must be greater than one";
pub const TAU2: &str = "tau2 must be greater than one";
pub const MU: &str = "mu must be in the interval [0, 1]";
pub const MAX_DEGREE: &str = "max_degree must be in the interval (0, n]";
pub const EXACTLY_ONE: &str = "Must assign exactly one of min_degree and average_degree";
pub const N_POSITIVE: &str = "n must be a positive integer";
pub const N_TOO_LARGE: &str = "n must be smaller than 2**32";
pub const MIN_DEGREE_POSITIVE: &str = "min_degree must be at least 1";
pub const MIN_DEGREE_ABOVE_MAX: &str = "min_degree must not exceed max_degree";
pub const MIN_COMMUNITY_POSITIVE: &str = "min_community must be at least 1";
pub const MAX_COMMUNITY_POSITIVE: &str = "max_community must be at least 1";
pub const COMMUNITY_BOUNDS: &str = "min_community must not exceed max_community";

/// `"<name> must be finite"`.
pub fn not_finite(name: &str) -> LfrError {
    LfrError::Value(format!("{name} must be finite"))
}

/// `"<name> must be an integer"`.
pub fn not_integer(name: &str) -> LfrError {
    LfrError::Value(format!("{name} must be an integer"))
}

/// A request in networkx's names and Python-facing types. Signed integers let
/// a negative value reach validation and fail there with a `ValueError`,
/// instead of failing the Python-to-Rust conversion with an `OverflowError`.
/// The degree and community bounds are floats because networkx takes floats
/// for them; one with an integer value means that integer.
#[derive(Clone, Debug)]
pub struct LfrParams {
    pub n: i64,
    pub tau1: f64,
    pub tau2: f64,
    pub mu: f64,
    pub average_degree: Option<f64>,
    pub min_degree: Option<f64>,
    pub max_degree: Option<f64>,
    pub min_community: Option<f64>,
    pub max_community: Option<f64>,
    pub tol: f64,
    pub max_iters: i64,
}

impl LfrParams {
    /// The four required parameters with networkx's defaults for the rest;
    /// one of `average_degree` and `min_degree` still has to be set.
    pub const fn new(n: i64, tau1: f64, tau2: f64, mu: f64) -> Self {
        Self {
            n,
            tau1,
            tau2,
            mu,
            average_degree: None,
            min_degree: None,
            max_degree: None,
            min_community: None,
            max_community: None,
            tol: 1e-7,
            max_iters: 500,
        }
    }
}

/// How the lower end of the degree law is fixed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum DegreeTarget {
    Min(usize),
    Average(f64),
}

/// A validated request in the types the stages use.
#[derive(Clone, Debug)]
pub struct Request {
    pub n: usize,
    pub tau1: f64,
    pub tau2: f64,
    pub mu: f64,
    pub degree: DegreeTarget,
    /// Highest degree a node can draw: `min(max_degree, n - 1)`. A simple graph
    /// has no node of degree `n`, so drawing from `[kmin, n - 1]` is networkx's
    /// default law `[kmin, n]` conditioned on what a graph can realize. Below
    /// `kmin` only for `min_degree = n`, which `generate` turns into an error.
    pub kmax: usize,
    pub min_community: Option<usize>,
    pub max_community: Option<usize>,
    pub tol: f64,
    pub max_iters: usize,
}

/// A count given as a float: finite, with an integer value. networkx draws
/// from a fractional bound as given (`int(xmin * u ** ...)` in `zipf_rv`
/// reaches below it), which is not the truncated power law.
fn whole(x: f64, name: &str) -> Result<f64, LfrError> {
    if !x.is_finite() {
        Err(not_finite(name))
    } else if x.fract() != 0.0 {
        Err(not_integer(name))
    } else {
        Ok(x)
    }
}

pub fn validate(p: &LfrParams) -> Result<Request, LfrError> {
    // networkx 3.4.2, generators/community.py:994-1009. NaN fails each test,
    // as it fails networkx's `not tau1 > 1`; an infinite max_degree fails
    // `0 < max_degree <= n`.
    if p.tau1.is_nan() || p.tau1 <= 1.0 {
        return Err(LfrError::value(TAU1));
    }
    if p.tau2.is_nan() || p.tau2 <= 1.0 {
        return Err(LfrError::value(TAU2));
    }
    if !(0.0..=1.0).contains(&p.mu) {
        return Err(LfrError::value(MU));
    }
    let n = p.n as f64;
    let max_degree = match p.max_degree {
        None => None,
        Some(d) if 0.0 < d && d <= n => Some(d),
        Some(_) => return Err(LfrError::value(MAX_DEGREE)),
    };
    if p.min_degree.is_some() == p.average_degree.is_some() {
        return Err(LfrError::value(EXACTLY_ONE));
    }

    // networkx takes an infinite exponent, and with average_degree samples a
    // law whose weights are NaN
    if p.tau1.is_infinite() {
        return Err(not_finite("tau1"));
    }
    if p.tau2.is_infinite() {
        return Err(not_finite("tau2"));
    }
    // networkx crashes on n = 0 (max() of an empty sequence, among others)
    if p.n < 1 {
        return Err(LfrError::value(N_POSITIVE));
    }
    // node ids are u32 throughout, and an edge key packs two of them in a u64
    if p.n > i64::from(u32::MAX) {
        return Err(LfrError::value(N_TOO_LARGE));
    }
    let max_degree = match max_degree {
        Some(d) => whole(d, "max_degree")? as usize,
        None => p.n as usize,
    };
    let degree = match (p.min_degree, p.average_degree) {
        (Some(k), _) => {
            // networkx raises "xmin < 1" from zipf_rv below 1, and loops
            // forever above max_degree. k = n = max_degree passes, as in
            // networkx, which then stops with ExceededMaxIterations or hangs.
            let k = whole(k, "min_degree")?;
            if k < 1.0 {
                return Err(LfrError::value(MIN_DEGREE_POSITIVE));
            }
            if k > max_degree as f64 {
                return Err(LfrError::value(MIN_DEGREE_ABOVE_MAX));
            }
            DegreeTarget::Min(k as usize)
        }
        // an unreachable average is networkx's "Could not match average_degree",
        // a RuntimeError, raised by the degree stage
        (None, Some(avg)) if avg.is_finite() => DegreeTarget::Average(avg),
        (None, Some(_)) => return Err(not_finite("average_degree")),
        (None, None) => unreachable!("exactly one is set"),
    };
    if !p.tol.is_finite() {
        return Err(not_finite("tol"));
    }
    // networkx raises "xmin < 1" below 1 and loops forever on crossed bounds;
    // a min_community above n is a size stage failure there, and here
    let min_community = p
        .min_community
        .map(|c| whole(c, "min_community"))
        .transpose()?;
    if min_community.is_some_and(|c| c < 1.0) {
        return Err(LfrError::value(MIN_COMMUNITY_POSITIVE));
    }
    let max_community = p
        .max_community
        .map(|c| whole(c, "max_community"))
        .transpose()?;
    if max_community.is_some_and(|c| c < 1.0) {
        return Err(LfrError::value(MAX_COMMUNITY_POSITIVE));
    }
    if let (Some(lo), Some(hi)) = (min_community, max_community)
        && lo > hi
    {
        return Err(LfrError::value(COMMUNITY_BOUNDS));
    }

    Ok(Request {
        n: p.n as usize,
        tau1: p.tau1,
        tau2: p.tau2,
        mu: p.mu,
        degree,
        kmax: max_degree.min(p.n as usize - 1),
        // a float beyond usize saturates, which the size stage reads as the
        // bound above n it is
        min_community: min_community.map(|c| c as usize),
        max_community: max_community.map(|c| c as usize),
        tol: p.tol,
        // networkx's range(max_iters) runs no attempt for a negative count
        max_iters: p.max_iters.max(0) as usize,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn ok() -> LfrParams {
        LfrParams {
            average_degree: Some(5.0),
            max_degree: Some(20.0),
            min_community: Some(20.0),
            max_community: Some(50.0),
            ..LfrParams::new(100, 3.0, 1.5, 0.1)
        }
    }

    fn err(p: &LfrParams) -> LfrError {
        validate(p).unwrap_err()
    }

    #[test]
    fn networkx_checks_and_messages() {
        let cases = [
            (LfrParams { tau1: 1.0, ..ok() }, TAU1),
            (
                LfrParams {
                    tau1: f64::NAN,
                    ..ok()
                },
                TAU1,
            ),
            (LfrParams { tau2: 0.5, ..ok() }, TAU2),
            (
                LfrParams {
                    tau2: f64::NAN,
                    ..ok()
                },
                TAU2,
            ),
            (LfrParams { mu: -0.01, ..ok() }, MU),
            (LfrParams { mu: 1.01, ..ok() }, MU),
            (
                LfrParams {
                    mu: f64::NAN,
                    ..ok()
                },
                MU,
            ),
            (
                LfrParams {
                    mu: f64::INFINITY,
                    ..ok()
                },
                MU,
            ),
            (
                LfrParams {
                    max_degree: Some(0.0),
                    ..ok()
                },
                MAX_DEGREE,
            ),
            (
                LfrParams {
                    max_degree: Some(-1.0),
                    ..ok()
                },
                MAX_DEGREE,
            ),
            (
                LfrParams {
                    max_degree: Some(101.0),
                    ..ok()
                },
                MAX_DEGREE,
            ),
            (
                LfrParams {
                    max_degree: Some(f64::NAN),
                    ..ok()
                },
                MAX_DEGREE,
            ),
            (
                LfrParams {
                    max_degree: Some(f64::INFINITY),
                    ..ok()
                },
                MAX_DEGREE,
            ),
            (
                LfrParams {
                    min_degree: Some(3.0),
                    ..ok()
                },
                EXACTLY_ONE,
            ),
            (
                LfrParams {
                    average_degree: None,
                    ..ok()
                },
                EXACTLY_ONE,
            ),
        ];
        for (p, msg) in cases {
            assert_eq!(err(&p), LfrError::value(msg), "{p:?}");
        }
    }

    #[test]
    fn networkx_order_is_kept() {
        // every argument is wrong: the tau1 message wins, as in networkx
        let p = LfrParams {
            tau1: 0.5,
            tau2: 0.5,
            mu: 2.0,
            max_degree: Some(0.0),
            ..ok()
        };
        assert_eq!(err(&p), LfrError::value(TAU1));
        // max_degree is checked before the min/average pair
        let p = LfrParams {
            max_degree: Some(0.0),
            average_degree: None,
            ..ok()
        };
        assert_eq!(err(&p), LfrError::value(MAX_DEGREE));
        // networkx's checks come before the finiteness of tau1
        let p = LfrParams {
            tau1: f64::INFINITY,
            mu: 2.0,
            ..ok()
        };
        assert_eq!(err(&p), LfrError::value(MU));
    }

    #[test]
    fn checks_networkx_lacks() {
        let md = |k| LfrParams {
            min_degree: Some(k),
            average_degree: None,
            ..ok()
        };
        let cases = [
            (
                LfrParams {
                    n: 0,
                    max_degree: None,
                    ..ok()
                },
                LfrError::value(N_POSITIVE),
            ),
            (md(0.0), LfrError::value(MIN_DEGREE_POSITIVE)),
            (md(21.0), LfrError::value(MIN_DEGREE_ABOVE_MAX)),
            (
                LfrParams {
                    min_community: Some(0.0),
                    ..ok()
                },
                LfrError::value(MIN_COMMUNITY_POSITIVE),
            ),
            (
                LfrParams {
                    max_community: Some(0.0),
                    ..ok()
                },
                LfrError::value(MAX_COMMUNITY_POSITIVE),
            ),
            (
                LfrParams {
                    min_community: Some(60.0),
                    ..ok()
                },
                LfrError::value(COMMUNITY_BOUNDS),
            ),
            (
                LfrParams {
                    n: 1 << 32,
                    max_degree: None,
                    ..ok()
                },
                LfrError::value(N_TOO_LARGE),
            ),
        ];
        for (p, e) in cases {
            assert_eq!(err(&p), e, "{p:?}");
        }
    }

    #[test]
    fn non_finite_and_fractional_values() {
        let md = |k| LfrParams {
            min_degree: Some(k),
            average_degree: None,
            ..ok()
        };
        let cases = [
            (
                LfrParams {
                    tau1: f64::INFINITY,
                    ..ok()
                },
                not_finite("tau1"),
            ),
            (
                LfrParams {
                    tau2: f64::INFINITY,
                    ..ok()
                },
                not_finite("tau2"),
            ),
            (
                LfrParams {
                    average_degree: Some(f64::NAN),
                    ..ok()
                },
                not_finite("average_degree"),
            ),
            (
                LfrParams {
                    average_degree: Some(f64::INFINITY),
                    ..ok()
                },
                not_finite("average_degree"),
            ),
            (
                LfrParams {
                    tol: f64::NAN,
                    ..ok()
                },
                not_finite("tol"),
            ),
            (
                LfrParams {
                    tol: f64::NEG_INFINITY,
                    ..md(5.0)
                },
                not_finite("tol"),
            ),
            (md(f64::NAN), not_finite("min_degree")),
            (md(f64::INFINITY), not_finite("min_degree")),
            (md(5.5), not_integer("min_degree")),
            (
                LfrParams {
                    max_degree: Some(19.5),
                    ..ok()
                },
                not_integer("max_degree"),
            ),
            (
                LfrParams {
                    min_community: Some(f64::NAN),
                    ..ok()
                },
                not_finite("min_community"),
            ),
            (
                LfrParams {
                    max_community: Some(f64::INFINITY),
                    ..ok()
                },
                not_finite("max_community"),
            ),
            (
                LfrParams {
                    max_community: Some(49.5),
                    ..ok()
                },
                not_integer("max_community"),
            ),
        ];
        for (p, e) in cases {
            assert_eq!(err(&p), e, "{p:?}");
        }
        // a float with an integer value is that integer
        let r = validate(&md(5.0)).unwrap();
        assert_eq!((r.degree, r.kmax), (DegreeTarget::Min(5), 20));
        assert_eq!((r.min_community, r.max_community), (Some(20), Some(50)));
    }

    #[test]
    fn degree_cap_and_defaults() {
        let r = validate(&LfrParams {
            max_degree: None,
            ..ok()
        })
        .unwrap();
        assert_eq!(r.kmax, 99, "max_degree = n is accepted and capped at n - 1");
        assert_eq!(r.degree, DegreeTarget::Average(5.0));
        // what networkx stops on with ExceededMaxIterations is not a ValueError
        let r = validate(&LfrParams {
            min_degree: Some(100.0),
            average_degree: None,
            max_degree: None,
            ..ok()
        })
        .unwrap();
        assert_eq!((r.degree, r.kmax), (DegreeTarget::Min(100), 99));
        let r = validate(&LfrParams {
            min_community: Some(101.0),
            max_community: Some(1e300),
            ..ok()
        })
        .unwrap();
        assert_eq!(r.min_community, Some(101));
        let r = validate(&LfrParams {
            max_iters: -3,
            ..ok()
        })
        .unwrap();
        assert_eq!(r.max_iters, 0);
    }
}
