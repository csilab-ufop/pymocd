//! The two power laws of LFR: node degrees (exponent `tau1`) and community sizes
//! (exponent `tau2`), both discrete and truncated to an integer range.
//!
//! networkx draws them with Devroye's Zipf rejection sampler and redraws any
//! value above the upper bound, which is exactly the law `p(k) ∝ k^-tau` on
//! `[lo, hi]`; inverting that law's CDF gives the same law in one draw. Only the
//! degree law for an `average_degree` departs from networkx, whose solver is
//! biased and whose integer minimum cannot hit most averages anyway (see
//! `law_for_average`).
//! This Source Code Form is subject to the terms of The GNU General Public License v3.0
//! Copyright 2026 - Guilherme Santos. If a copy of the MPL was not distributed with this
//! file, You can obtain one at https://www.gnu.org/licenses/gpl-3.0.html

use std::collections::VecDeque;

use rand::{Rng, RngExt};

use super::{LfrError, MATCH_AVERAGE, POWER_LAW};

/// Most products the exact size sampler spends on its table, one per
/// (remainder, size) pair: under a second of work. Reaching `hi`, past which a
/// prefix is accepted with a fair probability, costs about `width^2 / 2`, so
/// wider ranges than about 60,000 stop short and need more prefix draws.
const TABLE_BUDGET: u64 = 2_000_000_000;
/// The table stops once its values vary by at most this factor over `hi`
/// consecutive remainders; a prefix is then accepted with probability at least
/// `1 / FLAT`.
const FLAT: f64 = 2.0;
/// Most sizes the exact sampler draws for prefixes before giving up. A prefix
/// takes about `n / E[size]` draws, few when rejection has already failed,
/// since that takes a large `E[size]`.
const PREFIX_DRAWS: usize = 100_000_000;

/// The law `p(k) ∝ k^-tau` on the integers `[lo, hi]`, sampled by inverting its
/// CDF (the original LFR code's `lower_bound` on a cumulative table), or a
/// mixture of two such laws (see `mixture`).
#[derive(Clone, Debug)]
pub struct PowerLaw {
    lo: usize,
    pmf: Vec<f64>,
    cdf: Vec<f64>,
}

impl PowerLaw {
    /// `1 <= lo <= hi`.
    pub fn new(tau: f64, lo: usize, hi: usize) -> Self {
        Self::mixture(tau, lo, hi, 0.0)
    }

    /// `(1 - alpha)` times the law on `[lo, hi]` plus `alpha` times the law on
    /// `[lo + 1, hi]`, `0 <= alpha < 1`: still `p(k) ∝ k^-tau` above `lo`, with
    /// less mass at `lo`. Each law's weights are taken relative to its own
    /// lowest value, `lo^-tau` or `(lo + 1)^-tau`, so a steep law cannot
    /// underflow to an all-zero table.
    pub fn mixture(tau: f64, lo: usize, hi: usize, alpha: f64) -> Self {
        debug_assert!(1 <= lo && lo <= hi && (0.0..1.0).contains(&alpha));
        let mut pmf = normalized(tau, lo, hi);
        if alpha > 0.0 && lo < hi {
            let upper = normalized(tau, lo + 1, hi);
            pmf[0] *= 1.0 - alpha;
            for (p, q) in pmf[1..].iter_mut().zip(upper) {
                *p = (1.0 - alpha) * *p + alpha * q;
            }
        }
        let mut acc = 0.0;
        let cdf = pmf
            .iter()
            .map(|&p| {
                acc += p;
                acc
            })
            .collect::<Vec<f64>>();
        let mut law = Self { lo, pmf, cdf };
        // rounding must not leave a draw u in [acc, 1) without a value
        *law.cdf.last_mut().expect("lo <= hi") = 1.0;
        law
    }

    pub fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> usize {
        let u: f64 = rng.random();
        self.lo + self.cdf.partition_point(|&c| c < u)
    }

    /// `p(k)`; zero outside `[lo, hi]`.
    pub fn pmf(&self, k: usize) -> f64 {
        k.checked_sub(self.lo)
            .and_then(|i| self.pmf.get(i))
            .copied()
            .unwrap_or(0.0)
    }

    /// Mean and standard deviation.
    pub fn moments(&self) -> (f64, f64) {
        let mean: f64 = (self.lo..=self.hi()).map(|k| k as f64 * self.pmf(k)).sum();
        let var: f64 = (self.lo..=self.hi())
            .map(|k| (k as f64 - mean).powi(2) * self.pmf(k))
            .sum();
        (mean, var.max(0.0).sqrt())
    }

    pub const fn lo(&self) -> usize {
        self.lo
    }

    pub fn hi(&self) -> usize {
        self.lo + self.pmf.len() - 1
    }
}

/// `k^-tau / sum_{j=lo}^{hi} j^-tau` for `k` in `[lo, hi]`, each weight taken
/// relative to `lo^-tau`.
fn normalized(tau: f64, lo: usize, hi: usize) -> Vec<f64> {
    let base = (lo as f64).ln();
    let mut pmf: Vec<f64> = (lo..=hi)
        .map(|k| (-tau * ((k as f64).ln() - base)).exp())
        .collect();
    let total: f64 = pmf.iter().sum();
    for p in &mut pmf {
        *p /= total;
    }
    pmf
}

/// The degree law for an `average_degree`, as `(j, alpha)` for
/// `PowerLaw::mixture(tau, j, kmax, alpha)`: its mean is `avg` exactly.
///
/// `E(m)`, the mean of `p ∝ k^-tau` on `[m, kmax]`, rises strictly with `m`,
/// so an integer minimum degree hits only the averages `E(1), E(2), ...`, and
/// they are far apart at low averages (`tau = 2.5`, `kmax = 100`: `E(1) =
/// 1.80`, `E(2) = 4.15`, `E(3) = 6.46`). The request sits between `E(j) <= avg < E(j + 1)`,
/// and the mixture `(1 - alpha) law[j, kmax] + alpha law[j + 1, kmax]` with
/// `alpha = (E(j) - avg) / (E(j) - E(j + 1))` has mean `avg`. Above `j` it is
/// the power law with exponent `tau`; only the mass at `j` is reduced.
///
/// networkx 3.4.2 instead bisects a real minimum against a mean that divides a
/// truncated numerator by an untruncated Hurwitz zeta, then rounds it: for
/// `tau1 = 2.5`, `avg = 20`, `kmax = 50` it picks 14 (mean 22.71), and for
/// `tau1 = 2` it rejects `avg = 20` outright.
///
/// The request must lie in `[E(1), kmax]`, widened by `tol` on both sides, or
/// it fails with networkx's message; within `tol` outside, the nearer end
/// (`law[1, kmax]` or every degree `kmax`) is used. A negative `tol` narrows
/// the top to `kmax - |tol|` but not the bottom, which the first check sets
/// to `1 + |tol|` instead: the range is `[max(E(1), 1 + |tol|), kmax - |tol|]`.
/// The checks that need no table come first, so a request far above `kmax`
/// fails in constant time.
pub fn law_for_average(
    tau: f64,
    avg: f64,
    kmax: usize,
    tol: f64,
) -> Result<(usize, f64), LfrError> {
    // every law on [1, kmax] has its mean in [1, kmax]; NaN fails both tests
    if kmax < 1 || !(avg >= 1.0 - tol && avg <= kmax as f64 + tol) {
        return Err(LfrError::runtime(MATCH_AVERAGE));
    }
    if avg >= kmax as f64 {
        return Ok((kmax, 0.0));
    }
    // E(m) from the top down, with sums relative to m^-tau so that a steep
    // law cannot underflow: s0(m) = sum_{k>=m} (k/m)^-tau = 1 + r s0(m + 1)
    // and s1(m) = sum k (k/m)^-tau = m + r s1(m + 1), r = (1 + 1/m)^-tau
    let (mut s0, mut s1) = (1.0, kmax as f64);
    let mut above = kmax as f64;
    for m in (1..kmax).rev() {
        let r = (-tau * (1.0 / m as f64).ln_1p()).exp();
        s0 = 1.0 + r * s0;
        s1 = m as f64 + r * s1;
        let mean = s1 / s0;
        if mean <= avg {
            let alpha = if above > mean {
                ((avg - mean) / (above - mean)).clamp(0.0, 1.0)
            } else {
                0.0
            };
            // alpha = (E(j) - avg) / (E(j) - E(j + 1)), written with positive
            // differences; alpha = 1 is the law on [j + 1, kmax] itself
            return Ok(if alpha >= 1.0 {
                (m + 1, 0.0)
            } else {
                (m, alpha)
            });
        }
        above = mean;
    }
    // avg < E(1): within tol of it, or no law has this mean
    if avg >= above - tol {
        Ok((1, 0.0))
    } else {
        Err(LfrError::runtime(MATCH_AVERAGE))
    }
}

/// Whether every degree sequence of `n` draws from `law` has an odd sum: one
/// possible degree `d`, odd, and `n` odd. No draw is needed to know it.
pub fn always_odd(law: &PowerLaw, n: usize) -> bool {
    law.lo() == law.hi() && n % 2 == 1 && law.lo() % 2 == 1
}

/// `n` i.i.d. degrees, the whole sequence redrawn until its sum is even, at most
/// `max_iters` times (networkx's `_powerlaw_sequence` with its even-sum test).
pub fn degree_sequence<R: Rng + ?Sized>(
    law: &PowerLaw,
    n: usize,
    max_iters: usize,
    rng: &mut R,
) -> Result<Vec<u32>, LfrError> {
    if always_odd(law, n) {
        return Err(LfrError::runtime(POWER_LAW));
    }
    let mut deg = vec![0u32; n];
    for _ in 0..max_iters {
        let mut parity = 0;
        for d in &mut deg {
            *d = law.sample(rng) as u32;
            parity ^= *d & 1;
        }
        if parity == 0 {
            return Ok(deg);
        }
    }
    Err(LfrError::runtime(POWER_LAW))
}

/// Community sizes: i.i.d. draws from the truncated law, conditioned on
/// summing to exactly `n`.
///
/// First networkx's way: draw until the sum reaches `n` and keep the sequence
/// only if it lands on `n`, at most `max_iters` times. A try lands with
/// probability about `1 / E[size]`, so for communities of a few thousand nodes
/// every try can miss; then the same conditional law is sampled by
/// `exact_sizes`.
pub fn community_sizes<R: Rng + ?Sized>(
    tau: f64,
    lo: usize,
    hi: usize,
    n: usize,
    max_iters: usize,
    rng: &mut R,
) -> Result<Vec<usize>, LfrError> {
    // a size above n always overshoots, so cutting the law at n is the same
    // as networkx's rejection of it. With no size left, networkx rejects
    // every try (min_community above n) or, for bounds crossed by a default
    // (given and crossed is a ValueError), loops forever.
    let top = hi.min(n);
    if lo > top {
        return Err(LfrError::runtime(POWER_LAW));
    }
    let law = PowerLaw::new(tau, lo, top);
    let mut sizes = Vec::new();
    for _ in 0..max_iters {
        sizes.clear();
        let mut total = 0;
        while total < n {
            let s = law.sample(rng);
            sizes.push(s);
            total += s;
        }
        if total == n {
            return Ok(sizes);
        }
    }
    exact_sizes(&law, n, rng)
}

/// The size law conditioned on the sum, sampled exactly.
///
/// `h(r)`, the probability that the running sum of i.i.d. sizes ever equals
/// `r`, obeys `h(0) = 1` and `h(r) = sum_s p(s) h(r - s)`; it is the total
/// weight `prod p(s_i)` of the sequences summing to `r`. Given `r` nodes still
/// to cover, the next size is `s` with probability `p(s) h(r - s) / h(r)`,
/// which reproduces the sequence law `prod p(s_i)` over sequences summing to
/// `n` exactly.
///
/// Tabulating `h` up to `n` costs `n` times the width of the size range. Far
/// from zero `h` flattens towards `1 / E[size]`, so the table stops at the
/// first `t` past which `h` varies by at most `FLAT` over `hi` remainders.
/// Sizes are then drawn freely while more than `t` nodes remain, which stops
/// with `r` in `(t - hi, t]` remaining, and that prefix is kept with
/// probability `h(r) / max h` over the window: the weight of its completions,
/// which is what the conditional law gives it. A kept prefix is completed
/// from the table. `h(n) = 0`, or `h` zero over a whole window, means no
/// sequence exists.
fn exact_sizes<R: Rng + ?Sized>(
    law: &PowerLaw,
    n: usize,
    rng: &mut R,
) -> Result<Vec<usize>, LfrError> {
    let hi = law.hi();
    let h = renewal_table(law, n);
    let t = h.len() - 1;
    let mut sizes = Vec::new();
    if t == n {
        if h[n].is_nan() || h[n] <= 0.0 {
            return Err(LfrError::runtime(POWER_LAW));
        }
        complete(law, &h, n, &mut sizes, rng);
        return Ok(sizes);
    }
    let hmax = h[(t + 1).saturating_sub(hi)..]
        .iter()
        .copied()
        .fold(0.0, f64::max);
    if hmax.is_nan() || hmax <= 0.0 {
        return Err(LfrError::runtime(POWER_LAW));
    }
    let mut draws = 0;
    while draws < PREFIX_DRAWS {
        sizes.clear();
        let mut r = n;
        let mut over = false;
        while r > t {
            let s = law.sample(rng);
            draws += 1;
            // only when the table stopped short of hi, on its budget
            if s > r {
                over = true;
                break;
            }
            sizes.push(s);
            r -= s;
        }
        if !over && rng.random::<f64>() * hmax < h[r] {
            complete(law, &h, r, &mut sizes, rng);
            return Ok(sizes);
        }
    }
    Err(LfrError::runtime(POWER_LAW))
}

/// `h(0..=t)` for `exact_sizes`: `t = n`, or the first `t` at which `h` is
/// within a factor `FLAT` over `(t - hi, t]`, or where `TABLE_BUDGET` ran out.
fn renewal_table(law: &PowerLaw, n: usize) -> Vec<f64> {
    let (lo, hi) = (law.lo(), law.hi());
    // p(hi), ..., p(lo): h(r) is then a forward dot product of two slices
    let rev: Vec<f64> = (lo..=hi).rev().map(|s| law.pmf(s)).collect();
    let mut h = vec![1.0];
    let mut spent = 0u64;
    // indices of the window's running minimum and maximum
    let (mut mins, mut maxs) = (VecDeque::new(), VecDeque::new());
    for r in 1..=n {
        let v = if r < lo {
            0.0
        } else {
            let top = hi.min(r);
            spent += (top - lo + 1) as u64;
            dot(&rev[hi - top..], &h[r - top..=r - lo])
        };
        h.push(v);
        while mins.back().is_some_and(|&j| h[j] >= v) {
            mins.pop_back();
        }
        mins.push_back(r);
        while maxs.back().is_some_and(|&j| h[j] <= v) {
            maxs.pop_back();
        }
        maxs.push_back(r);
        if r >= hi {
            for q in [&mut mins, &mut maxs] {
                while q.front().is_some_and(|&j| j + hi <= r) {
                    q.pop_front();
                }
            }
            let (low, high) = (h[mins[0]], h[maxs[0]]);
            if low > 0.0 && high <= FLAT * low {
                break;
            }
        }
        if spent >= TABLE_BUDGET {
            break;
        }
    }
    h
}

/// `sum a_i b_i` over four running sums, which the compiler can vectorize; the
/// order of the additions is fixed, so the result is too.
fn dot(a: &[f64], b: &[f64]) -> f64 {
    let mut acc = [0.0f64; 4];
    let (mut ca, mut cb) = (a.chunks_exact(4), b.chunks_exact(4));
    for (x, y) in (&mut ca).zip(&mut cb) {
        for k in 0..4 {
            acc[k] += x[k] * y[k];
        }
    }
    let tail: f64 = ca
        .remainder()
        .iter()
        .zip(cb.remainder())
        .map(|(x, y)| x * y)
        .sum();
    (acc[0] + acc[1]) + (acc[2] + acc[3]) + tail
}

/// Appends sizes summing to `r`, drawn from the conditional law with the
/// table; `h(r) > 0`.
fn complete<R: Rng + ?Sized>(
    law: &PowerLaw,
    h: &[f64],
    mut r: usize,
    sizes: &mut Vec<usize>,
    rng: &mut R,
) {
    let (lo, hi) = (law.lo(), law.hi());
    while r > 0 {
        let top = hi.min(r);
        let total: f64 = (lo..=top).map(|s| law.pmf(s) * h[r - s]).sum();
        let mut u = rng.random::<f64>() * total;
        // the last size with positive weight absorbs any rounding residue
        let mut pick = 0;
        for s in lo..=top {
            let w = law.pmf(s) * h[r - s];
            if w > 0.0 {
                pick = s;
                if u < w {
                    break;
                }
                u -= w;
            }
        }
        debug_assert!(pick > 0, "h(r) > 0 implies a positive weight");
        sizes.push(pick);
        r -= pick;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    fn exact_mean(tau: f64, lo: usize, hi: usize) -> f64 {
        let law = PowerLaw::new(tau, lo, hi);
        (lo..=hi).map(|k| k as f64 * law.pmf(k)).sum()
    }

    #[test]
    fn average_law_brackets_the_request() {
        // brackets computed independently by the acceptance harness
        // (tpl_means_by_kmin): tau1 = 2.5 -> {11, 12}, tau1 = 3 -> {12, 13},
        // tau1 = 3 / avg 10 -> {6, 7}; tau1 = 2, where networkx raises, -> {10, 11}
        for &(tau, avg, j) in &[
            (2.5, 20.0, 11),
            (3.0, 20.0, 12),
            (3.0, 10.0, 6),
            (2.0, 20.0, 10),
        ] {
            let (lo, alpha) = law_for_average(tau, avg, 50, 1e-7).unwrap();
            assert_eq!(lo, j, "tau {tau}, avg {avg}");
            assert!(exact_mean(tau, j, 50) <= avg && avg < exact_mean(tau, j + 1, 50));
            assert!((0.0..1.0).contains(&alpha));
        }
    }

    #[test]
    fn average_law_has_the_requested_mean() {
        let cases = [
            (2.5, 20.0, 50),
            (2.5, 3.0, 100),
            (2.0, 7.3, 100),
            (2.5, 31.0, 50),
            (3.0, 4.2, 30),
            (1.5, 40.0, 999),
            (3.5, 5.0, 50),
            (2.0, 49.9, 50),
        ];
        for &(tau, avg, kmax) in &cases {
            let (j, alpha) = law_for_average(tau, avg, kmax, 1e-7).unwrap();
            // the owner's form: (1 - alpha) E(j) + alpha E(j + 1) = avg
            let mix =
                (1.0 - alpha) * exact_mean(tau, j, kmax) + alpha * exact_mean(tau, j + 1, kmax);
            assert!((mix - avg).abs() <= 1e-9 * avg, "{tau} {avg}: {mix}");
            // and the table that is sampled has that mean
            let law = PowerLaw::mixture(tau, j, kmax, alpha);
            let mean: f64 = (j..=kmax).map(|k| k as f64 * law.pmf(k)).sum();
            assert!((mean - avg).abs() <= 1e-9 * avg, "{tau} {avg}: {mean}");
            assert!((law.lo(), law.hi()) == (j, kmax));
            // above j, exactly the power law: p(k) / p(k + 1) = ((k + 1) / k)^tau
            for k in j + 1..kmax {
                let ratio = law.pmf(k) / law.pmf(k + 1);
                let want = ((k + 1) as f64 / k as f64).powf(tau);
                assert!((ratio / want - 1.0).abs() < 1e-9, "k {k}");
            }
            // at j, at most the power law's share
            let pure = PowerLaw::new(tau, j, kmax);
            assert!(law.pmf(j) / law.pmf(j + 1) <= pure.pmf(j) / pure.pmf(j + 1) * (1.0 + 1e-12));
        }
        // a request equal to E(j) is the plain law on [j, kmax]
        let (j, _) = law_for_average(3.0, 12.0, 50, 1e-7).unwrap();
        let (k, alpha) = law_for_average(3.0, exact_mean(3.0, j, 50), 50, 1e-7).unwrap();
        assert!(k == j && alpha < 1e-9, "{k} {alpha}");
    }

    #[test]
    fn steep_laws_do_not_underflow() {
        // with tau1 = 1e6 all the mass of law[m, kmax] sits at m, so E(m) = m
        assert_eq!(law_for_average(1e6, 5.0, 50, 1e-7), Ok((5, 0.0)));
        let (j, alpha) = law_for_average(1e6, 5.25, 50, 1e-7).unwrap();
        assert_eq!(j, 5);
        assert!((alpha - 0.25).abs() < 1e-12);
        let law = PowerLaw::mixture(1e6, 5, 50, alpha);
        assert!((law.pmf(5) - 0.75).abs() < 1e-12 && (law.pmf(6) - 0.25).abs() < 1e-12);
    }

    #[test]
    fn unreachable_average_is_an_error() {
        let e = Err(LfrError::runtime(MATCH_AVERAGE));
        assert_eq!(law_for_average(3.0, 60.0, 50, 1e-7), e);
        assert_eq!(law_for_average(3.0, 0.5, 50, 1e-7), e);
        assert_eq!(law_for_average(3.0, f64::NAN, 50, 1e-7), e);
        assert_eq!(law_for_average(3.0, 10.0, 50, f64::NAN), e);
        assert_eq!(law_for_average(3.0, 5.0, 0, 1e-7), e);
        // E(1) = 1.35 for tau1 = 3 on [1, 50]: 1.2 is below it, by more than tol
        assert!(exact_mean(3.0, 1, 50) > 1.3);
        assert_eq!(law_for_average(3.0, 1.2, 50, 1e-7), e);
        // within tol of an end, that end
        assert_eq!(law_for_average(3.0, 50.0, 50, 1e-7), Ok((50, 0.0)));
        assert_eq!(law_for_average(3.0, 50.0 + 5e-8, 50, 1e-7), Ok((50, 0.0)));
        let e1 = exact_mean(3.0, 1, 50);
        assert_eq!(law_for_average(3.0, e1 - 5e-8, 50, 1e-7), Ok((1, 0.0)));
        // a negative tol narrows the range
        assert_eq!(law_for_average(3.0, 49.5, 50, -1.0), e);
        assert!(law_for_average(3.0, 10.0, 50, -1.0).is_ok());
        // far above kmax: rejected before any O(kmax) work
        assert_eq!(law_for_average(2.5, 2f64.powi(40), (1 << 27) - 1, 1e-7), e);
    }

    #[test]
    fn sampler_matches_the_law() {
        for law in [
            PowerLaw::new(2.5, 3, 40),
            PowerLaw::mixture(2.5, 3, 40, 0.6),
        ] {
            let mut rng = StdRng::seed_from_u64(5);
            let draws = 200_000.0;
            let mut count = vec![0usize; 41];
            for _ in 0..200_000 {
                count[law.sample(&mut rng)] += 1;
            }
            assert!(count[..3].iter().all(|&c| c == 0));
            // every cell within 5 binomial standard deviations
            for (k, &c) in count.iter().enumerate().skip(3) {
                let p = law.pmf(k);
                let sd = (draws * p * (1.0 - p)).sqrt();
                let dev = (c as f64 - draws * p).abs();
                assert!(dev <= 5.0 * sd + 1.0, "k={k}: {c} vs {}", draws * p);
            }
        }
    }

    #[test]
    fn degree_sum_is_even() {
        let law = PowerLaw::new(2.0, 1, 30);
        let mut rng = StdRng::seed_from_u64(1);
        for _ in 0..20 {
            let deg = degree_sequence(&law, 101, 500, &mut rng).unwrap();
            assert_eq!(deg.iter().map(|&d| u64::from(d)).sum::<u64>() % 2, 0);
        }
    }

    #[test]
    fn sizes_sum_to_n_within_bounds() {
        let mut rng = StdRng::seed_from_u64(2);
        for _ in 0..20 {
            let s = community_sizes(1.5, 20, 100, 1000, 500, &mut rng).unwrap();
            assert_eq!(s.iter().sum::<usize>(), 1000);
            assert!(s.iter().all(|&x| (20..=100).contains(&x)));
        }
        // no 30s add up to 100, by either sampler
        assert!(community_sizes(1.5, 30, 30, 100, 500, &mut rng).is_err());
        // no size left: above n, or bounds crossed by a default
        let e = Err(LfrError::runtime(POWER_LAW));
        assert_eq!(community_sizes(1.5, 101, 120, 100, 500, &mut rng), e);
        assert_eq!(community_sizes(1.5, 30, 20, 100, 500, &mut rng), e);
    }

    #[test]
    fn one_odd_degree_and_odd_n_fail_at_once() {
        let mut rng = StdRng::seed_from_u64(1);
        let e = Err(LfrError::runtime(POWER_LAW));
        assert_eq!(
            degree_sequence(&PowerLaw::new(2.0, 3, 3), 11, 500, &mut rng),
            e
        );
        assert_eq!(
            degree_sequence(&PowerLaw::new(2.0, 3, 3), 12, 500, &mut rng),
            Ok(vec![3; 12])
        );
        assert_eq!(
            degree_sequence(&PowerLaw::new(2.0, 4, 4), 11, 0, &mut rng),
            e
        );
    }

    #[test]
    fn prefix_sampler_is_the_conditional_law() {
        // sizes 1..=3 summing to 40: the table flattens long before 40, so the
        // prefix path runs. Its law of the community count K must match the
        // exact one, W[r][k] = sum_s p(s) W[r - s][k - 1].
        let (n, law) = (40, PowerLaw::new(2.0, 1, 3));
        let t = renewal_table(&law, n).len() - 1;
        assert!(t < n, "prefix path: the table stops at {t}");
        let mut w = vec![vec![0.0f64; n + 1]; n + 1];
        w[0][0] = 1.0;
        for r in 1..=n {
            for k in 1..=r {
                w[r][k] = (1..=3.min(r)).map(|s| law.pmf(s) * w[r - s][k - 1]).sum();
            }
        }
        let total: f64 = w[n].iter().sum();
        let mut rng = StdRng::seed_from_u64(4);
        let draws = 40_000;
        let mut count = vec![0usize; n + 1];
        for _ in 0..draws {
            let sizes = exact_sizes(&law, n, &mut rng).unwrap();
            assert_eq!(sizes.iter().sum::<usize>(), n);
            count[sizes.len()] += 1;
        }
        for k in 1..=n {
            let p = w[n][k] / total;
            let sd = (draws as f64 * p * (1.0 - p)).sqrt();
            let dev = (count[k] as f64 - draws as f64 * p).abs();
            assert!(
                dev <= 5.0 * sd + 1.0,
                "K={k}: {} vs {}",
                count[k],
                draws as f64 * p
            );
        }
    }

    #[test]
    fn exact_sampler_covers_what_rejection_misses() {
        let mut rng = StdRng::seed_from_u64(3);
        // max_iters = 0 skips rejection: sizes come from the exact sampler
        for _ in 0..5 {
            let s = community_sizes(1.5, 500, 2000, 8000, 0, &mut rng).unwrap();
            assert_eq!(s.iter().sum::<usize>(), 8000);
            assert!(s.iter().all(|&x| (500..=2000).contains(&x)));
        }
        // wide ranges at large n stay cheap: the table stops where h flattens
        let s = community_sizes(1.5, 100, 1000, 1_000_000, 0, &mut rng).unwrap();
        assert_eq!(s.iter().sum::<usize>(), 1_000_000);
        // its sizes follow the law: mean of many [20, 100] sizes vs the law's
        // mean, which conditioning on the sum barely shifts
        let mut all = Vec::new();
        for _ in 0..200 {
            all.extend(community_sizes(1.5, 20, 100, 1000, 0, &mut rng).unwrap());
        }
        let mean = all.iter().sum::<usize>() as f64 / all.len() as f64;
        let want = exact_mean(1.5, 20, 100);
        assert!((mean - want).abs() < 2.0, "{mean} vs {want}");
    }
}
