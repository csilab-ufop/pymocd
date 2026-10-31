//! The Erdős–Gallai test, as a distance: how far a community's internal degree
//! sequence is from being realizable as a simple graph on its members.
//!
//! A pass/fail test is not enough for the repair that uses it. Exchanging one
//! member for an outsider rarely makes a bad sequence graphical in one step, but
//! it often shrinks the largest Erdős–Gallai violation, so the repair walks
//! that violation down to zero.
//!
//! The repair evaluates one exchange at a time, and a scan of the whole
//! community per evaluation made a community of 20,000 members cost seconds
//! per repair round. Only the first `h` inequalities can be violated, `h` the
//! number of `i` with `d_i >= i` (at most the square root of the degree sum),
//! so `Histogram` keeps a community's degrees by value and evaluates an
//! exchange in about `h` steps. Two cheaper tests decide most exchanges
//! without it: `zverovich` (sufficient) and `deficit` (exact, linear).
//! This Source Code Form is subject to the terms of The GNU General Public License v3.0
//! Copyright 2026 - Guilherme Santos. If a copy of the MPL was not distributed with this
//! file, You can obtain one at https://www.gnu.org/licenses/gpl-3.0.html

use std::collections::BTreeMap;

/// Buffers reused across calls; a call allocates only when a community is the
/// largest seen so far.
#[derive(Default)]
pub struct Scratch {
    count: Vec<u32>,
    sorted: Vec<u32>,
    prefix: Vec<i64>,
}

/// The Zverovich–Zverovich condition (Discrete Math. 1992): `count` positive
/// degrees with an even sum, all in `[low, high]`, are graphical when
/// `count >= (high + low + 1)^2 / (4 low)`. Sufficient only. A larger `high`
/// or a smaller `low` makes it harder to meet, so bounds on the true extremes
/// may be passed; with no positive degree the sequence is empty.
pub fn zverovich(count: u64, low: u32, high: u32) -> bool {
    if count == 0 {
        return true;
    }
    if low == 0 {
        return false;
    }
    let span = u128::from(high) + u128::from(low) + 1;
    4 * u128::from(low) * u128::from(count) >= span * span
}

/// A community's degrees by value, for the Erdős–Gallai deficit of a
/// sequence that differs from it in one exchanged member.
#[derive(Clone, Debug, Default)]
pub struct Histogram {
    count: BTreeMap<u32, u32>,
    len: u64,
    top: Vec<u32>,
    below: Vec<(u64, i64)>,
}

impl Histogram {
    pub fn new<I: IntoIterator<Item = u32>>(degrees: I) -> Self {
        let mut h = Self::default();
        for d in degrees {
            h.insert(d);
        }
        h
    }

    pub fn insert(&mut self, d: u32) {
        *self.count.entry(d).or_insert(0) += 1;
        self.len += 1;
    }

    /// Removes one degree `d`, which must be present.
    pub fn remove(&mut self, d: u32) {
        if let Some(c) = self.count.get_mut(&d) {
            *c -= 1;
            if *c == 0 {
                self.count.remove(&d);
            }
            self.len -= 1;
        }
    }

    /// Replaces one degree `out` by `inn`.
    pub fn replace(&mut self, out: u32, inn: u32) {
        self.remove(out);
        self.insert(inn);
    }

    /// The Erdős–Gallai deficit (see `deficit`) from the histogram, and the
    /// number of steps it took, about `h`, for the caller's work budget.
    ///
    /// With `d` sorted decreasingly, inequality `k` changes by
    /// `2 (d_k - (k - 1))` from inequality `k - 1`, so past `h` (where
    /// `d_k < k`) the violation never grows: the largest is among the first
    /// `h`. For `k <= h` every degree below `k` lies past position `k`, so the
    /// right-hand side is `k (k - 1) + k (#{d >= k} - k) + sum_{d < k} d`.
    pub fn deficit(&mut self) -> (u64, u64) {
        self.top.clear();
        'walk: for (&v, &c) in self.count.iter().rev() {
            for _ in 0..c {
                if (v as usize) < self.top.len() + 1 {
                    break 'walk;
                }
                self.top.push(v);
            }
        }
        let h = self.top.len();
        if h == 0 {
            return (0, 1);
        }
        // count and sum of the degrees below k, for k = 0..=h
        self.below.clear();
        self.below.resize(h + 1, (0, 0));
        for (&v, &c) in self.count.range(..h as u32) {
            let slot = &mut self.below[v as usize + 1];
            slot.0 += u64::from(c);
            slot.1 += i64::from(v) * i64::from(c);
        }
        for k in 1..=h {
            let (c, s) = self.below[k - 1];
            self.below[k].0 += c;
            self.below[k].1 += s;
        }
        let mut worst = 0i64;
        let mut lhs = 0i64;
        for k in 1..=h {
            lhs += i64::from(self.top[k - 1]);
            let (less, small) = self.below[k];
            let kk = k as i64;
            let at_least = (self.len - less) as i64;
            let rhs = kk * (kk - 1) + kk * (at_least - kk) + small;
            worst = worst.max(lhs - rhs);
        }
        (worst as u64, 2 * h as u64 + self.count.len().min(h) as u64)
    }

    /// The deficit with one degree `out` replaced by `inn`, leaving the
    /// histogram as it was.
    pub fn deficit_if(&mut self, out: u32, inn: u32) -> (u64, u64) {
        if out == inn {
            return self.deficit();
        }
        self.replace(out, inn);
        let result = self.deficit();
        self.replace(inn, out);
        result
    }
}

/// The largest Erdős–Gallai violation of a degree sequence:
///
/// `max(0, max_k [sum_{i<=k} d_i - k(k-1) - sum_{i>k} min(d_i, k)])`
///
/// over the sequence sorted in decreasing order. With an even sum, the sequence
/// is graphical exactly when this is zero. Linear in the length plus the largest
/// degree, which inside a community is below its size: a counting sort and one
/// pass over `k`, using that the `i > k` with `d_i >= k` form one block, of
/// length `#{d >= k} - k` when positive.
pub fn deficit<I: IntoIterator<Item = u32>>(degrees: I, s: &mut Scratch) -> u64 {
    s.sorted.clear();
    s.sorted.extend(degrees);
    let len = s.sorted.len();
    let Some(&max_d) = s.sorted.iter().max() else {
        return 0;
    };
    let max_d = max_d as usize;
    s.count.clear();
    s.count.resize(max_d + 1, 0);
    for &d in &s.sorted {
        s.count[d as usize] += 1;
    }
    let mut i = 0;
    for d in (0..=max_d).rev() {
        for _ in 0..s.count[d] {
            s.sorted[i] = d as u32;
            i += 1;
        }
    }
    s.prefix.clear();
    s.prefix.push(0);
    let mut acc = 0i64;
    for &d in &s.sorted {
        acc += i64::from(d);
        s.prefix.push(acc);
    }
    let total = acc;

    let mut worst = 0i64;
    // at_least = #{d >= k}, maintained as k rises
    let mut at_least = len;
    for k in 1..=len {
        at_least -= s.count.get(k - 1).copied().unwrap_or(0) as usize;
        let kk = k as i64;
        let block = at_least.saturating_sub(k) as i64;
        let rhs = kk * (kk - 1) + kk * block + (total - s.prefix[k.max(at_least)]);
        worst = worst.max(s.prefix[k] - rhs);
    }
    worst as u64
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::rngs::StdRng;
    use rand::{RngExt, SeedableRng};

    /// The textbook test, quadratic, for cross-checking.
    fn brute(seq: &[u32]) -> u64 {
        let mut d = seq.to_vec();
        d.sort_unstable_by(|a, b| b.cmp(a));
        let mut worst = 0i64;
        for k in 1..=d.len() {
            let lhs: i64 = d[..k].iter().map(|&x| i64::from(x)).sum();
            let rhs: i64 = (k * (k - 1)) as i64
                + d[k..]
                    .iter()
                    .map(|&x| i64::from(x.min(k as u32)))
                    .sum::<i64>();
            worst = worst.max(lhs - rhs);
        }
        worst as u64
    }

    #[test]
    fn known_sequences() {
        let mut s = Scratch::default();
        assert_eq!(
            deficit([2, 2, 0], &mut s),
            2,
            "two vertices cannot both have degree 2"
        );
        assert_eq!(deficit([3, 1, 1, 1], &mut s), 0, "star");
        assert_eq!(deficit([2, 2, 2], &mut s), 0, "triangle");
        assert_eq!(deficit([3, 3, 1, 1], &mut s), 2);
        assert_eq!(deficit([], &mut s), 0);
        // a few near-complete members and too few partners for them
        assert!(deficit([9, 9, 9, 8, 1, 1, 1, 1, 1, 0], &mut s) > 0);
    }

    #[test]
    fn matches_the_textbook_test() {
        let mut rng = StdRng::seed_from_u64(9);
        let mut s = Scratch::default();
        for _ in 0..5000 {
            let len = rng.random_range(1..30usize);
            let seq: Vec<u32> = (0..len).map(|_| rng.random_range(0..len as u32)).collect();
            assert_eq!(deficit(seq.iter().copied(), &mut s), brute(&seq), "{seq:?}");
        }
    }

    #[test]
    fn histogram_matches_the_textbook_test() {
        let mut rng = StdRng::seed_from_u64(10);
        for _ in 0..5000 {
            let len = rng.random_range(1..40usize);
            // skewed towards a few large degrees, as in a community with hubs
            let seq: Vec<u32> = (0..len)
                .map(|_| {
                    if rng.random_bool(0.2) {
                        rng.random_range(0..len as u32)
                    } else {
                        rng.random_range(0..4u32.min(len as u32))
                    }
                })
                .collect();
            let mut h = Histogram::new(seq.iter().copied());
            assert_eq!(h.deficit().0, brute(&seq), "{seq:?}");
            // one member exchanged, and the histogram left as it was
            let i = rng.random_range(0..len);
            let inn = rng.random_range(0..len as u32);
            let mut changed = seq.clone();
            changed[i] = inn;
            assert_eq!(
                h.deficit_if(seq[i], inn).0,
                brute(&changed),
                "{seq:?} [{i}] = {inn}"
            );
            assert_eq!(h.deficit().0, brute(&seq));
        }
    }

    #[test]
    fn zverovich_is_sufficient() {
        let mut rng = StdRng::seed_from_u64(11);
        let mut s = Scratch::default();
        let mut passed = 0;
        for _ in 0..20_000 {
            let len = rng.random_range(1..60usize);
            let lo = rng.random_range(1..5u32);
            let hi = rng.random_range(lo..lo + 12);
            let mut seq: Vec<u32> = (0..len).map(|_| rng.random_range(lo..=hi)).collect();
            if seq.iter().sum::<u32>() % 2 == 1 {
                seq[0] += 1;
            }
            let (low, high) = (*seq.iter().min().unwrap(), *seq.iter().max().unwrap());
            if zverovich(len as u64, low, high) {
                passed += 1;
                assert_eq!(deficit(seq.iter().copied(), &mut s), 0, "{seq:?}");
            }
        }
        assert!(passed > 1000, "the condition is met often enough to matter");
        assert!(zverovich(0, 0, 0) && !zverovich(3, 0, 2));
    }
}
