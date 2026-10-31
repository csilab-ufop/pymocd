//! Compensated summation shared by the metrics.
//! This Source Code Form is subject to the terms of The GNU General Public License v3.0
//! Copyright 2026 - Guilherme Santos. If a copy of the MPL was not distributed with this
//! file, You can obtain one at https://www.gnu.org/licenses/gpl-3.0.html

/// Neumaier's variant of Kahan summation (A. Neumaier, ZAMM 54(1):39-51, 1974).
/// The rounding error of every addition is collected in a second accumulator
/// and added back at the end, so the error does not grow with the number of
/// terms the way a plain left-to-right sum does. Entropy, MI and the expected
/// MI each add up to millions of similar terms on large inputs.
#[derive(Default, Clone, Copy)]
pub struct Neumaier {
    sum: f64,
    comp: f64,
}

impl Neumaier {
    pub fn add(&mut self, x: f64) {
        let t = self.sum + x;
        if self.sum.abs() >= x.abs() {
            self.comp += (self.sum - t) + x;
        } else {
            self.comp += (x - t) + self.sum;
        }
        self.sum = t;
    }

    pub fn value(self) -> f64 {
        self.sum + self.comp
    }
}

#[cfg(test)]
mod tests {
    use super::Neumaier;

    #[test]
    fn recovers_what_a_plain_sum_loses() {
        // 1 + 1e-16 rounds back to 1 in f64; a million of them must not vanish
        let mut s = Neumaier::default();
        s.add(1.0);
        for _ in 0..1_000_000 {
            s.add(1e-16);
        }
        assert!((s.value() - (1.0 + 1e-10)).abs() < 1e-15);
    }
}
