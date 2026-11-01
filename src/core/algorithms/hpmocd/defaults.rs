//! Default parameters for HP-MOCD (Santos et al. 2025).
//! This Source Code Form is subject to the terms of The GNU General Public License v3.0
//! Copyright 2025 - Guilherme Santos. If a copy of the MPL was not distributed with this
//! file, You can obtain one at https://www.gnu.org/licenses/gpl-3.0.html

pub const DEFAULT_DEBUG_LEVEL: i8 = 0;
pub const DEFAULT_POP_SIZE: usize = 100;
pub const DEFAULT_NUM_GENS: usize = 100;
// Crossover and mutation probabilities of the paper's experiments (Santos et al. 2025,
// Sec. 5.1 and Table 1): C_P = 0.8, M_P = 0.2.
pub const DEFAULT_CROSS_RATE: f64 = 0.8;
pub const DEFAULT_MUT_RATE: f64 = 0.2;
