//! Rayon global thread-pool sizing, exposed to Python.
//! This Source Code Form is subject to the terms of The GNU General Public License v3.0
//! Copyright 2026 - Guilherme Santos. If a copy of the MPL was not distributed with this
//! file, You can obtain one at https://www.gnu.org/licenses/gpl-3.0.html

use pyo3::prelude::*;
use pyo3_stub_gen::derive::gen_stub_pyfunction;
use rayon::ThreadPoolBuilder;
use std::sync::Once;

static INIT_RAYON: Once = Once::new();

/// Size Rayon's global thread pool.
///
/// Only RIMPSO (``rimpso``, ``rimpso_fronts``, ``rimpso_select``), HP-MOCD
/// (``hpmocd``, ``hpmocd_fronts``, ``HpMocd``), ``cdrme`` and ``gdpso`` use the
/// pool; MMCoMO, CCM, KRM, Shi-MOCD and MOGA-Net are single-threaded.
///
/// The pool can only be built once per process, so call this before any of
/// those detectors runs. The first call builds the pool and prints two lines
/// to standard output; later calls do nothing. If a parallel detector has
/// already run, Rayon has built its default pool, and this call raises
/// ``pyo3_runtime.PanicException`` (a ``BaseException``), as does every later
/// call; the detectors keep using the default pool. Tests that want a
/// different thread count need a fresh interpreter.
///
/// Args:
///     num_threads: number of threads in the global pool.
///
/// Returns:
///     ``None``.
#[gen_stub_pyfunction]
#[pyfunction]
pub fn max_cores(num_threads: usize) -> PyResult<()> {
    INIT_RAYON.call_once(|| {
        ThreadPoolBuilder::new()
            .num_threads(num_threads)
            .build_global()
            .unwrap();
        debug!(
            warn,
            "Global thread pool initialized with {} threads", num_threads
        );
        debug!(
            warn,
            "Using max_cores again has no effect, due to static ThreadPoolBuilder initialization"
        );
    });
    Ok(())
}
