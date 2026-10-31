//! Benchmark-graph generators: graphs with planted communities and their ground truth.
//! This Source Code Form is subject to the terms of The GNU General Public License v3.0
//! Copyright 2026 - Guilherme Santos. If a copy of the MPL was not distributed with this
//! file, You can obtain one at https://www.gnu.org/licenses/gpl-3.0.html

use std::ffi::CString;

use pyo3::exceptions::{PyRuntimeError, PyRuntimeWarning, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};
use pyo3_stub_gen::derive::gen_stub_pyfunction;

use crate::core::generators::lfr::{self, LfrError, LfrParams};

/// LFR benchmark graph with crisp communities, wired to the requested degrees
/// and mixing.
///
/// Takes the parameters, defaults and validation of networkx's
/// ``LFR_benchmark_graph`` (networkx 3.4.2) and builds undirected, unweighted
/// graphs whose communities do not overlap. Community sizes are drawn from a
/// power law with exponent ``tau2`` on ``[min_community, max_community]``,
/// conditioned on summing to ``n``. Degrees are drawn from a power law with
/// exponent ``tau1``: on ``[min_degree, max_degree]`` when ``min_degree`` is
/// given; for ``average_degree``, from the same law with less weight on its
/// lowest degree, so that the expected mean degree is ``average_degree``
/// exactly (see Args).
///
/// Every node gets exactly its drawn degree ``k``, of which
/// ``floor((1 - mu) * k)`` or ``ceil((1 - mu) * k)`` edges stay in its
/// community, rounded at random so that the expected number of internal edge
/// ends is ``(1 - mu)`` times the sum of the degrees. One exception: when a
/// community's internal degrees add up to an odd number and none of its
/// members can be rounded the other way (every ``(1 - mu) * k`` there is an
/// integer, for instance), one of its nodes gets one edge more or one fewer
/// inside; a later exchange of members may move that node to another
/// community. A node goes to a community larger than ``ceil((1 - mu) * k)``.
/// When no assignment exists that way, the fewest nodes that make one exist
/// may join a community that only holds ``floor((1 - mu) * k)`` (networkx
/// admits on ``round``, never below the floor), and those that rounded up
/// are rounded down there. For each of them another node that rounded down
/// is rounded up instead: a member of the same community when one can be,
/// otherwise a member of another community of the same size, otherwise any
/// node of the communities being placed (all of them, unless a repair is
/// placing a few again). When no node can be (communities too small for
/// their members' ceilings, as with every degree 2 at ``mu = 0.3`` and the
/// default ``min_community``), the internal total falls short and the
/// fraction of edges between communities ends above ``mu``; a
/// ``RuntimeWarning`` gives the realized fraction when the shortfall moves it
/// by more than 0.001. The graph has no self-loop and no repeated edge.
/// networkx's own wiring overshoots both the degrees and the mixing (networkx
/// issue #6809).
///
/// Some draws of degrees and community sizes cannot be wired this way, for
/// instance when a node's internal degree fits no community drawn, or when
/// one community holds more than half of all edge ends that leave a
/// community. Such a draw is rejected and replaced by a new one, up to
/// ``max_iters`` draws. The graphs returned are therefore conditioned on a
/// draw the generator could complete: where many draws are rejected, the
/// realized degrees lean towards draws without large hubs, and the community
/// sizes lean too. When ``average_degree`` is given and the mean degree of
/// the graph misses it by more than 4 standard errors of the degree law and
/// by more than 2%, a ``RuntimeWarning`` says so. A ``min_degree`` request
/// leans the same way but gets no such warning. The user guide measures
/// these effects.
///
/// The graph is generated without holding the GIL and without checking for
/// signals, so Ctrl-C takes effect only when the call returns.
///
/// ``min_degree``, ``max_degree``, ``min_community`` and ``max_community``
/// also take a float with an integer value, such as ``20.0``, and use that
/// integer. networkx takes floats for them too, except ``max_degree``
/// together with ``average_degree``, where it raises ``TypeError``. A
/// fractional or non-finite value raises ``ValueError``.
///
/// The same parameters and ``seed`` give the same graph on the same pymocd
/// build.
///
/// Args:
///     n: number of nodes; the nodes are ``0`` to ``n - 1``.
///     tau1: power-law exponent of the degree distribution; must be greater
///         than 1 and finite.
///     tau2: power-law exponent of the community-size distribution; must be
///         greater than 1 and finite.
///     mu: mixing parameter in [0, 1]: the fraction of each node's edges that
///         leave its community.
///     average_degree: requested mean degree. Give exactly one of
///         ``average_degree`` and ``min_degree``. With ``E(m)`` the mean of
///         the power law on ``[m, max_degree]`` (``max_degree`` capped at
///         ``n - 1``) and ``E(j) <= average_degree
///         < E(j + 1)``, degrees follow ``(1 - a)`` times the law on ``[j,
///         max_degree]`` plus ``a`` times the law on ``[j + 1, max_degree]``,
///         ``a = (E(j) - average_degree) / (E(j) - E(j + 1))``: the power law
///         above ``j``, with mean ``average_degree``.
///     min_degree: smallest degree a node can draw.
///     max_degree: largest degree a node can draw, in (0, n]; ``None`` means
///         ``n``. A degree is capped at ``n - 1``, the most a simple graph
///         allows.
///     min_community: smallest community size; ``None`` means the smallest
///         drawn degree.
///     max_community: largest community size; ``None`` means the largest
///         drawn degree.
///     tol: slack on ``average_degree``: a request up to ``tol`` outside the
///         reachable range, from ``E(1)`` up to the degree cap, is accepted and
///         uses the nearer end of the range. A negative ``tol`` refuses
///         requests above the cap minus ``abs(tol)`` or below ``1 + abs(tol)``;
///         ``E(1)`` itself is not moved. Not networkx's bisection tolerance:
///         the degree law is computed exactly here.
///     max_iters: draws of degrees, community sizes and assignment before
///         giving up, and, within a draw, attempts at a degree sequence with
///         an even sum and at a size sequence summing to ``n``, as in networkx.
///         Fewer draws are made once the repairs of the rejected ones have
///         used up the work of about 32 draws. That bounds the repair work
///         only: draws that fail before any repair (no placement) use none
///         of it, a draw that fails in the wiring has used only the repairs
///         it ran before, and ``max_iters`` alone bounds how many such draws
///         are made.
///     seed: integer seed in [0, 2**64); other integers raise ``OverflowError``.
///         ``None`` seeds from the operating system.
///
/// Returns:
///     ``(edges, ground_truth)``: ``edges`` lists every undirected edge once
///         as ``(u, v)`` with ``u < v``, in sorted order; ``ground_truth`` maps
///         every node to its community, numbered ``0`` to ``k - 1`` by first
///         appearance in node order.
///
/// Raises:
///     ValueError: an invalid parameter, with networkx's message for the
///         checks networkx makes.
///     RuntimeError: valid parameters for which no graph was found: at once
///         when every draw would fail the same way, otherwise when the draws
///         ran out (see ``max_iters``). The message is that of networkx's
///         ``ExceededMaxIterations`` for the stage that failed last, or "Could
///         not wire the inter-community edges; try decreasing max_community"
///         when that stage was the edges between communities. On some inputs
///         the outcome or the message differs from networkx's; the user guide
///         lists them.
///
/// Warns:
///     RuntimeWarning: at most once per cause and call. When nodes held to
///         their floor could not be compensated and that raises the fraction
///         of edges between communities by more than 0.001: the message gives
///         ``mu``, the realized fraction, the shift, and the remedy (a larger
///         ``min_community``, with ``max_community`` at least as large). When
///         ``average_degree`` was given and the mean
///         degree of the returned graph differs from it by more than 4
///         standard errors of the degree law (its standard deviation over
///         ``sqrt(n)``) and by more than 2%: the message gives both numbers
///         and the number of draws rejected.
#[gen_stub_pyfunction]
#[pyfunction]
#[pyo3(name = "lfr_graph", signature = (
    n,
    tau1,
    tau2,
    mu,
    average_degree = None,
    min_degree = None,
    max_degree = None,
    min_community = None,
    max_community = None,
    tol = 1e-7,
    max_iters = 500,
    seed = None,
))]
#[gen_stub(override_return_type(type_repr = "tuple[builtins.list[tuple[builtins.int, builtins.int]], builtins.dict[builtins.int, builtins.int]]", imports = ("builtins")))]
#[allow(clippy::too_many_arguments)]
pub fn lfr_graph_fn(
    py: Python<'_>,
    n: i64,
    tau1: f64,
    tau2: f64,
    mu: f64,
    average_degree: Option<f64>,
    #[gen_stub(override_type(type_repr = "typing.Optional[builtins.int]", imports = ("typing", "builtins")))]
    min_degree: Option<f64>,
    #[gen_stub(override_type(type_repr = "typing.Optional[builtins.int]", imports = ("typing", "builtins")))]
    max_degree: Option<f64>,
    #[gen_stub(override_type(type_repr = "typing.Optional[builtins.int]", imports = ("typing", "builtins")))]
    min_community: Option<f64>,
    #[gen_stub(override_type(type_repr = "typing.Optional[builtins.int]", imports = ("typing", "builtins")))]
    max_community: Option<f64>,
    tol: f64,
    max_iters: i64,
    seed: Option<u64>,
) -> PyResult<(Bound<'_, PyList>, Bound<'_, PyDict>)> {
    let params = LfrParams {
        n,
        tau1,
        tau2,
        mu,
        average_degree,
        min_degree,
        max_degree,
        min_community,
        max_community,
        tol,
        max_iters,
    };
    let graph = py
        .detach(|| lfr::generate(&params, seed))
        .map_err(|e| match e {
            LfrError::Value(m) => PyValueError::new_err(m),
            LfrError::Runtime(m) => PyRuntimeError::new_err(m),
        })?;
    for warning in graph.warnings {
        let message = CString::new(warning).map_err(|e| PyRuntimeError::new_err(e.to_string()))?;
        PyErr::warn(py, py.get_type::<PyRuntimeWarning>().as_any(), &message, 1)?;
    }
    let edges = PyList::new(py, graph.edges)?;
    let ground_truth = PyDict::new(py);
    for (v, c) in graph.communities.into_iter().enumerate() {
        ground_truth.set_item(v, c)?;
    }
    Ok((edges, ground_truth))
}
