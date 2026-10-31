//! Ground-truth agreement metrics between two node→community dicts.
//! This Source Code Form is subject to the terms of The GNU General Public License v3.0
//! Copyright 2026 - Guilherme Santos. If a copy of the MPL was not distributed with this
//! file, You can obtain one at https://www.gnu.org/licenses/gpl-3.0.html

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::{PyBool, PyDict, PyInt, PyType};
use pyo3_stub_gen::derive::gen_stub_pyfunction;

use crate::core::metrics::{ami, ari, contingency, entropy, f1, mutual_info, nmi};

/// Align two `{node: community}` dicts into label vectors.
///
/// Only nodes present in BOTH dicts are scored, silently: a partition that
/// omits nodes is compared on the intersection, not penalised for them. Nodes
/// are matched and labels told apart with Python's own hashing and equality,
/// so any hashable node id or label works, as with a dict key.
///
/// Raises `ValueError` when the intersection is empty but the dicts are not
/// both empty: nothing would be scored, and the usual cause is node ids that
/// differ in type (ground truth read from a file with `str` ids) or numbering.
/// Two empty dicts score 1.0, as sklearn scores empty labelings.
fn aligned_labels(
    partition: &Bound<'_, PyDict>,
    gt: &Bound<'_, PyDict>,
) -> PyResult<(Vec<i64>, Vec<i64>)> {
    let (a, b) = match aligned_int_labels(partition, gt)? {
        Some(aligned) => aligned,
        None => aligned_any_labels(partition, gt)?,
    };
    if a.is_empty() && !(partition.is_empty() && gt.is_empty()) {
        return Err(PyValueError::new_err(
            "partition and gt share no node; check that both use the same node ids \
             (the int 0 and the str '0' are different keys)",
        ));
    }
    Ok((a, b))
}

/// [`aligned_labels`] for any hashable node ids and labels, through Python
/// dicts: each label is replaced by a dense id, so equal labels (`1`, `1.0`,
/// `True`) get the same id.
fn aligned_any_labels(
    partition: &Bound<'_, PyDict>,
    gt: &Bound<'_, PyDict>,
) -> PyResult<(Vec<i64>, Vec<i64>)> {
    let py = partition.py();
    let ids_p = PyDict::new(py);
    let ids_g = PyDict::new(py);
    let mut a = Vec::with_capacity(partition.len().min(gt.len()));
    let mut b = Vec::with_capacity(a.capacity());
    for (node, c) in partition.iter() {
        if let Some(cg) = gt.get_item(&node)? {
            a.push(dense_id(&ids_p, &c)?);
            b.push(dense_id(&ids_g, &cg)?);
        }
    }
    Ok((a, b))
}

/// [`aligned_any_labels`] for the usual input, where every scored label is an
/// `int`, a `bool` or a NumPy integer scalar that fits `i32`, the types whose
/// hash and equality are those of their value: the labels are used as they
/// are, without the dense-id dicts. On 1e5 to 1e6 nodes this made `nmi` 1.0
/// to 1.5 times as fast with 1e3 communities per side and 1.6 to 3.3 times
/// with 3e5. Nodes are matched through `gt` itself, as in the generic path.
/// `None` as soon as a scored label is anything else, so the generic path
/// handles it.
fn aligned_int_labels(
    partition: &Bound<'_, PyDict>,
    gt: &Bound<'_, PyDict>,
) -> PyResult<Option<(Vec<i64>, Vec<i64>)>> {
    let mut ints = PlainInts::default();
    let mut a = Vec::with_capacity(partition.len().min(gt.len()));
    let mut b = Vec::with_capacity(a.capacity());
    for (node, c) in partition.iter() {
        if let Some(cg) = gt.get_item(&node)? {
            let (Some(c), Some(cg)) = (ints.get(&c), ints.get(&cg)) else {
                return Ok(None);
            };
            a.push(c.into());
            b.push(cg.into());
        }
    }
    Ok(Some((a, b)))
}

/// Reads the labels for [`aligned_int_labels`]. Only `int`, `bool` and
/// NumPy's own integer scalar types are admitted: an `int` subclass, or any
/// object with `__index__`, converts to an integer too but may hash and
/// compare differently from it, and would then share a community with labels
/// that a dict tells apart from it.
#[derive(Default)]
struct PlainInts<'py> {
    /// NumPy integer types already admitted.
    numpy: Vec<Bound<'py, PyType>>,
}

impl<'py> PlainInts<'py> {
    fn get(&mut self, obj: &Bound<'py, PyAny>) -> Option<i32> {
        let plain = obj.is_exact_instance_of::<PyInt>()
            || obj.is_exact_instance_of::<PyBool>()
            || self.is_numpy_int(obj);
        if plain { obj.extract().ok() } else { None }
    }

    fn is_numpy_int(&mut self, obj: &Bound<'py, PyAny>) -> bool {
        let ty = obj.get_type();
        if self.numpy.iter().any(|t| t.is(&ty)) {
            return true;
        }
        if is_numpy_integer_type(&ty) {
            self.numpy.push(ty);
            return true;
        }
        false
    }
}

/// True for NumPy's own integer scalar types (`numpy.int64`, `numpy.uint8`,
/// …), not for subclasses defined elsewhere. NumPy is looked up in
/// `sys.modules`, not imported: no NumPy scalar exists before NumPy is.
fn is_numpy_integer_type(ty: &Bound<'_, PyType>) -> bool {
    let check = || -> PyResult<bool> {
        let numpy = ty.py().import("sys")?.getattr("modules")?.get_item("numpy")?;
        Ok(ty.is_subclass(&numpy.getattr("integer")?)? && ty.getattr("__module__")?.eq("numpy")?)
    };
    check().unwrap_or(false)
}

/// The id of `label` in `ids`, adding it as `len(ids)` the first time it is seen.
fn dense_id(ids: &Bound<'_, PyDict>, label: &Bound<'_, PyAny>) -> PyResult<i64> {
    if let Some(id) = ids.get_item(label)? {
        return id.extract();
    }
    let id = ids.len() as i64;
    ids.set_item(label, id)?;
    Ok(id)
}

/// (NMI, AMI, ARI, F1) of a partition against ground truth (gt).
///
/// Only nodes present in both dicts are scored; a node missing from either is
/// ignored, not penalised. Node ids and community labels can be any hashable
/// values and are matched as dict keys are, so ``1``, ``1.0`` and ``True``
/// name the same community. ``-1`` is an ordinary label, so all ``-1`` nodes
/// count as one community.
///
/// The four values are those of ``nmi``, ``ami``, ``ari`` and ``f1``, bit for
/// bit, from one alignment of the dicts: NMI, AMI and ARI follow
/// scikit-learn's ``normalized_mutual_info_score``,
/// ``adjusted_mutual_info_score`` and ``adjusted_rand_score`` with their
/// default arithmetic-mean normalisation, and F1 is the pair-counting F1.
///
/// Args:
///     partition: ``dict[node, community]`` to score, e.g. a detector's output.
///     gt: ground-truth ``dict[node, community]``.
///
/// Returns:
///     ``(nmi, ami, ari, f1)`` as floats. ``(1.0, 1.0, 1.0, 1.0)`` when the two
///     dicts are the same partition up to renaming the communities, which
///     includes a single shared node and two empty dicts.
///
/// Raises:
///     ValueError: if the two dicts share no node and are not both empty,
///         which usually means that their node ids differ in type or numbering.
#[gen_stub_pyfunction]
#[pyfunction]
#[pyo3(name = "gt_metrics", signature = (partition, gt))]
pub fn gt_metrics_fn(
    partition: &Bound<'_, PyDict>,
    gt: &Bound<'_, PyDict>,
) -> PyResult<(f64, f64, f64, f64)> {
    let (a, b) = aligned_labels(partition, gt)?;
    Ok(crate::core::metrics::gt_metrics(&a, &b))
}

/// Normalised mutual information between two {node:community} dicts.
///
/// The mutual information divided by the arithmetic mean of the two
/// entropies, as scikit-learn's ``normalized_mutual_info_score`` computes it
/// by default. It is not corrected for chance.
///
/// Only nodes present in both dicts are scored; a node missing from either is
/// ignored, not penalised. Node ids and community labels can be any hashable
/// values and are matched as dict keys are, so ``1``, ``1.0`` and ``True``
/// name the same community. ``-1`` is an ordinary label, so all ``-1`` nodes
/// count as one community.
///
/// Args:
///     partition: ``dict[node, community]`` to score, e.g. a detector's output.
///     gt: ground-truth ``dict[node, community]``.
///
/// Returns:
///     NMI as a ``float`` in [0, 1]: ``1.0`` when the two dicts are the same
///     partition up to renaming the communities (which includes a single
///     shared node and two empty dicts), ``0.0`` when exactly one of them puts
///     every shared node in one community.
///
/// Raises:
///     ValueError: if the two dicts share no node and are not both empty,
///         which usually means that their node ids differ in type or numbering.
#[gen_stub_pyfunction]
#[pyfunction]
#[pyo3(name = "nmi", signature = (partition, gt))]
pub fn nmi_fn(partition: &Bound<'_, PyDict>, gt: &Bound<'_, PyDict>) -> PyResult<f64> {
    let (a, b) = aligned_labels(partition, gt)?;
    let ct = contingency(&a, &b);
    let mi = mutual_info(&ct);
    Ok(nmi::nmi(
        &ct,
        mi,
        entropy(&ct.rows, ct.n),
        entropy(&ct.cols, ct.n),
    ))
}

/// Adjusted mutual information between two {node:community} dicts.
///
/// The mutual information corrected for the agreement expected between random
/// labelings with the same community sizes (Vinh, Epps & Bailey 2010) and
/// normalised by the arithmetic mean of the two entropies, as scikit-learn's
/// ``adjusted_mutual_info_score`` computes it by default. Agreement at chance
/// level scores 0 and agreement worse than chance is negative. It is evaluated
/// in a form in which no sum cancels, so it stays accurate where
/// scikit-learn's rounding error is large, mainly when one or both dicts are
/// nearly all singletons; the two can then differ by more than 1e-10.
///
/// Only nodes present in both dicts are scored; a node missing from either is
/// ignored, not penalised. Node ids and community labels can be any hashable
/// values and are matched as dict keys are, so ``1``, ``1.0`` and ``True``
/// name the same community. ``-1`` is an ordinary label, so all ``-1`` nodes
/// count as one community.
///
/// Args:
///     partition: ``dict[node, community]`` to score, e.g. a detector's output.
///     gt: ground-truth ``dict[node, community]``.
///
/// Returns:
///     AMI as a ``float``: ``1.0`` when the two dicts are the same partition up
///     to renaming the communities (which includes a single shared node and
///     two empty dicts), ``0.0`` when exactly one of them puts every shared
///     node in one community or every shared node in a community of its own.
///
/// Raises:
///     ValueError: if the two dicts share no node and are not both empty,
///         which usually means that their node ids differ in type or numbering.
#[gen_stub_pyfunction]
#[pyfunction]
#[pyo3(name = "ami", signature = (partition, gt))]
pub fn ami_fn(partition: &Bound<'_, PyDict>, gt: &Bound<'_, PyDict>) -> PyResult<f64> {
    let (a, b) = aligned_labels(partition, gt)?;
    Ok(ami::ami(&contingency(&a, &b)))
}

/// Adjusted Rand index between two {node:community} dicts.
///
/// The share of node pairs on which the two dicts agree (together in both or
/// apart in both), corrected for chance (Hubert & Arabie 1985). The pair counts
/// are exact integers and the final ratio is rounded the way scikit-learn's
/// ``adjusted_rand_score`` rounds it, so the two return the same float.
///
/// Only nodes present in both dicts are scored; a node missing from either is
/// ignored, not penalised. Node ids and community labels can be any hashable
/// values and are matched as dict keys are, so ``1``, ``1.0`` and ``True``
/// name the same community. ``-1`` is an ordinary label, so all ``-1`` nodes
/// count as one community.
///
/// Args:
///     partition: ``dict[node, community]`` to score, e.g. a detector's output.
///     gt: ground-truth ``dict[node, community]``.
///
/// Returns:
///     ARI as a ``float``; ``1.0`` when no pair of shared nodes is together in
///     one dict and apart in the other, which includes a single shared node
///     and two empty dicts.
///
/// Raises:
///     ValueError: if the two dicts share no node and are not both empty,
///         which usually means that their node ids differ in type or numbering.
#[gen_stub_pyfunction]
#[pyfunction]
#[pyo3(name = "ari", signature = (partition, gt))]
pub fn ari_fn(partition: &Bound<'_, PyDict>, gt: &Bound<'_, PyDict>) -> PyResult<f64> {
    let (a, b) = aligned_labels(partition, gt)?;
    Ok(ari::ari(&contingency(&a, &b)))
}

/// Pairwise F1 between two {node:community} dicts.
///
/// A pair-counting score (harmonic mean of pair precision and pair recall over
/// same-community node pairs), not the community-matching F1 some papers
/// report. scikit-learn has no such function.
///
/// Only nodes present in both dicts are scored; a node missing from either is
/// ignored, not penalised. Node ids and community labels can be any hashable
/// values and are matched as dict keys are, so ``1``, ``1.0`` and ``True``
/// name the same community. ``-1`` is an ordinary label, so all ``-1`` nodes
/// count as one community.
///
/// Args:
///     partition: ``dict[node, community]`` to score, e.g. a detector's output.
///     gt: ground-truth ``dict[node, community]``.
///
/// Returns:
///     F1 as a ``float`` in [0, 1]; ``1.0`` when neither dict puts two shared
///     nodes in one community (a single shared node, two empty dicts, or every
///     shared node in a community of its own on both sides), where the ratio
///     is 0/0.
///
/// Raises:
///     ValueError: if the two dicts share no node and are not both empty,
///         which usually means that their node ids differ in type or numbering.
#[gen_stub_pyfunction]
#[pyfunction]
#[pyo3(name = "f1", signature = (partition, gt))]
pub fn f1_fn(partition: &Bound<'_, PyDict>, gt: &Bound<'_, PyDict>) -> PyResult<f64> {
    let (a, b) = aligned_labels(partition, gt)?;
    Ok(f1::f1(&contingency(&a, &b)))
}
