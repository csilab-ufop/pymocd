---
description: Release history of pymocd from v0.0.1 to v6.0.1, with the renamed, removed and re-pointed Python names in each release.
---

# Changelog

This page lists the changes in every tagged release of pymocd, newest first. It
starts with the work on `master` since the latest tag. Each entry is compiled
from the GitHub release notes, the commits between consecutive tags, and the
Python names that the extension module registers at each tag. Where an entry
says so, it also draws on the source distributions published on PyPI. Dates are
tag dates. Where the GitHub release was published more than a day later, the
entry says so. A version heading links to its GitHub release when there is one
and to the tagged source tree otherwise. "On PyPI" means that the
[PyPI release history](https://pypi.org/project/pymocd/) (its "Release history" tab) lists a
release with that version number. That release does not always contain the
tagged code; see v2.0.2 and v0.0.1.

## Upgrade notes

Some names were renamed, removed or re-pointed between releases. The table
below is derived from the module registration at each tag. "Replacement" is
filled in only where the release notes or the module source name one.

| Python name | Present in | Status in v6.0.1 | Replacement |
|---|---|---|---|
| `scale`, `scale_fronts` | v5.0.1 onward | Deprecated alias of `rimpso`, `rimpso_fronts` | [`rimpso`][pymocd.rimpso], [`rimpso_fronts`][pymocd.rimpso_fronts] |
| `smocc`, `smocc_fronts` | v5.1.0 | Removed in v6.0.0, no alias | [`rimpso`][pymocd.rimpso], [`rimpso_fronts`][pymocd.rimpso_fronts] |
| `mr_mocd`, `mr_mocd_fronts`, `mr_mocd_select` | v6.0.0 onward | Deprecated alias of the `rimpso*` functions | [`rimpso`][pymocd.rimpso], [`rimpso_fronts`][pymocd.rimpso_fronts], [`rimpso_select`][pymocd.rimpso_select] |
| `ariadne`, `ariadne_fronts` | v5.0.0 | Removed in v5.0.1 | `scale` or a baseline (v5.0.1 notes) |
| `Prism` class, `prism` | `Prism` v3.0.0 to v4.0.0, `prism` v4.0.0 | Removed in v5.0.0 | `scale` or a baseline (v5.0.1 notes) |
| `HpMocd` class | v0.0.2 to v5.0.0 (not registered at v1.0.1), v6.0.0 onward | Present | — |
| `Mocd` class | v5.0.0 | Removed in v5.0.1 | Plain functions such as [`mocd_q`][pymocd.mocd_q] (v5.0.1 notes) |
| `MOCD` class | v0.0.2 to v1.0.0, v1.1.3 to v1.2.5 | Removed in v1.3.0 (deprecated in v1.0.1) | — |
| `set_thread_count` | v1.0.0 to v5.0.0 | Renamed in v5.0.1 | [`max_cores`][pymocd.max_cores] |
| `fitness` | v0.0.1 to v5.0.0 | Removed in v5.0.1 | Native metrics such as [`gt_metrics`][pymocd.gt_metrics] (v5.0.1 notes) |
| `leoben` | v1.3.0 | Removed in v1.5.0 | — |
| `MMCoMO` class | v0.0.2 to v0.2.0 (an empty placeholder) | Removed in v1.0.0 | — |
| `nsga_ii`, `pesa_ii_maxq`, `pesa_ii_minimax` | v0.0.1 (module `pyevoea`) | Removed in v0.0.2 (algorithms became classes) | — |
| `pyevoea` (package and module name) | v0.0.1 | Renamed in v0.0.2 | `pymocd` |

Changes in behavior to watch for:

- **The name `scale` has pointed to different code over time.** In v5.0.1 it
  is SCALE. In v5.1.0 it is an alias of `smocc`, which the release notes
  describe as SCALE renamed, tuned and made deterministic. From v6.0.0 it is an
  alias of `rimpso`. The accepted keyword arguments differ in each of these
  versions.
- **`scale_fronts` changed its return shape in v6.0.0.** In v5.0.1 and v5.1.0
  it returns a list of partitions. From v6.0.0 it returns the 3-tuple
  `(partitions, points, selected)` of [`rimpso_fronts`][pymocd.rimpso_fronts].
- **`moga_net` changed its default `r` twice.** The default was `2.0` through
  v5.0.1, `1.5` in v5.1.0 and `2.0` again from v6.0.0. `moga_net_fronts` has
  existed since v5.1.0 and follows the same default.
- **The `HpMocd` class was removed in v5.0.1 and returned in v6.0.0.** The
  returning class is NSGA-II with pluggable Python objectives, according to the
  v6.0.0 notes.

## Unreleased

Compare: [v6.0.1...master](https://github.com/csilab-ufop/pymocd/compare/v6.0.1...master)

**Fixed**

- igraph input: every detector failed on an `igraph.Graph` because the node
  reader called `VertexSeq.__iter__()`, which igraph does not define. The
  reader now iterates the vertex sequence directly, and `tests/test_igraph.py`
  covers every detector.
- [`ami`][pymocd.ami] of two identical all-singleton labelings returned
  rounding noise instead of `1.0`: `0.0` for 2 and 4 nodes, `0.99982` for 14
  and `1.0071` for 309. Identical partitions now score exactly `1.0`.
- [`ami`][pymocd.ami] computed its log-factorials with an approximation
  accurate to about $2\times10^{-10}$ and summed the expected mutual
  information term by term. Its error grew with the number of nodes: for
  $10^4$ singletons against singletons plus one pair, whose exact AMI is 0, it
  returned $-2.0\times10^{-4}$. AMI is now evaluated in a form in which no sum
  cancels (see [Definitions](guide/evaluation.md#definitions)). Against 40- and
  60-digit reference values it was within $3\times10^{-16}$ on every input
  checked, up to $2\times10^6$ nodes, including labelings that are nearly all
  singletons on both sides, where scikit-learn is up to $10^{-4}$ off at $10^4$
  nodes (see
  [Agreement with scikit-learn](guide/evaluation.md#agreement-with-scikit-learn)).
- [`nmi`][pymocd.nmi] drifted from the exact value by up to $7\times10^{-11}$
  on labelings of $10^6$ to $10^7$ nodes because of plain summation, and could
  return `1.0000000000000002` for identical labelings. It now uses compensated
  summation and takes the logarithm of a ratio close to 1 from the exact
  integer difference; it was within $1.2\times10^{-16}$ of the exact value on
  those labelings and on nearly constant labelings of up to $2\times10^6$
  nodes.
- [`ari`][pymocd.ari] lost precision to cancellation when both labelings were
  nearly one community (an absolute error of $2\times10^{-12}$ at $1.4\times10^4$
  nodes and $2\times10^{-9}$ at $2.8\times10^7$), and returned `nan` for fewer
  than two shared nodes. It now uses exact integer pair counts, returns the
  same float as scikit-learn's `adjusted_rand_score`, and returns `1.0` for
  one shared node and for two empty dicts.
- [`gt_metrics`][pymocd.gt_metrics] returned `(0.0, 0.0, 0.0, 0.0)` when the
  two dicts share no node, while the single functions returned `1.0`, `1.0`,
  `nan` and `1.0`. Since nothing is scored then, all five now raise
  `ValueError`, unless both dicts are empty, which scores `1.0` as
  scikit-learn does for empty labelings. The usual cause is node ids of
  different types, such as `str` ids read from a file against a detector's
  `int` ids, which used to raise `TypeError` and would otherwise match no
  node now that any hashable id is accepted.

**Changed**

- The metric functions accept any hashable node ids and community labels,
  matched as dict keys are; integers outside the signed 32-bit range, strings
  and floats used to raise `OverflowError` or `TypeError`.
- [`ami`][pymocd.ami], and with it [`gt_metrics`][pymocd.gt_metrics], sums
  the expected mutual information over distinct community sizes instead of
  every pair of communities. Three inputs of $10^5$ and $10^6$ nodes that
  took 5.5 to 587 s now take under 0.1 s.
- NMI, AMI and ARI follow scikit-learn's special cases (for example, AMI is
  `0.0` when exactly one labeling is a single community), and AMI is exactly
  `0.0` when exactly one labeling is all singletons.
- The type stub `pymocd.pyi` now includes the `HpMocd` class and precise
  return types for `rimpso`, `rimpso_fronts`, `mmcomo`, `mmcomo_fronts` and
  `hpmocd_fronts`; the `hpmocd_fronts` docstring no longer says the `HpMocd`
  class is unreachable from Python.
- Documentation: the site was rebuilt with a user guide, one page per
  algorithm, an executed examples gallery and research pages.
- Build: the PyO3 requirement moves from 0.29.0 to 0.30.1, and `pyo3-stub-gen`
  0.23.2 is now required for PyO3 0.30.
- CI: workflows migrated from pip to uv. They then switched to `setup-python`
  for the interpreter and uv for installs.

## [v6.0.1](https://github.com/csilab-ufop/pymocd/tree/v6.0.1) — 2026-09-27

Compare: [v6.0.0...v6.0.1](https://github.com/csilab-ufop/pymocd/compare/v6.0.0...v6.0.1) · GitHub release: none · On PyPI

This is a maintenance release with no changes under `src/`. The Python API is
the same as in v6.0.0.

**Changed**

- Package metadata: the project URLs in `pyproject.toml` point to the
  `csilab-ufop` repository and to the documentation at csilab.ufop.br/pymocd.
- README and documentation: repository, clone and documentation-site URLs were
  updated, including the links in the algorithm pages. The old custom
  documentation domain was dropped.

## [v6.0.0](https://github.com/csilab-ufop/pymocd/releases/tag/v6.0.0) — 2026-09-08

Compare: [v5.1.0...v6.0.0](https://github.com/csilab-ufop/pymocd/compare/v5.1.0...v6.0.0) · GitHub release: yes · On PyPI

**Added**

- [`rimpso`][pymocd.rimpso]: RIMPSO. The release notes call it the "new
  recommended detector" and describe it as a "resolution-free CPM swarm". See
  [RIMPSO](algorithms/rimpso.md).
- [`rimpso_fronts`][pymocd.rimpso_fronts]: returns the whole archive, its
  points and the selected index.
- [`rimpso_select`][pymocd.rimpso_select]: runs the selector alone over
  external candidate partitions.
- [`cdrme`][pymocd.cdrme]: a re-implementation of CDRME, with no front
  accessor. A vendored copy of the CDRME reference implementation was added to
  the repository. See [CDRME](algorithms/cdrme.md).
- [`gdpso`][pymocd.gdpso]: a re-implementation of GDPSO. It optimizes modularity
  only and has no front accessor. See [GDPSO](algorithms/gdpso.md).
- [`HpMocd`][pymocd.HpMocd] class: NSGA-II with pluggable Python objectives,
  through [`set_objectives`][pymocd.HpMocd.set_objectives],
  [`set_on_generation`][pymocd.HpMocd.set_on_generation],
  [`generate_pareto_front`][pymocd.HpMocd.generate_pareto_front] and
  [`run`][pymocd.HpMocd.run].
- The aliases `mr_mocd`, `mr_mocd_fronts` and `mr_mocd_select` (see
  **Deprecated**).

**Changed**

- **Breaking:** `scale` now aliases `rimpso` instead of `smocc`.
- **Breaking:** `scale_fronts` now aliases `rimpso_fronts`. It returns a 3-tuple
  instead of a list of partitions.
- **Breaking:** `moga_net` and `moga_net_fronts` changed their default `r` from
  1.5 to 2.0.
- Internals: the CCM, HP-MOCD, KRM, MMCoMO, Shi-MOCD and MOGA-Net modules were
  each split into per-responsibility modules, and an algorithms module index was
  added. Clippy pedantic and nursery lints were fixed or silenced, and comments
  were trimmed across the graph, metrics, core and utility code.
- Documentation and README:
    - Every detector is documented together with its original implementation.
    - The assortative selector and the `seed` parameter are documented.
    - The node scope of the pair objective was corrected.
- Tests: the Python suite covers the new detectors.
- Benchmarks:
    - The hardened harness gained algorithm and size filters.
    - A Leiden-CPM baseline was added, also on the SNAP real networks.
    - MMCoMO is no longer capped on LFR.
    - Runs record ARI, HCV, `gt_k` and `mu_real` and save the partitions.
    - The scheduling of parallel and serial detectors was reordered.

**Removed**

- **Breaking:** `smocc` and `smocc_fronts`, with no alias. The release notes
  direct users to `rimpso` and `rimpso_fronts`.

**Deprecated**

- `mr_mocd`, `mr_mocd_fronts`, `mr_mocd_select`, `scale` and `scale_fronts`.
  Each is the same function object as the corresponding `rimpso*` function.
  They are not part of the type stubs.

The release notes list these names as unchanged:
[`hpmocd`][pymocd.hpmocd], [`hpmocd_fronts`][pymocd.hpmocd_fronts],
[`mocd_q`][pymocd.mocd_q], [`mocd_d`][pymocd.mocd_d], [`ccm`][pymocd.ccm],
[`ccm_fronts`][pymocd.ccm_fronts], [`krm`][pymocd.krm],
[`krm_fronts`][pymocd.krm_fronts], [`mmcomo`][pymocd.mmcomo],
[`mmcomo_fronts`][pymocd.mmcomo_fronts], [`max_cores`][pymocd.max_cores],
[`gt_metrics`][pymocd.gt_metrics], [`nmi`][pymocd.nmi], [`ami`][pymocd.ami],
[`ari`][pymocd.ari], [`f1`][pymocd.f1].

!!! note "Development history inside this cycle"
    Two detectors were developed during this cycle: MOPOTS and MR-MOCD.
    Commits removed SMOCC and MOPOTS as "superseded by mr_mocd". MR-MOCD was
    then renamed RIMPSO, and its front selector was replaced before the tag
    (commit `15ae58c`). In v6.0.0 the `mr_mocd*` names remain only as
    deprecated aliases of the `rimpso*` functions. MOPOTS is not part of any
    tagged release.

## [v5.1.0](https://github.com/csilab-ufop/pymocd/releases/tag/v5.1.0) — 2026-08-07

Compare: [v5.0.1...v5.1.0](https://github.com/csilab-ufop/pymocd/compare/v5.0.1...v5.1.0) · GitHub release: yes · On PyPI

The release notes call **SMOCC** (Sparse Multi-Objective Co-evolutionary
Community detection) the headline of this release. According to the notes, the
experimental `scale` detector was renamed, tuned, made deterministic and
optimized, and it now ships alongside HP-MOCD.

**Added**

- `smocc` and `smocc_fronts`, the renamed `scale` detector.
    - New keyword arguments `macro_cap` and `micro_mut`.
    - `smocc_fronts` also takes `refine` and `topo_mode`, as `scale_fronts`
      did in v5.0.1, and a new `obj_mode`.
- Pareto-front accessors [`ccm_fronts`][pymocd.ccm_fronts],
  [`krm_fronts`][pymocd.krm_fronts] and
  [`moga_net_fronts`][pymocd.moga_net_fronts]. The release notes describe them
  as matching `hpmocd_fronts` and `smocc_fronts`.
  [`hpmocd_fronts`][pymocd.hpmocd_fronts] is first registered at this tag.
- Correctness suites under `tests/`: determinism, partition validity and sanity
  smoke tests.
- Documentation: an automated MkDocs site with a generated API reference, pages
  for every detector's arguments and returns, examples pages and an SVG
  logo and favicon. The site was then served at `pymocd.guiolvr.com`, a domain
  dropped in v6.0.1.

**Changed**

- **Breaking:** `scale` and `scale_fronts` became aliases of `smocc` and
  `smocc_fronts`. The `adaptive_stop` and `conv_pval` arguments that `scale`
  took in v5.0.1 are not part of the `smocc` signature.
- **Breaking:** `moga_net` changed its default `r` from 2.0 to 1.5. The release
  notes call this "the paper's r = 1.5 setting". It also uses neighbour-only
  alleles. The default returned to 2.0 in v6.0.0.
- SMOCC defaults to a tuned configuration and a scalarisation selector. Commit
  subjects record these internal changes:
    - an O(n log n) two-objective non-dominated sort;
    - active-set label propagation in the decode step;
    - dense accumulators in the objectives;
    - no per-node hash maps in the micro operators;
    - component-split and wadj-weighted refinement.
- Baselines standardized (`mocd_q`/`mocd_d`, `moga_net`, `ccm`, `krm`,
  `mmcomo`):
    - Dense `Vec` label arrays replace hash-map genomes.
    - `FxHashMap` and `rayon` were dropped where, in the notes' words, "they
      bought nothing".
    - Shi-MOCD uses PESA-II region selection and an Erdős–Rényi null model.
    - MMCoMO uses deterministic tie-breaks. The notes state that it reproduces
      the paper's reported results.
- Build and CI:
    - Actions run on Node 24, and tests run on pull requests only.
    - Stub generation builds without the `extension-module` feature.
    - Rust doc-tests are disabled for the Python-facing docstrings.
- Benchmarks:
    - Moved out of `tests/` into per-experiment modules.
    - A hardened, resumable campaign harness was added, split into real and
      synthetic suites. The real suite includes SNAP `com-` networks with
      ground-truth communities.
    - The benchmark files from PR #54 were the first contribution by
      matheustrindade2002.

Some files were added and removed again within this cycle, so they are in
neither v5.0.1 nor v5.1.0: a benchmarks page in the new documentation site and
its CI job, a SCALE profiling binary, and a standalone SNAP benchmark script.

## [v5.0.1](https://github.com/csilab-ufop/pymocd/releases/tag/v5.0.1) — 2026-07-27

Compare: [v5.0.0...v5.0.1](https://github.com/csilab-ufop/pymocd/compare/v5.0.0...v5.0.1) · GitHub release: yes · On PyPI

The release notes cover every change since v4.0.0 (their compare link is
v4.0.0...v5.0.1). v5.0.0 has no release notes. According to the notes, every
detector returns a single partition as `dict[node, community]`, with isolated
nodes labelled `-1`.

**Added**

- `scale` and `scale_fronts`: SCALE, a sparse macro–micro co-evolutionary
  NSGA-II. Its objectives are KKM and ratio cut, and it selects a partition with
  a label-free SBM/MDL criterion. `adaptive_stop=True` enables statistically
  tested early stopping. `scale_fronts` returns a list of partitions.
- [`mmcomo`][pymocd.mmcomo] and [`mmcomo_fronts`][pymocd.mmcomo_fronts]:
  MMCoMO, a macro/micro co-evolutionary NSGA-II. The commit that made it a
  self-contained rewrite states that it reproduces Table IV of the paper. See
  [MMCoMO](algorithms/mmcomo.md).
- Native ground-truth metrics [`gt_metrics`][pymocd.gt_metrics],
  [`nmi`][pymocd.nmi], [`ami`][pymocd.ami], [`ari`][pymocd.ari] and
  [`f1`][pymocd.f1]. See [Evaluating partitions](guide/evaluation.md).
- The notes also list [`mocd_q`][pymocd.mocd_q], [`mocd_d`][pymocd.mocd_d],
  [`moga_net`][pymocd.moga_net], [`ccm`][pymocd.ccm] and [`krm`][pymocd.krm] as
  new. They were already exposed at the v5.0.0 tag.

**Changed**

- **Breaking:** `set_thread_count` was renamed
  [`max_cores`][pymocd.max_cores].
- HP-MOCD became self-contained, with its own NSGA-II, operators and
  objectives. Shi-MOCD, MOGA-Net and CCM each got their own objectives. The
  kernel ratio cut moved into KRM and the SBM/MDL code into SCALE.
- Build: the PyO3 requirement moved from 0.28.2 to 0.29.0, with
  `pyo3-stub-gen` 0.23 for compatibility. The generated `pymocd.pyi` is no
  longer tracked in git.
- The README and the benchmark scripts were updated to the current API.

**Removed**

- **Breaking:** `ariadne()` and `ariadne_fronts()`. The notes list PRISM and
  Ariadne as removed and suggest migrating to `scale()` or one of the
  baselines. `prism()` and `Prism` were already gone at the v5.0.0 tag.
- **Breaking:** the `Mocd` and `HpMocd` classes. In the notes' words, "all
  detectors are plain functions". `HpMocd` returned in v6.0.0.
- **Breaking:** `fitness`. The notes point to the new native metrics instead.
- The shared metaheuristics module.

## [v5.0.0](https://github.com/csilab-ufop/pymocd/tree/v5.0.0) — 2026-06-28

Compare: [v4.0.0...v5.0.0](https://github.com/csilab-ufop/pymocd/compare/v4.0.0...v5.0.0) · GitHub release: none · Not on PyPI

This version has no release notes. The v5.0.1 notes cover its changes.

**Added**

- Shared evolutionary engines (NSGA-II, NSGA-III with reference-point niching,
  and PESA-II) in a modular layout of algorithms and metaheuristics.
- Per-method crossover, mutation and generator operators, and shared
  objectives.
- Detector functions:
    - `ariadne` and `ariadne_fronts`, an auto-gamma NSGA-II detector with an
      SBM/MDL selector (removed in v5.0.1);
    - [`mocd_q`][pymocd.mocd_q] and [`mocd_d`][pymocd.mocd_d] (Shi-MOCD on
      PESA-II with the MOCD-Q and MOCD-D selectors);
    - [`moga_net`][pymocd.moga_net] (MOGA-Net);
    - [`ccm`][pymocd.ccm] (NSGA-III-CCM);
    - [`krm`][pymocd.krm] (NSGA-III-KRM).
- README: a detector table, then a rewrite with an algorithm table and
  references.

**Changed**

- [`hpmocd`][pymocd.hpmocd], a function since v4.0.0, was ported to the shared
  NSGA-II engine.

**Removed**

- **Breaking:** PRISM, meaning the `Prism` class and the `prism()` function.

At this tag the module still registers the `HpMocd` and `Mocd` classes,
`set_thread_count` and `fitness`.

## [v4.0.0](https://github.com/csilab-ufop/pymocd/releases/tag/v4.0.0) — 2026-05-03

Compare: [v3.0.0...v4.0.0](https://github.com/csilab-ufop/pymocd/compare/v3.0.0...v4.0.0) · GitHub release: yes · On PyPI

The release notes summarize this version as "Type stubs, public new `prism()`
algorithm (in beta testing for now) / `hpmocd()` functions, and crate
restructure". They state that there are no breaking changes to the existing
`Prism`/`HpMocd` Python API and that the new functions are additive.

**Added**

- Top-level functions `prism(graph)` and [`hpmocd(graph)`][pymocd.hpmocd]. They
  run each algorithm with library defaults, with no class instantiation. The
  `Prism` and `HpMocd` classes remain for tuning.
- Full type stubs (PEP 561). They are generated from the Rust annotations with
  `pyo3-stub-gen` by a `stub_gen` binary and a `make stubs` target, and CI
  generates them before building.
- Class-level docstrings: `help(Prism)` and `help(HpMocd)` show parameter
  tables.
- Usage examples.
- Node labels: the graph module gained a constructor that maps any hashable
  node label to an internal id (commit "support arbitrary hashable node
  labels"). At this tag only PRISM uses it, so `Prism` and `prism` accept
  string and negative-integer labels, while the HP-MOCD path still requires
  integer node ids. For the labels the current version accepts, see
  [Input graphs](guide/input-graphs.md).

**Changed**

- Python 3.9 was dropped from the release build matrix, which now builds for
  Python 3.10 to 3.13. The notes give problems with stub generation as the
  reason. `pyproject.toml` has declared `requires-python >= 3.10` since
  v1.0.1.
- Every Python-facing docstring was tightened and standardized.
- The crate was restructured. The algorithm modules moved under `src/core`, and
  PRISM was split into separate modules.
- PyPI metadata was enriched, the stubs ship with the package, and the
  dependencies were updated.
- Benchmarks: the real-world network benchmark is wired into the Makefile, and
  the HP-MOCD runs relabel non-integer nodes.

**Fixed**

- Example import paths after the move to `src/core`, and the README GIF.
- CI passes the Python interpreter to maturin and uses updated Python versions.

## [v3.0.0](https://github.com/csilab-ufop/pymocd/releases/tag/v3.0.0) — 2026-04-19

Compare: [v2.0.2...v3.0.0](https://github.com/csilab-ufop/pymocd/compare/v2.0.2...v3.0.0) · GitHub release: yes · On PyPI

**Added**

- PRISM (removed in v5.0.0), a discrete multi-objective particle swarm
  optimizer for community detection, registered as the `Prism` class. It was
  developed under the name UFOP. Commit subjects describe these components:
    - a non-dominated archive with crowding-distance pruning;
    - adaptive Clerc–Kennedy constriction;
    - a memetic Louvain step on the Pareto front;
    - LPA seeding;
    - a multi-resolution Q_γ objective (γ = 0.5, 1, 2) as the default
      evaluator.
- The release notes also report an AMI comparison of PRISM with Leiden,
  Louvain and HP-MOCD on LFR graphs with μ ≥ 0.5.
- Benchmarks:
    - a real-world network harness;
    - a two-dimensional sweep over μ and the number of nodes, on 10k–100k-node
      LFR graphs, with facet plots;
    - shared plotting and theming utilities.

**Changed**

- Release builds keep debug symbols for profiling.
- The release notes and a commit subject (`8faee62`) say that PRISM replaces a
  legacy SMPSO module, which was removed. The v2.0.2 tag contains no SMPSO
  source, so this removal is not visible between the two tags.

!!! warning "Beta notice from the release notes"
    PRISM was released in beta: "API, defaults, and objective formulation may
    change in a subsequent release". The notes advise pinning the exact version
    for reproducibility.

At the tagged commit, `pyproject.toml` and `Cargo.toml` still declare version
2.0.2. A later commit of the same day (`13ce3f8`) set both to 3.0.0.

## [v2.0.2](https://github.com/csilab-ufop/pymocd/releases/tag/v2.0.2) — 2026-03-09

Compare: [v2.0.1...v2.0.2](https://github.com/csilab-ufop/pymocd/compare/v2.0.1...v2.0.2) · GitHub release: yes, published 2026-04-05 · On PyPI, without the HP-MOCD additions below (see the note)

**Added**

- Custom objectives for HP-MOCD (PR #30). The notes describe this as the
  ability to replace HP-MOCD's default intra/inter objectives with your own
  optimization objectives. At this tag the class has:
    - an `objectives=None` argument to `HpMocd(...)`;
    - [`HpMocd.set_objectives`][pymocd.HpMocd.set_objectives];
    - a per-generation callback,
      [`HpMocd.set_on_generation`][pymocd.HpMocd.set_on_generation];
    - a `num_gens` property.
- A `make bump` target.

**Changed**

- The benchmarks were reorganized into modules for specific benchmarks.
- The documentation domain was updated.

The pymocd 2.0.2 files on PyPI were uploaded on 2026-03-05, four days before
this tag and before the custom-objectives commit (`4f8d51b`, 2026-03-09). In
the 2.0.2 source distribution, `HpMocd` has only `run` and
`generate_pareto_front`, without any of the additions above. The 3.0.0 source
distribution is the first on PyPI whose `HpMocd` has `set_objectives` and
`set_on_generation`.

## [v2.0.1](https://github.com/csilab-ufop/pymocd/tree/v2.0.1) — 2026-03-05

Compare: [v2.0.0...v2.0.1](https://github.com/csilab-ufop/pymocd/compare/v2.0.0...v2.0.1) · GitHub release: none · Not on PyPI

This tag sits on a rewritten history: v2.0.0 is not one of its ancestors. As a
result, the compare view also lists rewritten copies of older commits. The
entries below cover only the commits new in this tag.

**Changed**

- Dependencies: rand, PyO3, ndarray and rand_chacha were bumped to their latest
  versions.
- `Cargo.lock` was removed from the repository, since this is a library crate.
- The Python dependencies were moved.
- Version metadata: a commit titled "Update: 1.3.0 --> 2.0.0" sets the version
  to 2.0.0, so the files at this tag declare 2.0.0. A citation was also added.
- Benchmarks:
    - tqdm progress bars;
    - a shared `BENCHMARK_RUN_ID`;
    - an absolute, timestamped `SAVE_PATH`;
    - a reorganized `lfr_experiment`.

**Fixed**

- The `RngExt` import, for rand 0.10 compatibility.
- Benchmarks:
    - the Leiden igraph mapping, an empty-partition crash and a ddof warning;
    - scripts now run from any directory;
    - the GA parameter-sweep results are saved to `SAVE_PATH`.

## [v2.0.0](https://github.com/csilab-ufop/pymocd/releases/tag/v2.0.0) — 2026-03-04

Compare: [v1.5.2...v2.0.0](https://github.com/csilab-ufop/pymocd/compare/v1.5.2...v2.0.0) · GitHub release: yes · On PyPI

!!! warning "Regression notice from the release notes"
    The maintainer reported "a significant regression in algorithm quality".
    According to the notes, benchmark results showed that detection performance
    had degraded considerably compared to previous versions, with notably lower
    NMI, AMI and modularity scores against baselines such as Louvain and Leiden.
    The root cause was under investigation.

**Added**

- A Makefile, a `requirements.txt` file, a tests package `__init__` and a
  benchmark subfolder.

**Changed**

- The benchmarks were cleaned up and gained a Louvain/Leiden comparison.

**Removed**

- `setup.sh` and `test.py`.

The tagged commit does not descend from v1.5.0, v1.5.1 or v1.5.2. It branches
from a commit of 2025-08-23 (`246ab90`), which predates the v1.5.0 project
restructure. Its `pyproject.toml` and `Cargo.toml` declare version 1.3.0. No
later tag descends from v1.5.x either. The v1.5.0 rework (`df2e9c9`) and
restructure (`302aff0`) commits do not appear in the history of v2.0.0 or of
any later tag, under their own hashes or as rewritten copies with the same
subjects.

## [v1.5.2](https://github.com/csilab-ufop/pymocd/releases/tag/v1.5.2) — 2026-02-22

Compare: [v1.5.1...v1.5.2](https://github.com/csilab-ufop/pymocd/compare/v1.5.1...v1.5.2) · GitHub release: yes · Not on PyPI

**Changed**

- Dependencies:
    - PyO3 0.27.2 to 0.28.1;
    - rand 0.9.2 to 0.10.0;
    - rand_chacha 0.9.0 to 0.10.0.

**Fixed**

- Adapted to the rand 0.10 API, where `random()`, `random_bool()` and
  `random_range()` moved to `RngExt`.
- CI activates the virtual environment before running maturin and the tests.

At the tagged commit, `pyproject.toml` declares version 1.5.1 and `Cargo.toml`
declares 1.5.0.

## [v1.5.1](https://github.com/csilab-ufop/pymocd/releases/tag/v1.5.1) — 2026-01-25

Compare: [v1.5.0...v1.5.1](https://github.com/csilab-ufop/pymocd/compare/v1.5.0...v1.5.1) · GitHub release: yes · On PyPI

**Changed**

- CI compiles the Rust code on every change and runs the tests after
  compilation. The documentation is deployed only for changes inside `docs/`.
- Dependencies: PyO3 0.27.1 to 0.27.2 and ndarray 0.17.1 to 0.17.2.
- Documentation (then built with Hugo): a changelog page was added, and the
  version on the landing page was updated.

**Fixed**

- The citation formatting in the README.
- The license shield now reads the license from the repository instead of PyPI,
  which was failing.

## [v1.5.0](https://github.com/csilab-ufop/pymocd/releases/tag/v1.5.0) — 2025-11-20

Compare: [v1.3.0...v1.5.0](https://github.com/csilab-ufop/pymocd/compare/v1.3.0...v1.5.0) · GitHub release: yes · On PyPI

The commits of this version from `df2e9c9` (2025-09-02) on, and all commits of
v1.5.1 and v1.5.2, are not in the history of v2.0.0 or of any later tag, under
their own hashes or as rewritten copies with the same subjects; see v2.0.0.

**Changed**

- A commit described as a full rework (`df2e9c9`) changed 33 files.
- The project was restructured with separate Rust and Python sources and a
  simplified build system, including a Makefile.
- Dependencies:
    - rayon 1.10.0 to 1.11.0;
    - PyO3 0.25.1 to 0.26.0, 0.27.0 and then 0.27.1;
    - ndarray 0.17.0 to 0.17.1.
- Documentation (then built with Hugo):
    - a documentation update (PR #21);
    - a page on HP-MOCD selection strategies;
    - revised documentation pages, in a commit titled "blog: increase the
      quality a little bit, update outdated information";
    - an updated BibTeX citation.

**Removed**

- The `leoben` function added in v1.3.0, in a commit that removed "unused code,
  outdated functions".
- Duplicated and unused files.
- The Python unittest workflow (`python.yml`), in a commit titled "removed
  tests.yml".

**Fixed**

- The tests.

## [v1.3.0](https://github.com/csilab-ufop/pymocd/releases/tag/v1.3.0) — 2025-07-27

Compare: [v1.2.5...v1.3.0](https://github.com/csilab-ufop/pymocd/compare/v1.2.5...v1.3.0) · GitHub release: yes · On PyPI

The release notes list one pull request: "leobenalg + overall optimizations"
(PR #15).

**Added**

- A `leoben` function, from a commit titled "add leobenalg hpmocd".
- More algorithms for the benchmark comparisons.

**Changed**

- The dependencies were updated.

**Removed**

- The `MOCD` class is no longer registered in the Python module.

**Fixed**

- Benchmark experiments now save their CSV results.

Commits in this cycle also added a `csicea` function (`e4e70a4`) and removed it
again (`a600104`). It is not in the v1.3.0 tag.

## [v1.2.5](https://github.com/csilab-ufop/pymocd/releases/tag/v1.2.5) — 2025-06-08

Compare: [v1.1.3...v1.2.5](https://github.com/csilab-ufop/pymocd/compare/v1.1.3...v1.2.5) · GitHub release: yes · On PyPI

**Fixed**

- Building the internal graph no longer drops nodes without edges (isolated
  nodes) (PR #12).
- The Python unittest.

**Added**

- Handling of the community id of nodes with zero degree.
- A check that the input graph is in a valid format.
- Utility tests, and fixed documentation links (PR #13).

**Changed**

- The citation was updated.

## [v1.1.3](https://github.com/csilab-ufop/pymocd/tree/v1.1.3) — 2025-06-03

Compare: [v1.0.3...v1.1.3](https://github.com/csilab-ufop/pymocd/compare/v1.0.3...v1.1.3) · GitHub release: none · On PyPI

**Added**

- Support for `igraph.Graph` input (PR #8), with unit tests for igraph and
  NetworkX graphs. For the input formats the current version accepts, see
  [Input graphs](guide/input-graphs.md).
- A Python unittest workflow.

**Changed**

- The documentation comments in the Rust sources were expanded, and functions
  were moved out of `lib.rs` into separate files (`macros.rs` and `xfeats.rs`).
  The documentation was updated for igraph compatibility.
- The `MOCD` class is registered in the module again. Its registration had been
  commented out from v1.0.1 to v1.0.3.

**Fixed**

- The documentation's featured-image partial, and the README's igraph
  compatibility note.

## [v1.0.3](https://github.com/csilab-ufop/pymocd/releases/tag/v1.0.3) — 2025-06-01

Compare: [v1.0.2...v1.0.3](https://github.com/csilab-ufop/pymocd/compare/v1.0.2...v1.0.3) · GitHub release: yes · On PyPI

The release notes call this "the first stable release of pymocd". They also
summarize the changes from v0.2.0 to v1.0.2, which are listed under those
versions below.

**Fixed**

- The version number was corrected from 1.0.1 to 1.0.3 in `Cargo.toml`,
  `pyproject.toml` and the documentation.

## [v1.0.2](https://github.com/csilab-ufop/pymocd/tree/v1.0.2) — 2025-06-01

Compare: [v1.0.1...v1.0.2](https://github.com/csilab-ufop/pymocd/compare/v1.0.1...v1.0.2) · GitHub release: none · Not on PyPI

**Fixed**

- `HpMocd` is registered in the Python module again. It was missing in v1.0.1.

**Added**

- `CONTRIBUTING.md`, `CODE_OF_CONDUCT.md` and issue templates.
- Hugo-based Python documentation.
- `Cargo.lock`.

**Changed**

- README corrections and NetworkX format guidance (PR #5, by Lucas-S-Vieira,
  according to the v1.0.3 notes).

## [v1.0.1](https://github.com/csilab-ufop/pymocd/tree/v1.0.1) — 2025-05-31

Compare: [v1.0.0...v1.0.1](https://github.com/csilab-ufop/pymocd/compare/v1.0.0...v1.0.1) · GitHub release: none · On PyPI

**Changed**

- The v1.0.3 release notes report minor optimizations in HP-MOCD's
  evolutionary process for the changes since v0.2.0.
- The README gained an example, and the tests were made easier to reproduce.
  `setup.sh` moved to `unit_tests`, and `pyproject.toml` was fixed.

**Deprecated**

- MOCD. The v1.0.3 notes mark it deprecated, and its source carries a
  `#[deprecated]` attribute from this tag.

At this tag the Python module registers only `set_thread_count` and `fitness`.
The `HpMocd` registration is missing, which was fixed in v1.0.2, and the `MOCD`
registration is commented out.

## [v1.0.0](https://github.com/csilab-ufop/pymocd/tree/v1.0.0) — 2025-05-28

Compare: [v0.2.0...v1.0.0](https://github.com/csilab-ufop/pymocd/compare/v0.2.0...v1.0.0) · GitHub release: none · On PyPI

**Added**

- The `set_thread_count` function, first registered at this tag, and a global
  `debug!` macro.
- Experiments: a parallel experiment and a GA-parameter experiment.
- GNU GPLv3 headers in the files that were missing them.

**Changed**

- Pareto-front plotting now uses the NMI, AMI and modularity metrics.

**Removed**

- The `MMCoMO` placeholder class, in a commit that removed "unused files and
  algorithm".

**Fixed**

- The LFR experiment, the sequential saving of the parameter benchmark,
  `thread.py`, and a library missing from the setup.

## [v0.2.0](https://github.com/csilab-ufop/pymocd/releases/tag/v0.2.0) — 2025-05-20

Compare: [v0.1.0...v0.2.0](https://github.com/csilab-ufop/pymocd/compare/v0.1.0...v0.2.0) · GitHub release: yes · On PyPI

**Fixed**

- Issues when generating the Pareto frontier with `HpMocd`, according to the
  release notes. No file under `src/` changed between v0.1.0 and v0.2.0.

**Changed**

- The PyO3 requirement moved from 0.24.0 to 0.25.0.
- The main branch was renamed `master`, and the README was updated. The commit
  titled "changing main branch to master" also moved the Python scripts from
  `python/` to `unit_tests/` and added LFR and Pareto-front experiment scripts
  and a `setup.sh`.

## [v0.1.0](https://github.com/csilab-ufop/pymocd/releases/tag/v0.1.0) — 2025-04-01

Compare: [v0.0.2...v0.1.0](https://github.com/csilab-ufop/pymocd/compare/v0.0.2...v0.1.0) · GitHub release: yes · On PyPI

The only change since v0.0.2 is the version number: `pyproject.toml` goes from
0.0.1 to 0.1.0 and `Cargo.toml` from 0.1.1 to 0.1.0 (the commit subject reads
"update version 0.0.2 -> 0.1.0"). The release notes cover the changes since
v0.0.1, which are listed under v0.0.2.
They describe `HpMocd`, an NSGA-II-based algorithm "more optimized" in this
version, and `MOCD`, built on PESA-II.

## [v0.0.2](https://github.com/csilab-ufop/pymocd/tree/v0.0.2) — 2025-03-26

Compare: [v0.0.1...v0.0.2](https://github.com/csilab-ufop/pymocd/compare/v0.0.1...v0.0.2) · GitHub release: none · Not on PyPI

**Changed**

- **Breaking:** the project and repository were renamed, with a new logo,
  description and shields. The package and the extension module, both named
  `pyevoea` at v0.0.1, are named `pymocd` from this tag.
- **Breaking:** the algorithm structure was reworked so that every algorithm
  has its own class. The module now registers the classes `HpMocd`, `MOCD` and
  `MMCoMO` together with `fitness`. They replace the functions `nsga_ii`,
  `pesa_ii_maxq` and `pesa_ii_minimax`.
- The `MocdPesaII` class was renamed `MOCD`.
- HP-MOCD optimization, a local search in NSGA-II, and a folder for tests.
- Dependabot was enabled, and the PyO3 requirement moved from 0.23.5 to 0.24.0.

**Fixed**

- The `#[pymodule]` declaration, warnings and the `min_max` signature of
  `MocdPesaII`.

`MMCoMO` is a placeholder: the class has no fields or methods, only a `// todo`
comment. The commit that added it removed a `CoCoMi` class that had been added
after v0.0.1, so no tag contains `CoCoMi`. At the tagged commit,
`pyproject.toml` declares version 0.0.1 and `Cargo.toml` declares 0.1.1.

## [v0.0.1](https://github.com/csilab-ufop/pymocd/releases/tag/v0.0.1) — 2025-03-07

GitHub release: yes, published 2025-03-22 · On PyPI, not as this tagged code (see the note) · First tagged version

**Added**

- The code of the earlier repository, updated for the Rust 2024 edition:
  deprecated rand functions were replaced, outdated dependencies updated and
  the module renamed.
- CI workflows for Rust and for PyPI releases.
- At this tag the package (`pyproject.toml`, `Cargo.toml`) and the extension
  module are named `pyevoea`. The module registers the functions `nsga_ii`,
  `pesa_ii_maxq`, `pesa_ii_minimax` and `fitness`.

The pymocd 0.0.1 source distribution on PyPI, uploaded on 2025-03-21, does not
match this tag. Its module is named `pymocd` and registers `fitness` and the
classes `HpMocd`, `CoCoMi` and `MOCD`. The release notes, published 15 days
after the tag, describe `HpMocd` (NSGA-II) and `MOCD` (PESA-II), as in that
package. In the tagged history, both names appear as Python classes only from
v0.0.2.
