---
description: How to choose among pymocd's ten detector entry points, covering the recommended default, when to use HP-MOCD or the HpMocd class, when the re-implemented baselines make sense, and how determinism, front access and graph size differ between detectors.
---

# Choosing a detector

pymocd has ten detector entry points for nine algorithms. Two algorithms,
RIMPSO and HP-MOCD, are this library's own methods. The other seven are
re-implemented baselines, ported from their papers for reproduction and
benchmarking. This page explains which detector to use for which purpose. It
relies only on what the package docstrings, the module READMEs and the earlier
versions of these docs state. For a side-by-side table of objectives, engines
and selection rules, see the [Algorithms overview](../algorithms/index.md).

## Short answer

!!! tip "Decision path"
    1. **Are you reproducing a published baseline, or benchmarking against
       one?** Use that baseline's function (see
       [the baselines](#when-the-baselines-make-sense)). If the paper reports
       the best member of its Pareto front, collect the front with the
       matching `*_fronts` accessor. For Shi-MOCD, pass the paper's parameters
       explicitly.
    2. **Do you want to optimize objective functions of your own?** Use the
       [`HpMocd`][pymocd.HpMocd] class with `objectives=[...]`.
    3. **Do you need HP-MOCD
       ([Santos et al., 2025](https://doi.org/10.1007/s13278-025-01519-7))?**
       Use [`hpmocd`][pymocd.hpmocd], or
       `HpMocd(G, cross_rate=0.8, mut_rate=0.2).run()` for the paper's
       experimental setting.
    4. **Otherwise**, use [`rimpso`][pymocd.rimpso].
        - To see the partitions at every resolution, use
          [`rimpso_fronts`][pymocd.rimpso_fronts].
        - To score candidate partitions from elsewhere with RIMPSO's selection
          rule, use [`rimpso_select`][pymocd.rimpso_select].
        - To fly an independent run, change `seed`.

## Decision table

| If you want | Use | Why |
|---|---|---|
| One partition, with no resolution or community count to set | [`rimpso`][pymocd.rimpso] | recommended default; label-free selection with no parameter |
| The partitions at every resolution | [`rimpso_fronts`][pymocd.rimpso_fronts] | RIMPSO's Pareto front is the graph's resolution profile |
| The same output on every call, and independent runs on request | [`rimpso`][pymocd.rimpso] with `seed` | the only detector with a `seed` argument; deterministic at any thread count for a fixed node and edge order |
| HP-MOCD at the library's fixed defaults | [`hpmocd`][pymocd.hpmocd] | `pop_size` and `num_gens` match the paper (100/100); `cross_rate=0.7` and `mut_rate=0.5` differ from its CP = 0.8, MP = 0.2 (Sec. 5.1, Table 1), which [`HpMocd`][pymocd.HpMocd] accepts; maximum-modularity selection |
| HP-MOCD with another budget, per-member objective values, or a progress callback | [`HpMocd`][pymocd.HpMocd] | tunable constructor, `generate_pareto_front`, `set_on_generation` |
| Objective functions of your own | [`HpMocd`][pymocd.HpMocd] with `objectives=[...]` | Python callables replace the built-in objective pair |
| A baseline paper's best-in-front numbers | [`ccm_fronts`][pymocd.ccm_fronts], [`krm_fronts`][pymocd.krm_fronts], [`moga_net_fronts`][pymocd.moga_net_fronts], [`mmcomo_fronts`][pymocd.mmcomo_fronts] | those papers report the best member of the front |
| A selector-free comparison between detectors | any `*_fronts` accessor | compares the searches, not the selection rules |
| To test a selection rule on candidates from another source | [`rimpso_select`][pymocd.rimpso_select] | runs RIMPSO's rule without its search |
| A single-objective modularity baseline | [`gdpso`][pymocd.gdpso] | maximizes Newman–Girvan modularity directly |
| A random-walk baseline that returns one deterministic partition | [`cdrme`][pymocd.cdrme] | single linkage objective (Eq. 12 of its paper); deterministic for a fixed node and edge order, with no seed; see [CDRME](../algorithms/cdrme.md#when-to-use-it) |
| One search compared under two published selection rules | [`mocd_q`][pymocd.mocd_q] and [`mocd_d`][pymocd.mocd_d] | the same PESA-II search with Shi's two rules |

## The default: RIMPSO

The package docstring recommends `rimpso` or `hpmocd`, and the project README
and the earlier version of these docs call [`rimpso`][pymocd.rimpso] the
recommended default. The earlier docs and the function docstrings describe it
as follows.

- **One run covers every resolution.** RIMPSO splits the Constant Potts Model
  into two objectives, `cut` and `pair`. Maximizing CPM at any resolution
  $\gamma$ is minimizing a weighted sum of the two, so the Pareto front the
  swarm builds is the graph's
  whole resolution profile, and $\gamma$ "stops being a parameter the caller
  has to guess" (function docstring).
- **The partition is chosen for you.** The selection rule needs no ground
  truth, no target number of communities and no parameter. It returns the
  front member that a degree-corrected assortative block model fits best once
  its free densities are penalized. The single-community partition loses on
  the penalty and all-singletons is rejected by the same rule's assortativity
  check, so there is no degeneracy filter and no fallback stage.
- **The front is available.** [`rimpso_fronts`][pymocd.rimpso_fronts] returns
  every member with its $(\operatorname{cut}, \operatorname{pair})$ point and
  the index the selector picks.
- **Results are reproducible.** The function docstring states that the same
  graph and the same parameters, `seed` included, give the same partition on
  any number of threads. Here the same graph means the same nodes and edges,
  listed in the same order; listing them in another order can change the
  partition (see [Determinism and seeds](#determinism-and-seeds)).

The selector does not see ground truth, so it does not always return the front
member that best matches a reference partition. The RIMPSO module README traces
the remaining gap to cases where that member does not compress the graph,
such as LFR benchmarks at high mixing and Zachary's karate club. If you have a
selection rule of your own, apply it to the output of
[`rimpso_fronts`][pymocd.rimpso_fronts]. The [RIMPSO page](../algorithms/rimpso.md)
gives the figures the README reports.

## HP-MOCD and the `HpMocd` class

[`hpmocd`][pymocd.hpmocd] runs HP-MOCD, published by
[Santos et al. (2025)](https://doi.org/10.1007/s13278-025-01519-7), with fixed
defaults, and returns the member of the rank-1 front with the highest
modularity. It takes the graph and nothing else. The defaults (`pop_size=100`,
`num_gens=100`, `cross_rate=0.7`, `mut_rate=0.5`) match the population and
generation counts of the paper's experiments, but the paper used a crossover
probability of 0.8 and a mutation probability of 0.2 (Sec. 5.1 and Table 1);
`HpMocd(G, cross_rate=0.8, mut_rate=0.2).run()` runs that setting. The
[HP-MOCD page](../algorithms/hpmocd.md#parameters) lists the defaults and
compares them with the paper's setting.
[`hpmocd_fronts`][pymocd.hpmocd_fronts] runs the same configuration and returns
the whole rank-1 front instead.

The [`HpMocd`][pymocd.HpMocd] class runs the same NSGA-II with a tunable
budget (`pop_size`, `num_gens`, `cross_rate`, `mut_rate`). It also offers
three things the functions do not:

- **Custom objectives.** `objectives` takes a list of Python callables
  `(graph, partition) -> float`, all minimized, which replace the built-in
  `intra`/`inter` pair. The source describes the class as a generic NSGA-II
  with pluggable Python objectives, for objective research.
- **Objective values.** [`generate_pareto_front`][pymocd.HpMocd.generate_pareto_front]
  returns `[(partition, objectives), ...]`.
- **Progress reporting.** [`set_on_generation`][pymocd.HpMocd.set_on_generation]
  registers a callback `(gen, num_gens, front_size) -> None`.

!!! warning "With custom objectives, `run` does not select by modularity"
    [`HpMocd.run`][pymocd.HpMocd.run] returns the front member with the
    smallest sum of objective values. For the built-in pair that is the member
    with the highest modularity, because
    $Q = 1 - \mathrm{intra} - \mathrm{inter}$;
    for custom objectives the sum is not modularity. The module README advises calling
    `generate_pareto_front` and selecting outside. Python objectives are also
    evaluated sequentially (module README).

HP-MOCD is not deterministic: two calls on the same graph can return different
partitions, and there is no seed parameter. The gallery example
[Custom objectives with HpMocd](../examples/objectives/plot_custom_objectives.md) shows the
class in use.

## When the baselines make sense

The seven re-implemented baselines (CDRME, MMCoMO, NSGA-III-CCM, NSGA-III-KRM,
GDPSO, Shi-MOCD, MOGA-Net) are meant for papers and benchmarks. They were
ported from their papers so that benchmarks in this repository compare against
a faithful implementation rather than a paraphrase. Their pages are listed in
the [Algorithms overview](../algorithms/index.md#re-implemented-baselines).

### Reproducing a paper

The project README states that the shipped defaults follow each paper
wherever the paper states them. Where the paper is silent, the port fills the
gap and its module README says so: for example, MMCoMO's `beta`, the
`divisions` of NSGA-III-CCM and NSGA-III-KRM, and CDRME's `n_walk` and
`alpha_mut` are not given in their papers. Three cases need attention:

- **Shi-MOCD.** [`mocd_q`][pymocd.mocd_q] and [`mocd_d`][pymocd.mocd_d]
  default to this repository's HP-MOCD-parity benchmark budget
  (`cross_rate=0.9`, `mut_rate=0.1`). To reproduce
  [Shi et al. (2012)](https://doi.org/10.1016/j.asoc.2011.10.005), pass
  `cross_rate=0.6`, `mut_rate=0.4` and the per-graph sizes of their Table 1.
  The archive size is not a parameter: the module README states that it is
  derived as `min(pop_size, 100)`.
- **NSGA-III-KRM.** The paper gives a parameter row per dataset, and only the
  row for its first dataset (Zachary's karate club) is the default. The module
  README lists the other rows; see the [NSGA-III-KRM page](../algorithms/krm.md).
- **GDPSO.** The default `num_gens=250` (with `pop_size=100`) is the authors'
  reference code's setting. The module README states that the paper's
  published tables used 100 particles and 100 generations.

### Best-in-front protocols

Several baseline papers report the member of the front that best matches the
ground truth (the best-NMI member), alongside or instead of the member their
selection rule returns. The detector function returns only the selected
member, so reproducing those tables needs the front:

| Paper reports | Accessor |
|---|---|
| [Shaik et al. (2021)](https://doi.org/10.1007/s42979-020-00382-x), Tables 1–2: the best-NMI and the best-modularity member of the front | [`ccm_fronts`][pymocd.ccm_fronts], [`krm_fronts`][pymocd.krm_fronts] |
| [Pizzuti (2012)](https://doi.org/10.1109/TEVC.2011.2161090), Table 1: the best-NMI member of the front | [`moga_net_fronts`][pymocd.moga_net_fronts] |
| [Zhang et al. (2023)](https://doi.org/10.1109/MCI.2023.3277773), Table IV: the best-NMI member of the front | [`mmcomo_fronts`][pymocd.mmcomo_fronts] |

The same holds when you compare detectors. The
[`hpmocd_fronts`][pymocd.hpmocd_fronts] docstring warns that comparing
`hpmocd`'s single selected partition against another detector's best-in-front
score silently handicaps HP-MOCD, and the same reasoning applies to every
detector.
Compare selected partitions with selected partitions, and fronts with fronts.
A best-in-front score needs ground truth, so it measures what the search
found, not what the detector returns on unlabeled data.
[Working with Pareto fronts](fronts.md) and the gallery example
[Choosing a partition from a front](../examples/fronts/plot_choosing_from_front.md)
show the protocol.

!!! warning "Shi-MOCD on graphs without edges"
    Two of the known issues that the Shi-MOCD module README records, and
    deliberately leaves unfixed, concern graphs with no edges. On such a
    graph, `mocd_q` raises `RuntimeError("MOCD-Q produced an empty front")`,
    and `mocd_d` raises `pyo3_runtime.PanicException`, which
    `except Exception` does not catch.

## Single- or multi-objective

A multi-objective detector searches for a Pareto front of trade-offs and then
returns one member of it, so its answer depends on both the search and the
selection rule. Seven of the nine algorithms are multi-objective. Two are not:

- [`gdpso`][pymocd.gdpso] maximizes Newman–Girvan modularity alone.
- [`cdrme`][pymocd.cdrme] maximizes one scalar, the sum of the two linkage
  measures in its paper's Eq. (12). The module README notes that the
  "multi-objective" in CDRME's name refers to those two measures being summed.

Neither has a front. The selection rules also differ. Seven entry points
(`hpmocd`, `cdrme`, `mmcomo`, `ccm`, `krm`, `mocd_q` and `moga_net`) return the
candidate with the highest modularity, and `gdpso` maximizes modularity
directly. Modularity has a resolution limit
([Fortunato and Barthélemy, 2007](https://doi.org/10.1073/pnas.0605965104)):
modularity optimization may fail to identify communities smaller than a scale
that depends on the total size of the network and on how interconnected the
communities are. The GDPSO module README notes that GDPSO inherits this limit with no
mitigation, and the RIMPSO module README measured maximum modularity as a
selection rule on RIMPSO's own front and rejected it for this reason. The two
entry points that select differently are `rimpso` (a block-model fit) and
`mocd_d` (distance from random-graph control fronts).

## Determinism and seeds

| Detector | Repeated calls agree | `seed` argument | Threads |
|---|---|---|---|
| [`rimpso`][pymocd.rimpso] | yes | yes | parallel; result independent of the thread count |
| [`cdrme`][pymocd.cdrme] | yes | no | parallel; result independent of the thread count |
| [`gdpso`][pymocd.gdpso] | yes | no | parallel; result independent of the thread count |
| [`hpmocd`][pymocd.hpmocd], [`HpMocd`][pymocd.HpMocd] | no | no | parallel; scheduling is one source of variation |
| [`mmcomo`][pymocd.mmcomo] | no | no | single-threaded |
| [`ccm`][pymocd.ccm], [`krm`][pymocd.krm] | no | no | single-threaded |
| [`mocd_q`][pymocd.mocd_q], [`mocd_d`][pymocd.mocd_d] | no | no | single-threaded |
| [`moga_net`][pymocd.moga_net] | no | no | single-threaded |

The table follows the *Determinism* section of each module README. The
nondeterministic detectors draw from a thread-local generator seeded by the
operating system, and none of them exposes a seed. Of the three deterministic
detectors, only `rimpso` takes a `seed`. Its default, `0`, reproduces the
trajectory RIMPSO flew before `seed` was a parameter, and any other value
flies an independent one. `cdrme` and `gdpso` return the same partition for
the same graph and parameters, and offer no way to request a different run.

For all three, "the same graph" means the same graph object, or one that
lists the same nodes and edges in the same order. The partition can change if
the same nodes or edges are listed in a different order, for example when the
graph is loaded from another file or built by another library. When results
must be reproduced, keep both orders fixed, for example by adding the nodes in
sorted order and then the edges in sorted order (see
[Input order](performance.md#input-order)).

To characterize a nondeterministic detector, call it several times and report
the spread, as the module READMEs do for their karate-club reference
distributions (for example, 30 runs each for NSGA-III-CCM, NSGA-III-KRM,
MOGA-Net and Shi-MOCD). Two further points from those READMEs:

- On a small graph, a nondeterministic detector can return the same partition
  on every call. The READMEs attribute this to the objective landscape of the
  graph, not to determinism.
- `mocd_d` varies more than `mocd_q`, because its answer also depends on
  `rand_networks` freshly drawn random control graphs.

The block below checks which detectors have a front accessor and a `seed`
argument, whether two calls to `rimpso` agree, and whether `rimpso` returns
the member that `rimpso_fronts` marks as selected. It runs against the
installed package when this page is built.

```python exec="on" source="above" session="choosing" result="text"
import inspect

import networkx as nx
import pymocd

detectors = ["rimpso", "hpmocd", "cdrme", "mmcomo", "ccm",
             "krm", "gdpso", "mocd_q", "mocd_d", "moga_net"]
print(f"{'detector':<10}{'front accessor':<17}seed argument")
for name in detectors:
    params = inspect.signature(getattr(pymocd, name)).parameters
    front = name + "_fronts" if hasattr(pymocd, name + "_fronts") else "-"
    print(f"{name:<10}{front:<17}{'seed' in params}")

G = nx.karate_club_graph()
first = pymocd.rimpso(G)
second = pymocd.rimpso(G)
print()
print("rimpso, two calls agree:", first == second)

partitions, points, selected = pymocd.rimpso_fronts(G)
print("rimpso equals rimpso_fronts at the selected index:",
      partitions[selected] == first)
```

See [Performance and reproducibility](performance.md) for threads and seeds in
more detail.

## Access to the Pareto front

| Detector | Front accessor | What it returns |
|---|---|---|
| RIMPSO | [`rimpso_fronts`][pymocd.rimpso_fronts] | a tuple `(partitions, points, selected)`: the members, their `(cut, pair)` points, and the index the selector picks |
| RIMPSO | [`rimpso_select`][pymocd.rimpso_select] | for partitions you supply, the selected index and each candidate's `(cut, pair)` point |
| HP-MOCD | [`hpmocd_fronts`][pymocd.hpmocd_fronts] | a list of partitions |
| HP-MOCD | [`HpMocd.generate_pareto_front`][pymocd.HpMocd.generate_pareto_front] | a list of `(partition, objectives)` pairs |
| MMCoMO, NSGA-III-CCM, NSGA-III-KRM, MOGA-Net | [`mmcomo_fronts`][pymocd.mmcomo_fronts], [`ccm_fronts`][pymocd.ccm_fronts], [`krm_fronts`][pymocd.krm_fronts], [`moga_net_fronts`][pymocd.moga_net_fronts] | a list of partitions |
| GDPSO, CDRME | none | single objective, so there is no front |
| Shi-MOCD | none | the front exists internally but is not exposed |

Each call runs a new search. For a nondeterministic detector, the partition
that the detector function returns is therefore not necessarily a member of
the front that a separate `*_fronts` call returns. `rimpso` is deterministic,
and its docstring states that `rimpso_fronts` searches identically, so with
the same arguments `rimpso` returns the member that `rimpso_fronts` marks as
selected (the check in the previous section prints this). `rimpso_fronts`
returns the selected index together with the front, so one call gives both. See
[Working with Pareto fronts](fronts.md).

## Graph size, memory and threads

The module READMEs and docstrings state the following limits and costs.

- **MMCoMO** keeps a dense $n \times n$ diffusion kernel and a dense vote
  matrix, both $O(n^2)$ in memory, and computes the kernel with an $O(n^3)$
  eigendecomposition. The README calls this the cap on the detector's usable
  graph size, and the repository's synthetic benchmark registers MMCoMO with
  `max_nodes=2000`.
- **`mocd_d`** runs `rand_networks + 1` full searches: one on the graph and one
  on each random control graph, each with the full budget.
- **CDRME**'s cost is linear in `pop_size`.
- **Single-threaded baselines.** MMCoMO, NSGA-III-CCM, NSGA-III-KRM, Shi-MOCD
  and MOGA-Net run on one thread. For the last four, the module READMEs say
  this is deliberate, so that their cost tracks the published method rather
  than this library's infrastructure. RIMPSO, HP-MOCD, CDRME and GDPSO run on
  Rayon's global thread pool, which [`max_cores`][pymocd.max_cores] sizes.
  Its docstring states that only the first call in a process takes effect.
  Call it before any of these four detectors runs: on the current package, a
  first call to `max_cores` made after one of them has run raises
  `pyo3_runtime.PanicException`, because the pool already exists.
- **RIMPSO**'s module README reports selector measurements that include the
  five SNAP `com-` graphs, of up to 4 million vertices.

The other module READMEs do not state an overall time or memory bound; the
*Complexity* section of each algorithm page records what is known.
<!-- TODO(docs): overall time/memory bounds for RIMPSO, HP-MOCD, CDRME, NSGA-III-CCM, NSGA-III-KRM, GDPSO, Shi-MOCD and MOGA-Net are not stated in their module READMEs; add them here (and on each algorithm page) from the papers or a measured scaling study. -->

## See also

- [Algorithms overview](../algorithms/index.md): objectives, engines,
  selection rules and original code for all nine algorithms.
- [Core concepts](concepts.md): multi-objective detection, Pareto fronts and
  selection.
- [Evaluating partitions](evaluation.md): NMI, AMI, ARI and F1 against ground
  truth.
- Gallery: [Every detector on one graph](../examples/basics/plot_detector_tour.md)
  and [Comparing detectors on a planted partition](../examples/benchmarks/plot_compare_detectors.md).
