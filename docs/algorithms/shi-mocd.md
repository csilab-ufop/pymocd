---
description: Shi-MOCD in pymocd — a port of Shi et al.'s PESA-II community detector over the decomposed-modularity pair (intra, inter), with the MOCD-Q (max modularity) and MOCD-D (max–min distance to random-graph fronts) selection rules.
---

# Shi-MOCD (MOCD-Q / MOCD-D)

<span class="pm-pill pm-pill--baseline">Re-implemented baseline</span>
<span class="pm-pill pm-pill--neutral">Multi-objective</span>

Shi-MOCD is the multi-objective community detection method of
[Shi, Yan, Cai & Wu (2012)](https://doi.org/10.1016/j.asoc.2011.10.005). It
minimizes the two terms into which Newman's modularity $Q$ splits,
$\mathrm{intra}$ and $\mathrm{inter}$, with the PESA-II evolutionary algorithm
over a locus-based adjacency encoding, and keeps the non-dominated partitions
in an external archive. The paper proposes two rules for picking one partition
from that front, and pymocd exposes each as its own function:
[`mocd_q`][pymocd.mocd_q] returns the front member of highest modularity
(MOCD-Q), and [`mocd_d`][pymocd.mocd_d] returns the member that lies furthest
from the fronts of Erdős–Rényi control graphs (MOCD-D). Both return a single
partition; the front itself is not exposed.

<div class="pm-facts" markdown>

| | |
|---|---|
| **Entry points** | [`mocd_q`][pymocd.mocd_q] (MOCD-Q), [`mocd_d`][pymocd.mocd_d] (MOCD-D) |
| **Objectives** | $\mathrm{intra}$ and $\mathrm{inter}$ (Shi et al., 2012, Eqs. 3.5 and 3.6), both minimized; $Q = 1 - \mathrm{intra} - \mathrm{inter}$ (Eq. 3.7) |
| **Search** | PESA-II (Corne et al., 2001) with an adaptive hyper-grid, over the locus-based adjacency encoding (Park & Song, 1998) |
| **Selection rule** | MOCD-Q: maximum $Q$ (Eq. 3.8). MOCD-D: maximum of the minimum distance to Erdős–Rényi control fronts (Eqs. 3.9–3.11) |
| **Front accessor** | None; there is no `mocd_fronts` |
| **Deterministic** | No: there is no seed parameter, and the random generator is seeded from the operating system (module README, *Determinism*) |
| **Original code** | None published; the module README states that the authors' sites and GitHub were searched and no reference source exists |
| **Paper** | Shi et al. (2012), *Applied Soft Computing* 12(2):850–859 — [DOI](https://doi.org/10.1016/j.asoc.2011.10.005) |

</div>

!!! warning "The defaults are not Shi's configuration"
    [`mocd_q`][pymocd.mocd_q] and [`mocd_d`][pymocd.mocd_d] default to
    `pop_size=100`, `num_gens=100`, `cross_rate=0.9` and `mut_rate=0.1`. Their
    docstrings call this the repository's HP-MOCD-parity benchmark budget, and
    the module README attributes it to the HP-MOCD benchmark of
    [Santos et al. (2025)](https://doi.org/10.1007/s13278-025-01519-7): Table 1
    of that paper runs its MOCD baseline with 100 generations, a population of
    100, crossover probability 0.9 and mutation probability 0.1. These are the
    settings of that paper's MOCD baseline, not those of HP-MOCD itself.
    This budget is **not** the configuration published by Shi et al. (2012), which,
    as the module README quotes from the paper, uses crossover probability
    $p_c = 0.6$ and mutation probability $p_m = 0.4$ in every experiment, with
    population sizes and generation counts set per graph in the paper's
    Table 1.

    To run the paper's rates, pass them explicitly, for example
    `pymocd.mocd_q(G, cross_rate=0.6, mut_rate=0.4)`, and, for a graph that
    appears in the paper's Table 1, set `pop_size` and `num_gens` to its
    values there. The module README states that Table 1 assigns the internal
    population size, the archive size and the number of generations per
    graph; the values themselves are not recorded in the repository.
    <!-- TODO(docs): list Shi et al. (2012) Table 1 (per-graph ipsize / epsize / generations) here; the values are not in src/core/algorithms/mocd/README.md and the paper was not available to the docs writer -->
    `pop_size` sets both the internal population and, capped at 100, the
    archive size (see [Parameters](#parameters)), so an archive size that
    differs from `min(pop_size, 100)` cannot be set.

## Objectives

Let $G = (V, E)$ be the graph with $m = |E|$ edges, $C = \{c_1, \dots, c_k\}$
a partition of $V$, $l_c$ the number of edges with both endpoints in $c$ (each
counted once), and $d_c = \sum_{i \in c} \deg(i)$ the total degree of $c$, to
which every internal edge contributes twice. Shi-MOCD minimizes the pair of
Shi et al. (2012), Eqs. 3.5 and 3.6:

$$
\mathrm{intra}(C) = 1 - \frac{1}{m} \sum_{c \in C} l_c,
\qquad
\mathrm{inter}(C) = \sum_{c \in C} \left( \frac{d_c}{2m} \right)^{2}.
$$

The pair is an exact decomposition of the modularity of
[Newman & Girvan (2004)](https://doi.org/10.1103/PhysRevE.69.026113)
(Shi Eq. 3.7):

$$
Q(C) = 1 - \mathrm{intra}(C) - \mathrm{inter}(C).
$$

The module README warns that the names are easy to misread:
$\mathrm{intra}$ is *one minus* the fraction of internal edges, so minimizing
it maximizes internal edges, and $\mathrm{inter}$ is the degree-concentration
(null-model) term, which a single-community partition drives to 1.
Fragmenting a partition lowers $\mathrm{inter}$ and raises $\mathrm{intra}$;
coarsening does the reverse. The module README identifies this trade-off as
the reason the Pareto front spans the range of resolutions.

The objectives count edges. Edge weights are not read, and self-loops are
dropped when the graph is converted. These are the same two objectives that
[HP-MOCD](hpmocd.md) optimizes by default, with NSGA-II.

## Search

**Representation.** Shi-MOCD uses the locus-based adjacency encoding of
Park & Song (1998), which Shi et al. adopt in their §3.1.3. A genome holds one
gene per node; gene $i$ names one neighbor of node $i$, read as "link node $i$
to that neighbor". The communities are the connected components of the graph
formed by these links. A node without neighbors points at itself and becomes a
singleton. Initialization and mutation draw every gene from the node's
neighbors, and crossover copies each gene from one parent at the same
position, so the module README notes that every genome is valid by
construction and no repair step exists. Decoding is a
union-find over the $n$ links followed by a relabeling of the communities in
order of first appearance.

**PESA-II.** The search is PESA-II
([Corne, Jerram, Knowles & Oates, 2001](#reference)), reimplemented inside the
module. It keeps an internal population of size `pop_size` that is
regenerated every generation and an external archive of capacity
`min(pop_size, 100)` that holds the non-dominated partitions found
so far. The archive is both the result and the mating pool.

1. **Seeding.** `pop_size` random genomes are decoded and evaluated; each is
   offered to the archive, which is then truncated to its capacity.
2. **Generations.** For each of the `num_gens` generations (stopping early if
   the archive is empty):
    1. The archive is binned into a hyper-grid with 8 divisions per objective
       (64 cells), each axis normalized to the range the archive currently
       spans. The occupancy of a cell is its *squeeze factor*.
    2. `pop_size` children are bred. With probability `cross_rate` a child is
       the uniform crossover of two parents (each gene taken from either
       parent with probability 1/2); otherwise it is a copy of one parent.
       Every child is then mutated: each gene, with probability `mut_rate`, is
       redrawn uniformly from the node's neighbors (which may return the same
       neighbor).
    3. Each child is evaluated and offered to the archive. It is rejected only
       if an archive member dominates it; archive members it dominates are
       removed.
    4. If the archive exceeds its capacity, it is truncated.
3. **Result.** The final archive is the Pareto front.

**Region-based selection.** Parents are chosen by PESA-II's region-based
tournament: two occupied cells are drawn, the one with the lower squeeze
factor wins (ties by a coin flip), and a member is drawn uniformly from inside
the winning cell. The module README stresses that the tournament is over cells,
not members: a cell with one member beats a cell with twenty, which spreads the
archive along the front.

**Truncation.** While the archive is over capacity, a most crowded cell is
chosen (ties uniformly at random), a uniformly random member of it is removed,
and the grid is rebuilt before the next removal. Since only dominated
candidates are rejected, members with identical objective values can
accumulate; truncation is what bounds them.

## Selecting a partition

Both functions run the PESA-II search above on the input graph and apply one
of Shi's two model-selection rules to the resulting front $F$. Write
$f(C) = (\mathrm{intra}(C), \mathrm{inter}(C))$ for the objective point of a
partition.

**MOCD-Q** ([`mocd_q`][pymocd.mocd_q], Shi Eq. 3.8) returns the front member
of highest modularity. Since $Q = 1 - \mathrm{intra} - \mathrm{inter}$, this
is

$$
C^{\star} = \operatorname*{arg\,max}_{C \in F} Q(C)
          = \operatorname*{arg\,min}_{C \in F} \bigl( \mathrm{intra}(C) + \mathrm{inter}(C) \bigr).
$$

**MOCD-D** ([`mocd_d`][pymocd.mocd_d], Shi Eqs. 3.9–3.11) compares the front
with fronts evolved on random graphs. It generates $r$ = `rand_networks`
Erdős–Rényi $G(n, m)$ control graphs with the same node set and the same
number of edges as the input, evolves a full front $R_1, \dots, R_r$ on each,
and returns the real front member that is furthest from all of them:

$$
D(C) = \min_{1 \le j \le r} \; \min_{C' \in R_j} \bigl\lVert f(C) - f(C') \bigr\rVert_2,
\qquad
C^{\star} = \operatorname*{arg\,max}_{C \in F} D(C).
$$

The module README explains the idea: a random graph's front shows what "no
community structure" looks like in objective space, so the real point that
deviates most from it is taken as the most meaningful one. The equations above
describe the rule as implemented (`selection.rs`); the module README gives the
same rule in words and attributes it to Shi Eqs. 3.9–3.11.

The result is a `dict[int, int]` mapping each node to a community numbered
from 0. Isolated nodes are reported as community `-1`.

!!! note "Each call is a new search"
    Every call runs the search from scratch, and the search is not seeded. A
    call to `mocd_q` and a call to `mocd_d` on the same graph therefore select
    from the fronts of two independent searches, and `mocd_d` also draws new
    control graphs on every call. The module README reports that over 20 runs
    on the Florentine families graph, `mocd_d` returned between 2 and 7
    communities while `mocd_q` returned 3 every time, and attributes the
    spread to the control fronts rather than to the search.

!!! info "No front accessor"
    The Rust class behind both functions can produce the whole front, but the
    module README notes that it is not registered with the Python module, and
    there is no `mocd_fronts`. To study the $(\mathrm{intra}, \mathrm{inter})$
    front itself, [HP-MOCD](hpmocd.md) optimizes the same pair and exposes its
    front through [`hpmocd_fronts`][pymocd.hpmocd_fronts] and
    [`HpMocd.generate_pareto_front`][pymocd.HpMocd.generate_pareto_front]; see
    [Working with Pareto fronts](../guide/fronts.md).

## Parameters

One table covers both functions; `rand_networks` belongs to `mocd_d` only.

<!-- pymocd:defaults pymocd.mocd_q pymocd.mocd_d -->
| Parameter | Default | Description |
|---|---|---|
| `pop_size` | `100` | PESA-II internal population size: the number of random genomes in the seed population and of children bred per generation. Also sets the archive capacity to `min(pop_size, 100)`. Values below 1 are raised to 1. |
| `num_gens` | `100` | Number of generations. `0` returns a selection from the seeded archive. |
| `cross_rate` | `0.9` | Probability that a child is the uniform crossover of two parents rather than a copy of one. Must lie in $[0, 1]$. Shi et al. use `0.6`. |
| `mut_rate` | `0.1` | Per-gene probability of redrawing the gene uniformly from the node's neighbors. Must lie in $[0, 1]$. Shi et al. use `0.4`. |
| `rand_networks` | `3` | `mocd_d` only. Number of Erdős–Rényi control graphs; each one costs a full evolutionary run with the same budget. The module README attributes the value 3 to Shi §3.2, which uses three control fronts for the artificial networks (and one for the timed real-network runs). |

The archive capacity (`epsize` in PESA-II) is not a parameter: it is derived
as `min(pop_size, 100)` with a floor of 1. The module README states that the
cap of 100 encodes the largest archive in Shi's Table 1, and that the archive
of 200 used in the paper's §4.1.2 cannot be reached through the public API.
The number of grid divisions per objective is fixed at 8.

The `cross_rate` and `mut_rate` defaults differ from the paper; see the
warning at the top of this page.

### Degenerate inputs

The module README records the first two behaviors below as known issues that
are left unfixed on purpose. Every row was checked against the current
package.

| Input | Result |
|---|---|
| A graph with no edges (including an empty or single-node graph) | `mocd_q` raises `RuntimeError("MOCD-Q produced an empty front")`. `mocd_d` panics with "Real Pareto front is empty."; in Python this is `pyo3_runtime.PanicException`, which derives from `BaseException`, so `except Exception` does not catch it. |
| `rand_networks=0` | There are no control fronts, every member gets the same score, and `mocd_d` returns the first member of the archive: an arbitrary partition, not a selected one. |
| `cross_rate` or `mut_rate` outside $[0, 1]$ | Both functions panic (`pyo3_runtime.PanicException`, message `p=… is outside range [0.0, 1.0]`). |
| `num_gens=0` | The seeded archive is returned as the front, and the selection rule is applied to it. |

## Complexity

The module README states that both functions run single-threaded, with no
parallelism anywhere in the module, so the runtime of this baseline tracks the
published single-threaded method rather than the repository's optimizations.
It also states that `mocd_d` costs `rand_networks + 1` full searches: the real
one plus one per control graph, each with the same `pop_size` and `num_gens`.

Within one search, each generation breeds and evaluates `pop_size` children.
Evaluating a child decodes its genome (a union-find over $n$ links) and then
visits every node and each of its neighbors once to compute
$\mathrm{intra}$ and $\mathrm{inter}$. Each parent draw scans the archive for
the members of the winning cell, inserting a child into the archive tests it
for dominance against the archive members, and truncation rebuilds the grid
after each removal. Neither the module README nor the docstrings state an overall
bound on time or memory.
<!-- TODO(docs): overall time/memory complexity of Shi-MOCD — take it from Shi et al. (2012), doi:10.1016/j.asoc.2011.10.005, if the paper states one, or have the maintainer state a bound for this port in src/core/algorithms/mocd/README.md -->

## When to use it

- **As a published baseline.** Shi-MOCD is one of the re-implemented
  baselines (ports of other authors' published methods), and the two
  functions implement Shi's two selection rules on the same search. For new
  analyses, start from [Choosing a detector](../guide/choosing.md).
- **To follow the paper, pass its configuration.** The defaults are a
  benchmark budget (the docstrings' HP-MOCD-parity budget), not Shi's
  settings; see the warning at the top of this page.
- **MOCD-Q or MOCD-D.** MOCD-Q needs one search and selects by modularity.
  MOCD-D selects by distance from a random-graph reference rather than by
  modularity, but costs `rand_networks + 1` searches, and the module README
  reports it as the noisier of the two rules (see
  [Selecting a partition](#selecting-a-partition)).
- **Not when you need the front.** Neither function returns it. HP-MOCD
  optimizes the same objective pair and exposes its front.
- **Not when you need reproducible output.** There is no seed parameter, so
  repeated calls can return different partitions. [RIMPSO](rimpso.md), one of
  the library's own methods, returns identical output on repeated calls; see
  [Performance and reproducibility](../guide/performance.md).
- **Not on graphs without edges**, where both functions fail (see
  [Degenerate inputs](#degenerate-inputs)).

## Divergences from the paper

No reference implementation exists, so the module README lists the following
as deliberate readings of, or departures from, the paper text:

- **Uniform crossover.** Shi names the operator "uniform two-point crossover"
  but describes it as an independent per-gene choice between the two parents.
  The port follows the description, which is uniform crossover, not
  two-segment crossover. A gene inherited at its own position is valid in
  either parent, so either reading keeps genomes valid.
- **Union-find decoding** instead of the paper's backtracking. The connected
  components, and therefore the partition, are the same.
- **Erdős–Rényi $G(n, m)$ null model.** MOCD-D's control graphs are uniformly
  random simple graphs with the same node set and edge count, Shi's "random
  networks with the same scale" (§3.2). The module README records that a
  degree-preserving null model was considered and rejected because it is not
  what the paper describes. All $n$ nodes, including isolated ones, are
  inserted into each control graph, so its node count matches the input.
- **Full budget for the control fronts.** Each control search uses the same
  `pop_size` and `num_gens` as the real one. The module README explains that
  an under-evolved random front never reaches the fragmented region, which
  makes the real front's fragmented extreme look spuriously deviant and leads
  MOCD-D to the wrong answer.
- **A separate, single-threaded engine.** PESA-II is reimplemented inside the
  module instead of reusing the repository's shared NSGA-II/NSGA-III code.
- **Derived archive size and benchmark rates.** The archive size is not
  exposed (see [Parameters](#parameters)), and the Python functions ship the
  HP-MOCD benchmark rates instead of Shi's.

## Example

The example runs both selection rules on Zachary's karate club at the
defaults and MOCD-Q again at the paper's rates. For each result it prints the
number of communities $k$, the two objectives computed from the partition with
NetworkX, and $1 - \mathrm{intra} - \mathrm{inter}$ next to NetworkX's own
modularity, which should agree (Eq. 3.7). Edge weights are ignored, matching
pymocd.

```python exec="on" source="above" session="shi-mocd" result="text"
import networkx as nx
import pymocd

G = nx.karate_club_graph()
m = G.number_of_edges()

def describe(name, partition):
    groups = {}
    for node, community in partition.items():
        groups.setdefault(community, set()).add(node)
    intra = 1 - sum(G.subgraph(c).number_of_edges() for c in groups.values()) / m
    inter = sum((sum(d for _, d in G.degree(c)) / (2 * m)) ** 2 for c in groups.values())
    q = nx.community.modularity(G, groups.values(), weight=None)
    print(f"{name:<24} k = {len(groups)}   intra = {intra:.4f}   inter = {inter:.4f}"
          f"   1 - intra - inter = {1 - intra - inter:.4f}   Q = {q:.4f}")

describe("mocd_q (defaults)", pymocd.mocd_q(G))
describe("mocd_d (defaults)", pymocd.mocd_d(G))
# Shi et al. (2012) rates; the paper's Table 1 also sets pop_size and num_gens per graph.
describe("mocd_q (pc=0.6, pm=0.4)", pymocd.mocd_q(G, cross_rate=0.6, mut_rate=0.4))
```

MOCD-D draws new control graphs on every call. The next block repeats it five
times and prints the number of communities of each result; the counts need
not agree.

```python exec="on" source="above" session="shi-mocd" result="text"
print([len(set(pymocd.mocd_d(G).values())) for _ in range(5)])
```

## Reference

C. Shi, Z. Yan, Y. Cai and B. Wu, "Multi-objective community detection in
complex networks", *Applied Soft Computing*, vol. 12, no. 2, pp. 850–859,
2012.
[doi:10.1016/j.asoc.2011.10.005](https://doi.org/10.1016/j.asoc.2011.10.005)

```bibtex
@article{Shi_2012,
  title={Multi-objective community detection in complex networks},
  volume={12},
  ISSN={1568-4946},
  url={http://dx.doi.org/10.1016/j.asoc.2011.10.005},
  DOI={10.1016/j.asoc.2011.10.005},
  number={2},
  journal={Applied Soft Computing},
  publisher={Elsevier BV},
  author={Shi, Chuan and Yan, Zhenyu and Cai, Yanan and Wu, Bin},
  year={2012},
  month=Feb,
  pages={850–859}
}
```

Other works cited on this page:

- D. W. Corne, N. R. Jerram, J. D. Knowles and M. J. Oates, "PESA-II:
  Region-based selection in evolutionary multiobjective optimization", in
  *Proc. Genetic and Evolutionary Computation Conference (GECCO 2001)*,
  pp. 283–290, 2001.
- Y. Park and M. Song, "A genetic algorithm for clustering problems", in
  *Proc. 3rd Annual Conference on Genetic Programming*, pp. 568–575, 1998.
- M. E. J. Newman and M. Girvan, "Finding and evaluating community structure
  in networks", *Physical Review E*, vol. 69, no. 2, 026113, 2004.
  [doi:10.1103/PhysRevE.69.026113](https://doi.org/10.1103/PhysRevE.69.026113)
- G. O. Santos, L. S. Vieira, G. Rossetti, C. H. G. Ferreira and
  G. J. P. Moreira, "A high-performance evolutionary multiobjective community
  detection algorithm", *Social Network Analysis and Mining*, vol. 15, no. 1,
  110, 2025.
  [doi:10.1007/s13278-025-01519-7](https://doi.org/10.1007/s13278-025-01519-7).
  Preprint: [arXiv:2506.01752](https://arxiv.org/abs/2506.01752).

## See also

- API: [`mocd_q`][pymocd.mocd_q] and [`mocd_d`][pymocd.mocd_d]
- [HP-MOCD](hpmocd.md), which optimizes the same objective pair with NSGA-II
  and exposes its front
- [Algorithms overview](index.md)
- [Core concepts](../guide/concepts.md), [Choosing a detector](../guide/choosing.md)
  and [Working with Pareto fronts](../guide/fronts.md)
- [Input graphs](../guide/input-graphs.md) and
  [Performance and reproducibility](../guide/performance.md)
- Gallery: [Detector tour](../examples/basics/plot_detector_tour.md) and
  [Comparing detectors](../examples/benchmarks/plot_compare_detectors.md)
- The module README, with the implementation notes, known issues and file
  layout:
  [`src/core/algorithms/mocd/README.md`](https://github.com/csilab-ufop/pymocd/blob/master/src/core/algorithms/mocd/README.md)
