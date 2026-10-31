---
description: HP-MOCD, the library's NSGA-II community detector over Shi's decomposed modularity (intra, inter), with max-modularity selection and pluggable Python objectives.
---

# HP-MOCD

<span class="pm-pill pm-pill--own">pymocd method</span>
<span class="pm-pill pm-pill--neutral">Multi-objective</span>

HP-MOCD (High-Performance Multi-Objective Community Detection) searches for
partitions that trade off the two terms of Newman's modularity $Q$ in the
formulation of [Shi et al. (2010)](https://doi.org/10.1109/CEC.2010.5585987):
$\mathrm{intra}$, the fraction of edges that run between communities, and
$\mathrm{inter}$, the sum of the communities' squared degree shares. Both are
minimised by NSGA-II ([Deb et al., 2002](https://doi.org/10.1109/4235.996017))
over a label-based encoding, with a four-parent majority-vote crossover and a
neighbour-majority mutation. [`hpmocd`][pymocd.hpmocd] returns the member of
the first Pareto front with the highest $Q$;
[`hpmocd_fronts`][pymocd.hpmocd_fronts] returns the whole front. The method is
described in [Santos et al. (2025)](https://doi.org/10.1007/s13278-025-01519-7).

<div class="pm-facts" markdown>

| | |
|---|---|
| **Entry points** | [`hpmocd`][pymocd.hpmocd], [`hpmocd_fronts`][pymocd.hpmocd_fronts], [`HpMocd`][pymocd.HpMocd] (tunable budget, custom objectives) |
| **Objectives** | $\mathrm{intra}$ and $\mathrm{inter}$ (Santos et al., 2025, Eq. 6, after Shi et al., 2010), both minimised; $Q = 1 - \mathrm{intra} - \mathrm{inter}$ |
| **Search** | NSGA-II (Deb et al., 2002) on a label-per-node encoding |
| **Selection rule** | the rank-1 member of highest modularity $Q$ |
| **Front accessor** | [`hpmocd_fronts`][pymocd.hpmocd_fronts] (partitions); [`HpMocd.generate_pareto_front`][pymocd.HpMocd.generate_pareto_front] (partitions with objective values) |
| **Deterministic** | no: there is no seed parameter, and two calls on the same graph can return different partitions (module README, *Determinism*) |
| **Original code** | this library; the module README states there is no external implementation |
| **Paper** | Santos et al. (2025), *Social Network Analysis and Mining* 15(1):110 — [DOI](https://doi.org/10.1007/s13278-025-01519-7) |

</div>

## Objectives

Let $C$ be a partition of the graph's nodes into communities $c$, $m$ the
number of undirected edges, $l_c$ the number of edges with both endpoints in
$c$ (each counted once), and $d_c = \sum_{v \in c} \deg(v)$ the total degree of
$c$, to which every internal edge contributes twice. HP-MOCD minimises the
objective pair that Santos et al. (2025, Eq. 6) take from the formulation of
Shi et al. (2010). The module README writes the same equations as Eqs. 3.5 and
3.6 of Shi et al.'s later journal paper,
[Shi et al. (2012)](https://doi.org/10.1016/j.asoc.2011.10.005):

$$
\mathrm{intra}(C) = 1 - \frac{1}{m} \sum_{c \in C} l_c,
\qquad
\mathrm{inter}(C) = \sum_{c \in C} \left( \frac{d_c}{2m} \right)^{2}.
$$

The pair is an exact decomposition of Newman's modularity:

$$
Q(C) = 1 - \mathrm{intra}(C) - \mathrm{inter}(C).
$$

$\mathrm{intra}$ is lowest when few edges are cut, which favours a few large
communities (a single community gives $\mathrm{intra} = 0$ and
$\mathrm{inter} = 1$); $\mathrm{inter}$ is lowest when the degree is spread
over many small communities. The module README describes the resulting Pareto
front as the trade-off between few cut edges and many small communities.

Two details of the implementation matter when you compare values:

- The objectives count edges. The graph is read without edge weights, so a
  `weight` attribute on a NetworkX graph plays no role.
- A graph with no edges ($m = 0$) evaluates to $(0, 0)$ for every partition.

## Search

**Representation.** An individual holds one community label per node; there is
no decoding step. The initial population draws each node's label uniformly
from $\{0, \dots, n-1\}$, where $n$ is the number of nodes, so the search starts
at the fine end of the front. Isolated nodes carry a label through the search
(mutation cannot move them, since it needs a neighbour) and are returned as
community `-1`.

**Loop.** With $N$ = `pop_size`, each of the `num_gens` generations performs:

1. **Survivor selection.** Fast non-dominated sorting, then crowding distance
   within each front, then a sort by rank (ascending) and crowding distance
   (descending), truncating the pool to $N$ individuals.
2. **Offspring.** $N$ children. Each child draws four distinct parents by
   binary tournament (lower rank wins; equal ranks go to the larger crowding
   distance). With probability `cross_rate` the child is the *ensemble
   crossover* of the four parents: each node takes the community that most
   parents give it, with ties drawn uniformly at random. Otherwise the child is
   a copy of one of the four parents, drawn uniformly.
3. **Mutation.** Every node of the child is drawn independently with
   probability `mut_rate`, and each drawn node moves to the community held by
   most of its neighbours.
4. **Evaluation.** The children are evaluated and appended to the population,
   and the per-generation callback (if any) is called.

After the last generation the search returns the individuals of rank 1. Because
the last generation's children are never ranked (see
[Implementation notes](#implementation-notes)), these are survivors of the
final selection step, so the front holds at most $N$ partitions. Nothing removes
duplicates, so a front can contain the same partition more than once.

## Selecting a partition

[`hpmocd`][pymocd.hpmocd] and [`HpMocd.run`][pymocd.HpMocd.run] return the
rank-1 member that maximises $M - \sum_i f_i$, where $f_1, \dots, f_M$ are the
member's objective values. With the built-in pair this is
$2 - \mathrm{intra} - \mathrm{inter} = 1 + Q$, so the rule picks the
member of highest modularity. No ground truth or target number of communities
is involved.

The whole front is available in two forms:

- [`hpmocd_fronts`][pymocd.hpmocd_fronts] returns every rank-1 partition as a
  `list[dict]`, with the default parameters listed below. Its docstring recommends it for
  comparing HP-MOCD with other detectors on the same footing (best member of
  the front), since comparing `hpmocd`'s single selected partition against
  another detector's best front member handicaps HP-MOCD.
- [`HpMocd.generate_pareto_front`][pymocd.HpMocd.generate_pareto_front] returns
  `[(partition, objectives), ...]`, with the objective values in the configured
  order (`[intra, inter]` by default).

With custom objectives, $M - \sum_i f_i$ is no longer modularity, but `run`
still uses it. In that case take the front from `generate_pareto_front` and
apply your own rule; see [Working with Pareto fronts](../guide/fronts.md) and
[Core concepts](../guide/concepts.md).

!!! note "Each call is a new search"
    `hpmocd`, `hpmocd_fronts`, `HpMocd.run` and `HpMocd.generate_pareto_front`
    each run the search from scratch. Since the search is not seeded, calling
    `run()` after `generate_pareto_front()` on the same object does not, in
    general, return a member of the front you already have. To get both, call
    `generate_pareto_front()` once and select from its result.

## Parameters

[`hpmocd`][pymocd.hpmocd] and [`hpmocd_fronts`][pymocd.hpmocd_fronts] take the
graph and nothing else, and run with the fixed defaults below. To change the
budget, or to supply your own objectives, construct [`HpMocd`][pymocd.HpMocd]:

<!-- pymocd:defaults pymocd.HpMocd -->
| Parameter | Default | Description |
|---|---|---|
| `debug_level` | `0` | `0` is silent; `1` or more prints two lines when the object is constructed (the debug level and a graph summary) and logs every 10th generation and the last. The lines are written by Rust to the process's standard output, bypassing `sys.stdout`. |
| `pop_size` | `100` | NSGA-II population size $N$. Must be at least `4` (see [Implementation notes](#implementation-notes)). |
| `num_gens` | `100` | Number of generations. |
| `cross_rate` | `0.7` | Probability that a child is the ensemble crossover of its four parents rather than a copy of one. Must lie in $[0, 1]$. |
| `mut_rate` | `0.5` | Per-node probability of a neighbour-majority move. Must lie in $[0, 1]$. |
| `objectives` | `None` | List of Python callables `(graph, partition) -> float`, all minimised. A non-empty list replaces the built-in `intra`/`inter` pair. |

!!! note "Defaults and the paper's experimental setting"
    The paper (Santos et al., 2025, Sec. 5.1 and Table 1) reports, for its
    experiments, 100 generations, a population of 100, four parents per
    crossover, a crossover probability of 0.8 and a mutation probability of
    0.2. The shipped `pop_size` and `num_gens` defaults match; the shipped
    `cross_rate` and `mut_rate` defaults are 0.7 and 0.5. To run the paper's
    setting, pass the two rates:
    `pymocd.HpMocd(G, cross_rate=0.8, mut_rate=0.2)`. The arXiv preprint
    ([arXiv:2506.01752](https://arxiv.org/abs/2506.01752)) reports the same
    setting.

The class also provides:

[`set_objectives(objectives)`][pymocd.HpMocd.set_objectives]
:   Replaces the objectives after construction; an empty list reverts to the
    built-in pair.

[`set_on_generation(callback)`][pymocd.HpMocd.set_on_generation]
:   Registers `callback(generation, num_gens, front_size)`, called once per
    generation after the children are merged into the population; `None`
    clears it. `generation` counts from `0`, and `front_size` is the number of
    rank-1 survivors (the new children are not ranked yet).

[`num_gens`][pymocd.HpMocd.num_gens]
:   The number of generations, as a read-only property.

Input graphs are NetworkX graphs with integer node ids, or igraph graphs; see
[Input graphs](../guide/input-graphs.md).

## Complexity

Santos et al. (2025, Sec. 4.5) derive a theoretical time complexity of
$O(G \cdot N_p |V|)$, which the abstract also states; the module README states
no bound. In the paper, $G$ is the number of generations (`num_gens`), $N_p$
the population size (`pop_size`, written $N$ elsewhere on this page) and $|V|$
the number of nodes. The paper derives a per-generation cost of
$O(N_p |V| + N_p \log N_p)$ for sparse graphs ($|E| = O(|V|)$, with a constant
mutation rate), takes the $O(N_p |V|)$ term as dominant when $|V|$ is much
larger than $\log N_p$, and counts non-dominated sorting of the two objectives
as $O(N_p \log N_p)$. The arXiv preprint
([arXiv:2506.01752](https://arxiv.org/abs/2506.01752), v2, Sec. 4.5) gives the
same derivation and refers to its Appendix A for the space complexity.

The shipped non-dominated sort builds the full domination matrix of the
merged pool of up to $2N_p$ individuals, comparing every pair of individuals
(`nsga2/sorting.rs`).

The implementation runs these steps in parallel on the Rayon thread pool:
evaluation of the built-in objectives (across the population, and across
communities when a partition has more than 8), offspring production, the
domination matrix of the non-dominated sort, and mutation sweeps over more than
128 drawn nodes. Python objectives are evaluated sequentially. The pool size is
set with [`max_cores`][pymocd.max_cores]; see
[Performance and reproducibility](../guide/performance.md).

## When to use it

- **To run the method with its shipped defaults.** `hpmocd(G)` runs the fixed
  defaults (see the note under [Parameters](#parameters) on how they compare
  with the paper's experimental setting) and returns the max-$Q$ member of the
  front.
- **To study the modularity trade-off.** The front from `hpmocd_fronts` or
  `HpMocd.generate_pareto_front` spans partitions from few, large communities
  (low $\mathrm{intra}$) to many, small ones (low $\mathrm{inter}$).
- **For objective research.** The source describes `HpMocd` as a generic
  NSGA-II with pluggable Python objectives, for objective research: any list of
  minimised callables replaces the built-in pair while the operators stay the
  same.
- **For selector-free benchmarks.** Compare against other detectors' fronts
  with `hpmocd_fronts`, as its docstring recommends.

Consider another detector when you need the same output on every call:
HP-MOCD has no seed, whereas [RIMPSO](rimpso.md) returns identical output on
repeated calls. Edge weights are not used by HP-MOCD's objectives. For a
side-by-side comparison of the detectors, see
[Choosing a detector](../guide/choosing.md).

## Implementation notes

These points come from the code and the module README and describe the
shipped behaviour.

- **The last generation's children never reach the front.** Children are
  created unranked, and the search returns the merged population without a
  final non-dominated sort, so the rank-1 filter sees only the survivors of the
  last selection step. The module README keeps this as is because it is the
  behaviour the published results were produced with.
- **An early majority settles crossover ties.** The vote for a node stops as
  soon as one community reaches $\lceil 4/2 \rceil = 2$ parents. An even 2–2
  split is therefore decided by which community reached 2 first, in parent
  order, rather than by the uniform tie draw.
- **Mutation has two sweeps that are not equivalent.** With at most 128 drawn
  nodes, the sweep is sequential (Gauss–Seidel): a node sees the moves made
  earlier in the same sweep. Above that, the sweep is parallel (Jacobi): every
  node is scored against the labels from before the sweep.
- **Parameters are not validated.** A `pop_size` of 1, 2 or 3 makes the
  four-parent draw loop forever, and a `cross_rate` or `mut_rate` outside
  $[0, 1]$ raises a `PanicException`. With `pop_size=0` the population is
  empty, and with `num_gens=0` no individual is ever ranked; either way the
  front is empty and `run` raises a `PanicException` ("Empty population").
- **Python objectives share one dict.** A single partition dict is cleared and
  refilled for every individual, so an objective must not keep a reference to
  it between calls. The dict holds the search's raw labels: they are not
  renumbered, and isolated nodes have not yet been mapped to `-1`.
- **Output labels are renumbered.** Returned partitions number communities from
  `0` and assign `-1` to isolated nodes.
- **The search is not deterministic.** The module README traces this to five
  sources, none of which can be switched off: unseeded random number
  generators, Rayon's work distribution, the unsorted order in which fronts are
  peeled during non-dominated sorting, hash-map iteration order, and the
  floating-point summation order of the parallel objective evaluation. It
  reports that, over 1000 runs on the karate club graph, the returned partition
  had 3, 4 or 5 communities.

The module README in the repository,
[`src/core/algorithms/hpmocd/README.md`](https://github.com/csilab-ufop/pymocd/blob/master/src/core/algorithms/hpmocd/README.md),
has the full list, with file and function names.

## Example

The one-line call runs with the default parameters, and the class with the
paper's crossover and mutation rates runs the setting of its experiments. The
class then runs the same search with a smaller budget and returns the front
with its objective values.

```python exec="on" source="above" session="hpmocd" result="text"
import networkx as nx
import pymocd

G = nx.karate_club_graph()
print("hpmocd:", len(set(pymocd.hpmocd(G).values())), "communities")

paper = pymocd.HpMocd(G, cross_rate=0.8, mut_rate=0.2).run()   # Santos et al. (2025), Sec. 5.1
print("paper's setting:", len(set(paper.values())), "communities")

model = pymocd.HpMocd(G, pop_size=40, num_gens=30)   # smaller budget than the default
front = model.generate_pareto_front()                # [(partition, [intra, inter]), ...]
print(len(front), "rank-1 members after", model.num_gens, "generations")
print(" k  intra  inter       Q")
for partition, (intra, inter) in sorted(front, key=lambda member: member[1][0])[::6]:
    k = len(set(partition.values()))
    print(f"{k:2d}  {intra:.3f}  {inter:.3f}  {1 - intra - inter:+.3f}")

best = model.run()   # a new search; returns its highest-Q rank-1 member
communities = [{v for v in best if best[v] == c} for c in set(best.values())]
print("run(): Q =", round(nx.community.modularity(G, communities, weight=None), 3))
```

Down the table, $\mathrm{intra}$ rises while $\mathrm{inter}$ falls, and the
number of communities $k$ tends to grow: the front runs from coarse to fine
partitions. The table prints every sixth member; applied to this front, the
selection rule would pick the member of largest $Q$, which need not be among
the printed rows. `run()` applies the rule to the front of its own, new
search. The last line
recomputes $Q$ with NetworkX, unweighted, because HP-MOCD counts edges rather
than weights.

### Custom objectives

Objectives are Python callables `(graph, partition) -> float`, all minimised.
Below, `cut_fraction` restates the built-in $\mathrm{intra}$ in Python, and
`largest_share` replaces $\mathrm{inter}$ with the share of nodes in the largest
community. The callback records the size of the first front at every
generation.

```python exec="on" source="above" session="hpmocd" result="text"
from collections import Counter

def cut_fraction(graph, partition):
    return sum(partition[u] != partition[v] for u, v in graph.edges()) / graph.number_of_edges()

def largest_share(graph, partition):
    return max(Counter(partition.values()).values()) / len(partition)

custom = pymocd.HpMocd(G, pop_size=20, num_gens=10, objectives=[cut_fraction, largest_share])
sizes = []
custom.set_on_generation(lambda gen, num_gens, front_size: sizes.append(front_size))
custom_front = custom.generate_pareto_front()   # select from this yourself
print("first-front size per generation:", sizes)
for partition, (cut, share) in sorted(custom_front, key=lambda member: member[1][0])[::5]:
    print(f"k={len(set(partition.values())):2d}  cut={cut:.3f}  largest={share:.3f}")
```

Python objectives run sequentially, so keep them cheap or the budget small. The
gallery example [Custom objectives](../examples/objectives/plot_custom_objectives.md)
develops this further.

## Reference

Santos, G. O., Vieira, L. S., Rossetti, G., Ferreira, C. H. G., & Moreira,
G. J. P. (2025). A high-performance evolutionary multiobjective community
detection algorithm. *Social Network Analysis and Mining*, 15(1), 110.
[https://doi.org/10.1007/s13278-025-01519-7](https://doi.org/10.1007/s13278-025-01519-7)

```bibtex
@article{Santos2025,
  author    = {Santos, Guilherme O. and Vieira, Lucas S. and Rossetti, Giulio and Ferreira, Carlos H. G. and Moreira, Gladston J. P.},
  title     = {A high-performance evolutionary multiobjective community detection algorithm},
  journal   = {Social Network Analysis and Mining},
  year      = {2025},
  volume    = {15},
  number    = {1},
  pages     = {110},
  doi       = {10.1007/s13278-025-01519-7},
  url       = {https://doi.org/10.1007/s13278-025-01519-7},
  issn      = {1869-5469},
  date      = {2025-11-18}
}
```

The preprint of this paper is
[arXiv:2506.01752](https://arxiv.org/abs/2506.01752); its arXiv record
carries the DOI above.

The objective pair and the search engine come from:

- Shi, C., Zhong, C., Yan, Z., Cai, Y., & Wu, B. (2010). A multi-objective
  approach for community detection in complex network. *IEEE Congress on
  Evolutionary Computation*, 1–8.
  [https://doi.org/10.1109/CEC.2010.5585987](https://doi.org/10.1109/CEC.2010.5585987)
  (the formulation the HP-MOCD paper cites)
- Shi, C., Yan, Z., Cai, Y., & Wu, B. (2012). Multi-objective community
  detection in complex networks. *Applied Soft Computing*, 12(2), 850–859.
  [https://doi.org/10.1016/j.asoc.2011.10.005](https://doi.org/10.1016/j.asoc.2011.10.005)
  (the equations as numbered in the module README)
- Deb, K., Pratap, A., Agarwal, S., & Meyarivan, T. (2002). A fast and elitist
  multiobjective genetic algorithm: NSGA-II. *IEEE Transactions on Evolutionary
  Computation*, 6(2), 182–197.
  [https://doi.org/10.1109/4235.996017](https://doi.org/10.1109/4235.996017)

See [Citing pymocd](../citing.md) for how to cite the library.

## See also

- API: [`hpmocd`][pymocd.hpmocd], [`hpmocd_fronts`][pymocd.hpmocd_fronts],
  [`HpMocd`][pymocd.HpMocd].
- Guide: [Core concepts](../guide/concepts.md),
  [Working with Pareto fronts](../guide/fronts.md),
  [Choosing a detector](../guide/choosing.md),
  [Evaluating partitions](../guide/evaluation.md),
  [Performance and reproducibility](../guide/performance.md).
- Algorithms: [RIMPSO](rimpso.md), the library's other method;
  [Shi-MOCD](shi-mocd.md), which optimises the same objective pair with
  PESA-II; [Algorithms overview](index.md).
- Gallery: [Custom objectives](../examples/objectives/plot_custom_objectives.md),
  [Choosing from a front](../examples/fronts/plot_choosing_from_front.md),
  [Detector tour](../examples/basics/plot_detector_tour.md).
