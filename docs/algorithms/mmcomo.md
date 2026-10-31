---
description: MMCoMO, the macro-micro co-evolutionary multi-objective community detector of Zhang et al. (2023), re-implemented in pymocd from the paper.
---

# MMCoMO

<span class="pm-pill pm-pill--baseline">Re-implemented baseline</span>
<span class="pm-pill pm-pill--neutral">Multi-objective</span>

MMCoMO (Macro-Micro population-based Co-evolutionary Multi-Objective community
detection) minimizes two objectives, kernel *k*-means (KKM) and ratio cut (RC).
Two NSGA-II populations search side by side: a *micro* population that stores
one community label per node, and a *macro* population that stores a set of
community centres and assigns every other node to its most similar centre
through a dense diffusion-kernel similarity matrix. Every `gap` generations
the two populations exchange elite solutions, a modularity local search
refines the micro elites, and the similarity matrix is pulled toward the
co-membership pattern of those elites. At the end, the non-dominated members
of the union of both populations form one Pareto front;
[`mmcomo`][pymocd.mmcomo] returns its member of highest Newman modularity, and
[`mmcomo_fronts`][pymocd.mmcomo_fronts] returns the whole front. pymocd's
port was written from the paper alone; the module README reports that no
original implementation has been published.

<div class="pm-facts" markdown>

| | |
|---|---|
| **Entry points** | [`mmcomo`][pymocd.mmcomo], [`mmcomo_fronts`][pymocd.mmcomo_fronts] |
| **Objectives** | kernel *k*-means (KKM) and ratio cut (RC), both minimized (paper Eq. 1) |
| **Search** | two co-evolving NSGA-II populations, label-based (micro) and centre-based (macro), coupled through a diffusion-kernel similarity matrix |
| **Selection rule** | maximum Newman modularity $Q$ over the merged rank-1 front (the paper's Table III rule) |
| **Front accessor** | [`mmcomo_fronts`][pymocd.mmcomo_fronts] |
| **Deterministic** | no: two calls on the same graph can return different partitions; there is no seed parameter, and every random draw comes from an OS-seeded thread-local generator (module README) |
| **Original code** | none published; the module README reports that neither the authors' site nor GitHub carries code for it |
| **Paper** | Zhang et al. (2023), *IEEE Computational Intelligence Magazine* 18(3):69–86 — [DOI](https://doi.org/10.1109/MCI.2023.3277773) |

</div>

## Objectives

Let $G = (V, E)$ be the graph with $n = |V|$ nodes, $m = |E|$ edges, adjacency
matrix $A$ and degree $\deg(v)$. A partition $C = \{V_1, \dots, V_k\}$ has $k$
communities, where $k$ counts the distinct labels present. For a community
$V_i$ define

$$
L(V_i, V_i) = \sum_{u, v \in V_i} A_{uv},
\qquad
L(V_i, \bar V_i) = \sum_{v \in V_i} \deg(v) - L(V_i, V_i),
$$

so $L(V_i, V_i)$ is twice the number of edges inside $V_i$ and
$L(V_i, \bar V_i)$ is the number of edges leaving it. The two search
objectives (the paper's Eq. 1) are

$$
\mathrm{KKM}(C) = 2(n - k) - \sum_{i=1}^{k} \frac{L(V_i, V_i)}{|V_i|},
\qquad
\mathrm{RC}(C) = \sum_{i=1}^{k} \frac{L(V_i, \bar V_i)}{|V_i|},
$$

both minimized. KKM falls as communities become internally dense per member
and as $k$ rises; RC rises as $k$ rises. That opposition produces the front.
The module README states that the paper attributes this pair to
Shi et al. (2012) and to Gong et al. (2014), and that the exact formulation
implemented is the one of Gong et al. (2014). Both are in the paper's
reference list as indexed by Crossref, as references [10]
([DOI](https://doi.org/10.1016/j.asoc.2011.10.005)) and [11]
([DOI](https://doi.org/10.1109/TEVC.2013.2260862)).

Both objectives are computed over all $n$ nodes, isolated nodes included. An
isolated node adds nothing to either numerator, but it can change $k$ or a
community size: a micro individual starts with it as its own community, which
adds one to $k$; in a macro individual, unless it is itself a centre, decoding
attaches it to an existing centre and it adds one to that community's
$|V_i|$. Isolated nodes therefore perturb the objective values, although the
reported partition always maps them to `-1`.

Newman modularity,

$$
Q(C) = \sum_{i=1}^{k} \left[ \frac{e_i}{m} - \left( \frac{d_i}{2m} \right)^2 \right],
$$

with $e_i$ the number of edges inside $V_i$ and $d_i$ the sum of its degrees,
is maximized but is **not** a search objective. It is the ascent criterion of
the local search and the rule that selects one partition from the final front.

## Search

**Similarity matrix.** Before the search starts, MMCoMO computes the dense
$n \times n$ diffusion kernel of Kondor and Lafferty (2002),

$$
SM = \exp\big(\beta (A - D)\big),
$$

where $D$ is the diagonal degree matrix and $\beta$ is the parameter `beta`.
pymocd computes it with a self-contained cyclic Jacobi eigendecomposition.
$SM$ changes during the run (see the exchange step below).

**Two representations.**

Micro population
:   One community label per node; every label is itself a node index. A new
    individual gives each node the index of a uniformly drawn neighbour (an
    isolated node keeps its own index).

Macro population
:   One bit per node marking it as a community centre. *Decoding* (the paper's
    Eqs. 3–5) labels each centre with itself and assigns every other node $i$
    to the centre $c$ that maximizes $SM_{ic}$. *Encoding* (Eq. 8) inverts
    this: in each community, the member with the greatest summed similarity to
    the other members becomes the centre. A new individual draws a centre count
    $c$ uniformly from $[1, \lceil \sqrt{n} \rceil]$; the first half of the
    population picks its $c$ centres at random among the $\min(3c, n)$
    highest-degree nodes, the second half picks them uniformly.

Each population has `pop_size` members, so `2 * pop_size` individuals are
carried in total.

**Generational loop.** For $t = 1, \dots, T$ with $T$ = `num_gens`:

1. *Micro offspring.* Each child is a copy of a binary-tournament winner
   (lower rank wins, larger crowding distance breaks ties). With probability
   `cross_rate` a second tournament winner donates one community: a random
   node is drawn, and every node that the donor puts in that node's community
   takes the donor's label. Then each node, with probability $1/n$, adopts the
   current label of a random neighbour.
2. *Macro offspring.* Uniform crossover of two tournament winners, then each
   bit flips with probability `mut_rate`. An all-zero child receives one
   uniformly drawn centre. Children are decoded and evaluated.
3. *Survival.* If $t$ is not a multiple of `gap`, each population merges
   parents and children and truncates back to `pop_size` by NSGA-II
   environmental selection (Deb et al., 2002,
   [DOI](https://doi.org/10.1109/4235.996017)), independently.
4. *Exchange*, when $t$ is a multiple of `gap`:
    - **Guidance** (the paper's Algorithm 2). Every rank-1 macro member is
      decoded again with the current $SM$, evaluated, and added to the micro
      pool before micro truncation.
    - **Local search** (Algorithm 1, line 11). The surviving rank-1 micro
      members are improved in place by modularity ascent, then re-evaluated.
    - **Influence** (Algorithm 3). The rank-1 micro members are the elites. A
      vote matrix $SM^{v}$ holds, for each node pair, the fraction of elites
      that place the pair in the same community, and the similarity matrix is
      relaxed toward it (the paper's Eq. 7):

        $$
        SM \leftarrow (1 - \rho)\, SM + \rho\, SM^{v},
        \qquad
        \rho = 0.5\, \frac{t}{T}.
        $$

        Each elite is then encoded, decoded with the updated $SM$, and added to
        the macro pool before macro truncation.

**Mergence.** After the last generation, the union of both populations is
ranked by non-dominated sorting and the rank-1 members form the returned
front. If the front were empty, the all-singletons partition would be returned
instead; the module README notes that this happens only when `pop_size = 0`.

The paper defers the local-search move set to its own reference [38]. pymocd
uses the first phase of the Louvain method (Blondel et al., 2008,
[DOI](https://doi.org/10.1088/1742-5468/2008/10/P10008)): repeated sweeps that
move single nodes to the neighbouring community of largest modularity gain,
capped at 64 sweeps, accepting a move only when its gain exceeds staying by
more than $10^{-12}$.

!!! note "Input handling"
    `mmcomo` and `mmcomo_fronts` read node ids and edges only; edge weights
    are not used. The graph is treated as undirected, and self-loops and
    duplicate edges are dropped when the internal graph is built
    (`graph.rs`). See [Input graphs](../guide/input-graphs.md) for accepted
    input types.

## Selecting a partition

[`mmcomo`][pymocd.mmcomo] returns the front member of highest Newman
modularity $Q$, which the module README identifies as the rule behind the
paper's Table III. [`mmcomo_fronts`][pymocd.mmcomo_fronts] returns the whole
merged front as a `list[dict[int, int]]` instead; the module README notes
that the paper's Table IV uses a best-NMI rule over the front, which needs
this list and a ground truth. Any other rule can be applied to the same list;
see [Working with Pareto fronts](../guide/fronts.md).

Two properties of the returned list matter when you post-process it:

- It is not deduplicated. The front is the rank-1 set of the union of both
  final populations, and the same partition can appear more than once.
  Canonicalize labels before counting distinct candidates (the
  [example](#example) shows one way).
- It can include extreme partitions. On a connected graph, for example, the
  one-community partition has zero ratio cut, so no other partition dominates
  it, and it is on the front whenever a final population contains it.

Since the search is not seeded, a call to `mmcomo` and a call to
`mmcomo_fronts` are two independent runs. To obtain the selected partition
together with the front it came from, call `mmcomo_fronts` once and take the
max-$Q$ member yourself.

## Parameters

<!-- pymocd:defaults pymocd.mmcomo pymocd.mmcomo_fronts -->
| Parameter | Default | Description |
|---|---|---|
| `pop_size` | `100` | Size of **each** population, so `2 * pop_size` individuals are carried. |
| `num_gens` | `50` | Number of generations. All of them always run; there is no convergence stop. |
| `cross_rate` | `0.1` | Probability of the micro one-way community graft. Macro crossover is uniform and always applied. |
| `mut_rate` | `0.1` | Per-bit flip probability of the macro centre genome. Micro mutation is fixed at $1/n$ and is not exposed. |
| `gap` | `10` | Co-evolution interval: the populations exchange, and the local search runs, every `gap` generations. `0` is raised to `1`, so it means every generation. |
| `beta` | `0.05` | Diffusion-kernel exponent $\beta$ in $SM = \exp(\beta(A - D))$. The paper does not specify it; the value is a choice of this re-implementation. `beta = 0` makes $SM$ the identity. |

The source file `defaults.rs` attributes these defaults to the paper's
Table II and Section IV-A, except `beta`, which the paper does not specify.
The parameter descriptions come from the module README.

Some constants are fixed in the code and not exposed: the local search stops
after 64 sweeps and requires a modularity gain above $10^{-12}$; the Jacobi
eigensolver runs at most 100 cyclic sweeps; and macro initialization uses the
centre ceiling $\lceil \sqrt{n} \rceil$ and the $3c$ sampling pool described
above.

### Degenerate inputs

The module README documents a defined result for each trivial case; none of
them raises.

| Input | Result |
|---|---|
| Empty graph | `{}` from `mmcomo` and `[]` from `mmcomo_fronts`; the search is skipped. |
| Single isolated node | `{0: -1}` from `mmcomo`; `mmcomo_fronts` returns `2 * pop_size` identical copies of it, not a one-member front. On a larger graph with no edges, every node is `-1` and the front again repeats that one partition, up to `2 * pop_size` times. |
| `pop_size = 0` | Both populations start empty, so the front is empty and the all-singletons fallback fires: every node is its own community. |
| `num_gens = 0` | The generational loop never runs; the initial populations are merged as they stand. |
| `gap = 0` | Raised to `1`: the populations exchange every generation. |

## Complexity

The module README states that the similarity matrix $SM$ and the vote matrix
built at each exchange both take $O(n^2)$ memory, and that the
eigendecomposition used to compute $SM$ is $O(n^3)$. It identifies these
dense costs as the limit on the graph sizes MMCoMO can handle; the
repository's synthetic-network benchmark script
(`benchmarks/_exp_synt_net/algorithms.py`) registers MMCoMO with
`max_nodes=2000`, and the comment beside that registration gives the dense
kernel and its eigendecomposition as the reason. The module README does not
state a time bound per generation or for a whole run.
<!-- TODO(docs): per-generation / whole-run time complexity of MMCoMO is not stated in src/core/algorithms/mmcomo/README.md; check whether Zhang et al. (2023) give one (the paper was not available to the docs writer) -->

The implementation runs single-threaded: the module README states that the
module contains no parallelism.

## When to use it

- **Reproducing or benchmarking against MMCoMO.** pymocd ships the
  re-implemented baselines for papers and benchmarks. `defaults.rs`
  attributes the defaults to the paper (except `beta`), and the module README
  lists the points where the port fills gaps the paper leaves open (see
  [Divergences from the paper](#divergences-from-the-paper)).
- **Small graphs.** Because $SM$ is dense, memory grows with $n^2$ and the
  kernel computation with $n^3$; the repository's own benchmark script caps
  MMCoMO at 2000 nodes (see [Complexity](#complexity)).
- **When you want the (KKM, RC) front itself.**
  [`mmcomo_fronts`][pymocd.mmcomo_fronts] exposes every non-dominated
  candidate for your own selection rule. NSGA-III-KRM
  ([`krm`](krm.md)) also uses kernel *k*-means and ratio cut, with modularity
  as a third objective.
- **Not when you need reproducible output.** There is no seed parameter, so
  repeated runs are not guaranteed to agree. [RIMPSO](rimpso.md), one of the
  library's own methods, takes a `seed` and returns the same partition on
  repeated calls; see [Choosing a detector](../guide/choosing.md) and
  [Performance and reproducibility](../guide/performance.md).

## Divergences from the paper

The module README lists the following points, each traceable to the code; the
paper is silent or defers to another reference on all of them.

- **`beta` is not from the paper.** The kernel exponent is never specified
  there, so it is a public parameter with default `0.05`.
- **Macro initialization is not specified by the paper,** which defers it to
  its reference [46]. The centre-count range $[1, \lceil \sqrt{n} \rceil]$,
  the high-degree half of the population sampling from the $\min(3c, n)$
  highest-degree nodes, and the uniform second half are this implementation's
  choices.
- **The local-search move set is not specified by the paper,** which fixes
  only its target (rank-1 micro members) and objective ($Q$) and defers to its
  reference [38]. pymocd uses the Louvain first phase with the cap and
  threshold given above.
- **Eq. 4's denominator is skipped.** Decoding takes $\arg\max_c SM_{ic}$
  directly; the normalization in Eq. 4 is constant in $c$, so the argmax is
  unchanged.
- **Two degenerate genomes are repaired.** An all-zero genome decodes with
  the single highest-degree node as its centre; an all-zero *offspring*
  genome instead receives one uniformly drawn centre, so that the repair does
  not collapse the macro population onto one node.
- **The eigendecomposition is self-contained.** $SM$ is computed by cyclic
  Jacobi rotations written in the module rather than through a LAPACK
  binding; the module README gives avoiding that dependency as the reason.
  $A - D$ is explicitly re-symmetrized before the solver runs.
- **Isolated nodes follow a repository convention.** They take part in the
  search and in Eq. 1's sums, but are reported as community `-1`, as in every
  pymocd detector. Their in-search treatment differs between the two
  populations (see [Objectives](#objectives)).
- **`cross_rate` and `mut_rate` are not symmetric.** `cross_rate` applies
  only to the micro graft and `mut_rate` only to macro bit flips; micro
  mutation is fixed at $1/n$, and macro crossover has no probability.
- **Guidance re-decodes.** Algorithm 2, line 5 decodes each rank-1 macro
  elite again with the *current* $SM$ instead of reusing the labels stored on
  the individual. The two can differ, because $SM$ has moved since the
  individual was built.

## Example

Run MMCoMO on Zachary's karate club and inspect the selected partition. The
modularity is computed with `weight=None` because `nx.karate_club_graph()`
carries edge weights, while MMCoMO counts edges (see the *Input handling* note
under [Search](#search)), so the printed $Q$ is the one MMCoMO maximizes.

```python exec="on" source="above" session="mmcomo" result="text"
import networkx as nx
import pymocd

G = nx.karate_club_graph()
part = pymocd.mmcomo(G)  # max-Q member of the merged front


def communities(p):
    groups = {}
    for node, c in p.items():
        groups.setdefault(c, set()).add(node)
    return list(groups.values())


print("communities:", len(set(part.values())))
print(f"modularity Q: {nx.community.modularity(G, communities(part), weight=None):.3f}")
```

The front holds the candidates the selection rule chose from. Count its
distinct members, the range of community counts it covers, and recover the
max-$Q$ member from the same run:

```python exec="on" source="above" session="mmcomo" result="text"
def canonical(p):
    relabel = {}
    return tuple(relabel.setdefault(p[v], len(relabel)) for v in sorted(p))


front = pymocd.mmcomo_fronts(G)
sizes = [len(set(p.values())) for p in front]
best = max(front, key=lambda p: nx.community.modularity(G, communities(p), weight=None))

print("front members:", len(front), "| distinct:", len({canonical(p) for p in front}))
print("communities per member:", min(sizes), "to", max(sizes))
print("max-Q member:", len(set(best.values())), "communities,",
      f"Q = {nx.community.modularity(G, communities(best), weight=None):.3f}")
```

The second block is a separate run, so its max-$Q$ member need not equal the
partition from the first block.

## Reference

Zhang, L., Yang, H., Yang, S., & Zhang, X. (2023). A macro-micro
population-based co-evolutionary multi-objective algorithm for community
detection in complex networks [Research Frontier]. *IEEE Computational
Intelligence Magazine*, 18(3), 69–86.
<https://doi.org/10.1109/MCI.2023.3277773>

```bibtex
@article{Zhang_2023,
  title     = {A Macro-Micro Population-Based Co-Evolutionary Multi-Objective Algorithm for Community Detection in Complex Networks [Research Frontier]},
  volume    = {18},
  ISSN      = {1556-6048},
  url       = {http://dx.doi.org/10.1109/MCI.2023.3277773},
  DOI       = {10.1109/mci.2023.3277773},
  number    = {3},
  journal   = {IEEE Computational Intelligence Magazine},
  publisher = {Institute of Electrical and Electronics Engineers (IEEE)},
  author    = {Zhang, Lei and Yang, Haipeng and Yang, Shangshang and Zhang, Xingyi},
  year      = {2023},
  month     = Aug,
  pages     = {69–86}
}
```

Works the implementation draws on, as cited in the module README:

- Shi, C., Yan, Z., Cai, Y., & Wu, B. (2012). Multi-objective community
  detection in complex networks. *Applied Soft Computing*, 12(2), 850–859.
  <https://doi.org/10.1016/j.asoc.2011.10.005> (the paper's reference [10];
  one of the two sources it attributes the KKM/RC pair to, per the module
  README).
- Gong, M., Cai, Q., Chen, X., & Ma, L. (2014). Complex network clustering by
  multiobjective discrete particle swarm optimization based on decomposition.
  *IEEE Transactions on Evolutionary Computation*, 18(1), 82–97.
  <https://doi.org/10.1109/TEVC.2013.2260862> (the paper's reference [11];
  the exact KKM/RC formulation).
- Deb, K., Pratap, A., Agarwal, S., & Meyarivan, T. (2002). A fast and elitist
  multiobjective genetic algorithm: NSGA-II. *IEEE Transactions on
  Evolutionary Computation*, 6(2), 182–197.
  <https://doi.org/10.1109/4235.996017> (the selection engine).
- Kondor and Lafferty (2002). Diffusion kernels on graphs and other discrete
  input spaces. *Proceedings of ICML 2002*, 315–322 (the similarity kernel; the
  paper's reference [68]).
- Blondel, V. D., Guillaume, J.-L., Lambiotte, R., & Lefebvre, E. (2008). Fast
  unfolding of communities in large networks. *Journal of Statistical
  Mechanics: Theory and Experiment*, 2008(10), P10008.
  <https://doi.org/10.1088/1742-5468/2008/10/P10008> (the local-search move
  set).

## See also

- API: [`mmcomo`][pymocd.mmcomo], [`mmcomo_fronts`][pymocd.mmcomo_fronts]
- [Core concepts](../guide/concepts.md): multi-objective detection, Pareto
  fronts and selection
- [Working with Pareto fronts](../guide/fronts.md)
- [Choosing a detector](../guide/choosing.md)
- [Algorithms overview](index.md) and the related baseline
  [NSGA-III-KRM](krm.md)
- Gallery: [Every detector on one graph](../examples/basics/plot_detector_tour.md),
  [Choosing a partition from a front](../examples/fronts/plot_choosing_from_front.md),
  [Comparing detectors on a planted partition](../examples/benchmarks/plot_compare_detectors.md)
- Module README with the full derivation:
  [`src/core/algorithms/mmcomo/README.md`](https://github.com/csilab-ufop/pymocd/blob/master/src/core/algorithms/mmcomo/README.md)
