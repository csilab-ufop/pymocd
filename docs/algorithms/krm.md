---
description: NSGA-III-KRM in pymocd, a re-implementation of Shaik, Ravi and Deb (2021) that runs NSGA-III over kernel k-means, ratio cut and modularity and returns the max-modularity member of its Pareto front.
---

# NSGA-III-KRM

<span class="pm-pill pm-pill--baseline">Re-implemented baseline</span>
<span class="pm-pill pm-pill--neutral">Multi-objective</span>

NSGA-III-KRM, from Shaik, Ravi and Deb (2021)
([doi:10.1007/s42979-020-00382-x](https://doi.org/10.1007/s42979-020-00382-x)),
optimizes three objectives together: kernel k-means (KKM) and ratio cut (RC),
both minimized, and Newman's modularity $Q$, maximized. Its three letters name
those three objectives. KKM is minimized by splitting the graph into
singletons and RC by merging it into one community, so the two pull in
opposite directions and are optimized against each other. The search is
the reference-point-based NSGA-III of Deb and Jain (2014)
([doi:10.1109/TEVC.2013.2281535](https://doi.org/10.1109/TEVC.2013.2281535)),
run over a locus-based genome in which every genome decodes to a valid
partition.
[`krm`][pymocd.krm] returns the member of the final rank-1 Pareto front with
the highest modularity, and [`krm_fronts`][pymocd.krm_fronts] returns the whole
front. The same paper also defines [NSGA-III-CCM](ccm.md), which keeps the
engine and the modularity objective but optimizes Community Score and
Community Fitness instead of KKM and RC.

<div class="pm-facts" markdown>

| | |
|---|---|
| **Entry points** | [`krm`][pymocd.krm], [`krm_fronts`][pymocd.krm_fronts] |
| **Objectives** | Kernel k-means $\mathrm{KKM}$ and ratio cut $\mathrm{RC}$, minimized; modularity $Q$, maximized |
| **Search** | NSGA-III (Deb and Jain, 2014) over a locus-based genome, in a private, sequential engine |
| **Selection rule** | Highest modularity $Q$ on the rank-1 front |
| **Front accessor** | [`krm_fronts`][pymocd.krm_fronts], a `list[dict[int, int]]` |
| **Deterministic** | No. The engine draws from an OS-seeded generator and the API has no seed parameter (module README) |
| **Original code** | None released; the module README reports that the authors' sites and GitHub were searched |
| **Paper** | Shaik, Ravi and Deb (2021), *SN Computer Science* 2(1):1 — [doi:10.1007/s42979-020-00382-x](https://doi.org/10.1007/s42979-020-00382-x) |

</div>

## Objectives

The module README states that the $(\mathrm{KKM}, \mathrm{RC})$ pair the paper
adopts is the objective pair of MODPSO. The [MMCoMO](mmcomo.md) module README
attributes the exact formulation of the same pair to Gong, Cai, Chen and Ma
(2014), "Complex Network Clustering by Multiobjective Discrete Particle Swarm
Optimization Based on Decomposition",
[doi:10.1109/TEVC.2013.2260862](https://doi.org/10.1109/TEVC.2013.2260862);
MMCoMO optimizes the same two formulas with a different engine and without
modularity as a third objective.

Let the partition be $C = \{V_1, \dots, V_k\}$ of the $n$ nodes of the graph,
with $A$ the adjacency matrix, $\bar V_i$ the complement of $V_i$, $m$ the
number of edges, $l_c$ the number of edges inside community $c$ (each counted
once) and $d_c$ the sum of the degrees of the nodes of $c$. Define

$$
L(V_a, V_b) = \sum_{i \in V_a} \sum_{j \in V_b} A_{ij} .
$$

$L(V_i, V_i)$ counts each internal edge from both of its endpoints, so it is
twice the number of edges inside $V_i$, and
$L(V_i, \bar V_i) = \sum_{v \in V_i} \deg(v) - L(V_i, V_i)$ is the cut of
$V_i$. The three objectives are

$$
\mathrm{KKM} = 2(n - k) - \sum_{i=1}^{k} \frac{L(V_i, V_i)}{|V_i|}, \qquad
\mathrm{RC} = \sum_{i=1}^{k} \frac{L(V_i, \bar V_i)}{|V_i|}, \qquad
Q = \sum_{c} \left[ \frac{l_c}{m} - \left( \frac{d_c}{2m} \right)^{2} \right].
$$

Denser communities lower KKM, and fewer edges between communities lower RC.
$Q$ is the modularity of Newman and Girvan (2004),
[doi:10.1103/PhysRevE.69.026113](https://doi.org/10.1103/PhysRevE.69.026113).
Inside the search, $n$ counts every node of the graph; an isolated node is a
singleton community there.

**Opposite degenerate extremes.** KKM is minimized by fragmentation: with every
node in its own community, $k = n$ and every $L(V_i, V_i) = 0$, so
$\mathrm{KKM} = 0$. RC is minimized by coarsening: a single community has no
cut, so $\mathrm{RC} = 0$. Neither objective alone rules out its own
degenerate extreme, which is why they are optimized against each other. A unit
test in the module pins both extremes, and a two-community split lying between
them, on a graph of two triangles joined by one edge; the
[example](#example) recomputes those values.

The engine minimizes, so the vector each individual stores is
$(\mathrm{KKM}, \mathrm{RC}, -Q)$, and every rank, dominance test and
reference-point distance is computed in that space. The sign of $Q$ is undone
only when the public function reads the selected member's modularity.

## Search

**Representation.** The genome is the locus-based encoding of Pizzuti (2012),
[doi:10.1109/TEVC.2011.2161090](https://doi.org/10.1109/TEVC.2011.2161090),
which the module README traces to Park and Song (1998). It has one cell per
node, in the graph's stable node order. Cell $u$ holds either $u$ itself or
the position of one of $u$'s neighbors. Decoding is a union-find that joins
$u$ with the cell's value; each connected component is one community. Every
genome that an operator can produce is valid, so there is no repair step
anywhere in the module. Initialization draws each cell uniformly from
$\{u\} \cup N(u)$, where $N(u)$ is the set of neighbors of $u$. Isolated nodes
have no neighbor to point at, decode to singletons, and are reported as
community `-1`.

**Loop.** With `pop_size` $= N$:

1. Draw $N$ random genomes, decode and evaluate them, and sort them into
   non-dominated fronts.
2. Build the structured reference points of Das and Dennis (1998),
   [doi:10.1137/S1052623496307510](https://doi.org/10.1137/S1052623496307510),
   for $M = 3$ objectives and $p$ = `divisions` subdivisions:
   $H = \binom{M + p - 1}{p}$ points, which is 91 at the default $p = 12$.
3. For each of `num_gens` generations:
    1. Breed $N$ offspring. Two parents are drawn uniformly at random, with no
       fitness-based mating selection. With probability `cross_rate` the child
       is their uniform, per-gene crossover; otherwise it is a clone of one of
       the two, chosen with probability 1/2 each. Every gene of the child is
       then independently resampled from $\{u\} \cup N(u)$ with probability
       `mut_rate`.
    2. Run NSGA-III environmental selection over the $2N$ parents and
       offspring (Deb and Jain, 2014, Algorithms 1 to 4). Whole non-dominated
       fronts are kept while they fit. The remaining places go to members of
       the first front that does not fit: the kept members and that front
       together ($S_t$) are translated by the ideal point and scaled by the
       hyperplane intercepts, each is associated with its nearest reference
       line by perpendicular distance, and niching picks the survivors from
       that front. The least crowded reference point wins, with ties among
       equally crowded points broken uniformly at random; an empty niche takes
       its closest member and a non-empty one takes a random member.
    3. Apply the paper's two customizations, in this order: any individual
       whose canonical labeling repeats one already seen is replaced by a
       fresh random genome, and then so is any individual that puts every node
       in one community. Replacements are decoded and evaluated at once.
    4. Re-rank the population. The customizations overwrite individuals after
       selection has ranked them: a replacement starts unranked, and the
       other members' ranks were computed against the individuals it
       displaced.
4. Rank the final population once more. Its rank-1 members are the result.

The engine is that of Deb and Jain (2014),
[doi:10.1109/TEVC.2013.2281535](https://doi.org/10.1109/TEVC.2013.2281535),
with the fast non-dominated sort of Deb et al. (2002),
[doi:10.1109/4235.996017](https://doi.org/10.1109/4235.996017).

### Randomness

NSGA-III-KRM is not deterministic: two consecutive calls can return different
results. The engine and the niching step draw from the thread-local generator,
which is seeded from the operating system, and no seed is exposed. Every
source of variation is a random draw: the initial population, both parent
picks, the crossover coin flips, the mutation resamples, the replacement
genomes of the two customizations, and the two random tie-breaks inside
niching. The search is sequential, so the thread count is not a further source
of variation.

To show the spread, [`krm`][pymocd.krm] was run 30 times with its defaults on
`nx.karate_club_graph()` (pymocd 6.1.0 from PyPI, Python 3.10, on a shared
48-core Linux server). It returned 4 communities in 25 runs and 5 in the other
5. Scored with [`ami`][pymocd.ami] against the two
[karate ground-truth conventions](../guide/evaluation.md#the-karate-ground-truth-convention),
the AMI was 0.545 ± 0.042 (mean ± sample standard deviation; range 0.431 to
0.578) against NetworkX's `club` attribute, and 0.649 ± 0.045 against the
convention that puts node 8 with the Officer faction. A second batch of 300
runs gave 0.533 ± 0.045 and 0.638 ± 0.046, with a mean of 4.19 communities.
See [Performance and reproducibility](../guide/performance.md#nondeterministic-detectors)
for how to handle nondeterministic detectors in experiments.

## Selecting a partition

[`krm`][pymocd.krm] returns the rank-1 member with the highest modularity $Q$.
The docstring calls this the paper's recommended decision rule when no ground
truth is available.

[`krm_fronts`][pymocd.krm_fronts] returns every rank-1 member as a list of
partitions. The paper reports both the best-NMI and the best-modularity member
of the front (its Tables 1 and 2), and reproducing that protocol needs the
whole candidate set; score each member against the ground truth with
[`nmi`][pymocd.nmi].

Each call runs its own search. Because the search is not deterministic, the
partition returned by `krm` need not be a member of the front returned by a
separate `krm_fronts` call. To obtain a front and its selected member from the
same run, call `krm_fronts` and apply the max-modularity rule yourself, as the
[example](#example) below does. See
[Working with Pareto fronts](../guide/fronts.md) for more on fronts.

## Parameters

<!-- pymocd:defaults pymocd.krm pymocd.krm_fronts -->
| Parameter | Default | Description |
|---|---|---|
| `pop_size` | `100` | Population size $N$. The paper's best combination on its dataset D1, Zachary's karate club. The module README requires it to stay at or above the number of reference points, $H = 91$ at the default `divisions`; the code does not check this. |
| `num_gens` | `100` | Number of generations. From the same combination. |
| `cross_rate` | `0.8` | Probability that a child is the uniform crossover of its two parents rather than a clone of one of them. From the same combination. |
| `mut_rate` | `1/34` | Per-gene probability of resampling a cell from $\{u\} \cup N(u)$. From the same combination: $1/n$ on the 34-node karate graph. |
| `divisions` | `12` | Das and Dennis subdivisions $p$; $H = \binom{p + 2}{2}$ reference points for the three objectives, 91 at $p = 12$. The paper does not state a value; the module README takes 12 from pymoo's convention for three objectives. |

Every parameter after `graph` is optional, and [`krm`][pymocd.krm] and
[`krm_fronts`][pymocd.krm_fronts] take the same ones with the same defaults.

!!! note "The budget defaults are the paper's row for one 34-node graph"
    The default `mut_rate` is the constant $1/34 \approx 0.0294$ (the stub
    prints it as `0.029411764705882353`), which the module README describes as
    $1/n$ on the 34-node karate graph. It is not recomputed from the graph you
    pass. Because mutation acts on each gene independently, a child of an
    $n$-node graph has on average $n \cdot$ `mut_rate` resampled genes, which
    is $1$ on karate.

    The paper tunes the parameters separately for each of its four datasets.
    Its Section 8.1 (Parameter Setting) gives the best NSGA-III-KRM
    combination for each: the one, among the sensitivity grid of its Table 3,
    with the highest mean over 10 runs of the product of the best modularity
    and the best NMI. Table 3 fixes 100 generations for every dataset. Only
    the D1 row is shipped as the default; to reproduce another dataset, pass
    its row explicitly:

    | Paper dataset | Nodes | `pop_size` | `num_gens` | `cross_rate` | `mut_rate` |
    |---|---|---|---|---|---|
    | D1, Zachary's karate club (the default) | 34 | 100 | 100 | 0.8 | $1/34$ |
    | D2, bottlenose dolphins | 62 | 200 | 100 | 0.85 | $1/124$ |
    | D3, American college football | 115 | 500 | 100 | 0.9 | $1/230$ |
    | D4, books about US politics | 105 | 400 | 100 | 0.9 | $2/105$ |

    Table 3 offers each dataset the mutation rates $1/n$, $2/n$ and
    $1/(2n)$; the chosen rate is $1/n$ on D1, $1/(2n)$ on D2 and D3, and
    $2/n$ on D4. The rows, the dataset names and sizes, and the section and
    table numbers are taken from the arXiv preprint of the paper,
    [arXiv:2005.03181v1](https://arxiv.org/abs/2005.03181v1) (May 2020). The
    published article is paywalled and was not read, so its numbering may
    differ.

??? info "Internal constants"
    These are fixed in the module and are not exposed. In
    `src/core/algorithms/krm/nsga3/normalize.rs`:
    `NEAR_ZERO = 1e-10` is the threshold below which a Gauss–Jordan pivot, an
    intercept or a fallback scale counts as zero; `MIN_INTERCEPT = 1e-6` is the
    smallest accepted hyperplane intercept, below which a fallback scaling is
    used; and `ASF_OFF_AXIS_WEIGHT = 1e-6` is the off-axis weight of the
    achievement scalarizing function that picks the extreme points. In
    `nsga3/niching.rs`, `DEGENERATE_REF_NORM2 = 1e-30`: a reference direction
    with a smaller squared norm is treated as the origin, and the distance
    becomes the point's own norm.

### Degenerate inputs

The budget arguments are not validated. Every row was checked against the
current package.

| Input | Result |
|---|---|
| `pop_size=0` | `krm` panics with "empty Pareto front"; `krm_fronts` returns `[]`. |
| `cross_rate` or `mut_rate` outside $[0, 1]$ | Both functions panic with `p=… is outside range [0.0, 1.0]`. |
| An empty graph | `krm` returns `{}`; `krm_fronts` returns `pop_size` empty dicts. |
| A graph with no edges | `krm` maps every node to `-1`; `krm_fronts` returns `pop_size` copies of that partition. |

A panic reaches Python as `pyo3_runtime.PanicException`, which derives from
`BaseException`, so `except Exception` does not catch it.

## Complexity

The module README states no bound, and the arXiv preprint of the paper
([arXiv:2005.03181v1](https://arxiv.org/abs/2005.03181v1)) has no complexity
analysis; the paywalled published version was not read. The bounds below are
derived from the source of this implementation. Write $N$ for `pop_size`, $G$ for `num_gens`, $H$ for the
number of reference points (91 at the default `divisions`), and $n$ and $m$
for the numbers of nodes and edges.

**Time.** This implementation takes $O(N(n + m) + N^2 + NH)$ time per
generation and $O(G\,[N(n + m) + N^2 + NH])$ per run. With `pop_size` at least
$H$, as the module README asks, $NH \le N^2$ and the bound per run is
$O(G\,[N(n + m) + N^2])$. The search runs on one thread. Three steps make up
a generation:

- **Evaluation, $O(N(n + m))$.** `make_individual` (`nsga3/engine.rs`) decodes
  a genome with `Locus::decode` (`locus.rs`) and then calls `community_stats`
  (`objectives.rs`) twice, once for KKM and RC and once for $Q$; each call
  walks every node's neighbor list once. A generation evaluates its $N$
  offspring and every individual the two customizations replace, at most
  $2N - 1$ more: up to $N - 1$ duplicates, then up to $N$ single-community
  members. Breeding (`operators.rs`), the canonical labelings the duplicate
  filter hashes, and the single-community test add $O(Nn)$.
- **Non-dominated sorting, $\Theta(N^2)$.** `fast_non_dominated_sort`
  (`nsga3/individual.rs`) tests every ordered pair of individuals for
  dominance. It runs on the $2N$ parents and offspring inside
  `environmental_selection` (`nsga3/survival.rs`), and again on the $N$
  survivors after the customizations.
- **Niching, $O(NH)$.** `niche_select` (`nsga3/niching.rs`) measures the
  distance from each of the at most $2N$ members of $S_t$ to all $H$ reference
  lines, and each of its at most $N$ picks scans all $H$ niche counts. The
  normalization in `nsga3/normalize.rs` is linear in the size of $S_t$ for
  three objectives.

Outside the loop, `Locus::build` costs $O(n + m)$, the initial population $N$
evaluations and one sort, and `api.rs` ranks the final population once more.
Converting the result costs $O(n)$ for `krm` and up to $O(Nn)$ for
`krm_fronts`, which converts every rank-1 member. The decode is a union-find
with path compression but no union by rank or size, so the linear cost of an
evaluation is not guaranteed in the worst case.

**Memory.** $O(Nn + m + N^2 + H)$. Each individual stores its genome ($n$
`usize`) and its labels ($n$ `i32`), 12 bytes per node on a 64-bit build, and
the $2N$ parents and offspring are alive together during environmental
selection: about $24Nn$ bytes. `Locus` stores the $n + 2m$ candidate positions.
The domination lists of the sort hold at most one entry per pair of
individuals, under $2N^2$ for the $2N$-member pool, and the duplicate filter's
hash set holds up to $N$ canonical labelings ($4Nn$ bytes) after the pool is
cut back to $N$. Nothing accumulates across generations.

**Measured scaling.** `krm` was timed on a shared 48-core Linux server
(pymocd 6.1.0 from PyPI, Python 3.10), three runs per setting, each run in a
fresh process. The table gives the median wall time of the `krm` call, the
ratio to the previous row, and the largest peak resident memory of the three
runs above the process's size just before the call. The LFR rows use
`pymocd.lfr_graph(n, 2.5, 1.5, 0.3, average_degree=10, max_degree=50,
min_community=20, max_community=100, seed=1)`; the `pop_size` rows use
Zachary's karate club, where $n + m$ is small enough for the $N^2$ term to
dominate.

| Varied | Graph | Setting | Median time | Ratio | Added peak memory |
|---|---|---|---|---|---|
| $n$ | LFR, $N = G = 100$ | $n$ = 2,000 ($m$ = 9,795) | 1.49 s | | 5.8 MiB |
| | | $n$ = 4,000 ($m$ = 19,667) | 3.20 s | 2.14 | 11.2 MiB |
| | | $n$ = 8,000 ($m$ = 40,227) | 6.91 s | 2.16 | 21.8 MiB |
| | | $n$ = 16,000 ($m$ = 79,847) | 14.3 s | 2.08 | 43.7 MiB |
| | | $n$ = 32,000 ($m$ = 160,335) | 29.3 s | 2.05 | 86.2 MiB |
| | | $n$ = 64,000 ($m$ = 319,998) | 66.5 s | 2.27 | 176.5 MiB |
| $G$ | LFR $n$ = 4,000, $N = 100$ | $G$ = 50 | 1.80 s | | 11.1 MiB |
| | | $G$ = 100 | 3.24 s | 1.80 | 10.8 MiB |
| | | $G$ = 200 | 6.36 s | 1.96 | 11.2 MiB |
| | | $G$ = 400 | 13.4 s | 2.11 | 11.2 MiB |
| $N$ | karate, $G = 100$ | $N$ = 100 | 0.052 s | | 0.4 MiB |
| | | $N$ = 200 | 0.184 s | 3.51 | 0.5 MiB |
| | | $N$ = 400 | 0.470 s | 2.56 | 1.2 MiB |
| | | $N$ = 800 | 1.95 s | 4.16 | 3.0 MiB |
| | | $N$ = 1,600 | 7.80 s | 3.99 | 10.7 MiB |

The measurements agree with the bound, with one drift. Doubling $n$, which
also doubles $m$ on these graphs, roughly doubled the time (ratios 2.05 to
2.27), and memory grew linearly at about 2.8 KiB per node, of which the
population accounts for $24N$ = 2,400 bytes at $N = 100$. The drift: time per
unit of $G\,N\,(n + m)$ rose from 12.7 ns at $n$ = 2,000 to 15.3 ns at
32,000 and 17.3 ns at 64,000. The bound does not predict this, and its cause
was not investigated; the server's load average also rose during the sweep,
from about 12 to 26. Time was linear in $G$ and memory did not depend on it.
On karate, doubling $N$ from 400 to 800 and from 800 to 1,600 quadrupled the
time, the $N^2$ term of non-dominated sorting, and memory grew faster than
linearly, as the $O(N^2)$ domination lists predict.

## When to use it

- **Reproducing or comparing with Shaik, Ravi and Deb (2021).** The
  `pop_size`, `num_gens`, `cross_rate` and `mut_rate` defaults are the paper's
  best combination on karate, the paper's rows for its other three datasets
  are in the [note above](#parameters), and [`krm_fronts`][pymocd.krm_fronts] exposes
  the front for the paper's best-NMI and best-modularity protocol.
- **Comparing engines on the same objective pair.** KRM and
  [MMCoMO](mmcomo.md) optimize the same KKM and RC formulas; KRM adds
  modularity as a third objective and uses NSGA-III. Both return the
  max-modularity member of their rank-1 front.
- **As a baseline, not as a default.** NSGA-III-KRM is one of the library's
  re-implemented baselines, which exist for papers and benchmarks. For general
  use, start from [Choosing a detector](../guide/choosing.md);
  [RIMPSO](rimpso.md) is the library's recommended default.

## Divergences from the paper

The module README lists these deliberate differences:

- The NSGA-III engine is private to the module and sequential, instead of a
  shared engine, and uses no thread parallelism, so that the baseline's cost
  and behavior track the paper's rather than this library's infrastructure.
- Mating draws two parents uniformly at random, with no fitness-based
  selection. The README identifies this as pymoo's NSGA-III for unconstrained
  problems, whose mating comparison ranks only constraint violation and so
  always falls through to a coin flip; there is no binary tournament.
- The ideal point and the extreme points are recomputed at every
  normalization, whereas pymoo carries them across generations. The README
  calls this a known minor difference from the reference the paper's tooling
  used.
- `divisions = 12` resolves an omission in the paper; see
  [Parameters](#parameters).
- Modularity is computed inside the module rather than by the library's shared
  modularity function, so the inner loop never builds a partition object. The
  formula is the same.
- The duplicate filter compares canonical labelings, so it also removes
  partitions that differ only by a permutation of community labels. The README
  notes that without this the encoding's many genomes per partition would
  defeat the filter almost entirely.
- The public functions return normalized partitions in which isolated nodes
  are community `-1`; the paper leaves them as singletons.

## Example

The first block evaluates KKM and RC, written directly from the formulas
above, on the module's unit-test graph: two triangles joined by one edge. It
shows the two degenerate extremes and the split between them.

```python exec="on" source="above" session="krm" result="text"
import networkx as nx
import pymocd


def communities(partition):
    groups = {}
    for node, label in partition.items():
        groups.setdefault(label, set()).add(node)
    return list(groups.values())


def kkm_rc(G, partition):
    """KKM and RC as defined above (assumes no isolated nodes)."""
    groups = communities(partition)
    n, k = G.number_of_nodes(), len(groups)
    within = cut = 0.0
    for V in groups:
        l_in = sum(1 for u in V for v in G[u] if v in V)  # L(V, V)
        within += l_in / len(V)
        cut += (sum(d for _, d in G.degree(V)) - l_in) / len(V)  # L(V, V-bar)
    return 2 * (n - k) - within, cut


T = nx.Graph([(0, 1), (1, 2), (0, 2), (3, 4), (4, 5), (3, 5), (2, 3)])
cases = {
    "singletons": {v: v for v in T},
    "two triangles": {v: int(v >= 3) for v in T},
    "one community": {v: 0 for v in T},
}
for name, p in cases.items():
    kkm, rc = kkm_rc(T, p)
    print(f"{name:>13}: KKM = {kkm:.3f}  RC = {rc:.3f}")
```

The second block runs [`krm`][pymocd.krm] and
[`krm_fronts`][pymocd.krm_fronts] on Zachary's karate club and prints three
members of the front: the coarsest, the max-modularity one, and the finest.
The two calls are independent runs. The modularity is computed with
`weight=None` because `nx.karate_club_graph()` carries edge weights, while the
modularity inside KRM counts edges.

```python exec="on" source="above" session="krm" result="text"
G = nx.karate_club_graph()

part = pymocd.krm(G)        # max-Q member of one run's rank-1 front
print("krm:", len(set(part.values())), "communities")

front = pymocd.krm_fronts(G)  # the rank-1 front of a second, independent run
print("krm_fronts:", len(front), "partitions on the front")

rows = []
for p in front:
    kkm, rc = kkm_rc(G, p)
    q = nx.community.modularity(G, communities(p), weight=None)
    rows.append((len(set(p.values())), kkm, rc, q))

coarsest = min(rows, key=lambda r: r[0])
best_q = max(rows, key=lambda r: r[3])
finest = max(rows, key=lambda r: r[0])
for name, (k, kkm, rc, q) in [("coarsest", coarsest), ("max Q", best_q), ("finest", finest)]:
    print(f"{name:>8}: k = {k:2d}  KKM = {kkm:7.2f}  RC = {rc:6.2f}  Q = {q:.4f}")
```

Compare the KKM and RC columns across the three rows to see how members at
different numbers of communities $k$ trade the two objectives, and where the
max-modularity member falls between the coarsest and the finest.

## Reference

Tanveerulhuq Shaik, Vadlamani Ravi and Kalyanmoy Deb. "Evolutionary
Multi-Objective Optimization Algorithm for Community Detection in Complex
Social Networks." *SN Computer Science* 2(1), article 1, 2021 (published online
7 November 2020).
[doi:10.1007/s42979-020-00382-x](https://doi.org/10.1007/s42979-020-00382-x).

The volume, issue, article number and dates follow the Crossref record of the
DOI. A preprint of the paper by the same three authors is
[arXiv:2005.03181v1](https://arxiv.org/abs/2005.03181v1) (May 2020); this
page's parameter rows and dataset names come from it, because the published
article is paywalled.

```bibtex
@article{Shaik_2020,
  title     = {Evolutionary Multi-Objective Optimization Algorithm for Community Detection in Complex Social Networks},
  volume    = {2},
  ISSN      = {2661-8907},
  url       = {http://dx.doi.org/10.1007/s42979-020-00382-x},
  DOI       = {10.1007/s42979-020-00382-x},
  number    = {1},
  journal   = {SN Computer Science},
  publisher = {Springer Science and Business Media LLC},
  author    = {Shaik, Tanveerulhuq and Ravi, Vadlamani and Deb, Kalyanmoy},
  year      = {2020},
  month     = Nov,
  pages     = {1}
}
```

The entry is the publisher's record from doi.org, which dates the article by
its online publication (2020); the journal issue is dated 2021. The `pages`
field is added: it holds the article number 1 from Crossref, which the doi.org
record leaves out.

Works the implementation builds on:

- Clara Pizzuti. "A Multiobjective Genetic Algorithm to Find Communities in
  Complex Networks." *IEEE Transactions on Evolutionary Computation* 16(3):
  418–430, 2012.
  [doi:10.1109/TEVC.2011.2161090](https://doi.org/10.1109/TEVC.2011.2161090).
- Kalyanmoy Deb and Himanshu Jain. "An Evolutionary Many-Objective
  Optimization Algorithm Using Reference-Point-Based Nondominated Sorting
  Approach, Part I: Solving Problems With Box Constraints." *IEEE
  Transactions on Evolutionary Computation* 18(4): 577–601, 2014.
  [doi:10.1109/TEVC.2013.2281535](https://doi.org/10.1109/TEVC.2013.2281535).
- Kalyanmoy Deb, Amrit Pratap, Sameer Agarwal and T. Meyarivan. "A fast and
  elitist multiobjective genetic algorithm: NSGA-II." *IEEE Transactions on
  Evolutionary Computation* 6(2): 182–197, 2002.
  [doi:10.1109/4235.996017](https://doi.org/10.1109/4235.996017).
- Indraneel Das and J. E. Dennis. "Normal-Boundary Intersection: A New Method
  for Generating the Pareto Surface in Nonlinear Multicriteria Optimization
  Problems." *SIAM Journal on Optimization* 8(3): 631–657, 1998.
  [doi:10.1137/S1052623496307510](https://doi.org/10.1137/S1052623496307510).
- M. E. J. Newman and M. Girvan. "Finding and evaluating community structure
  in networks." *Physical Review E* 69(2): 026113, 2004.
  [doi:10.1103/PhysRevE.69.026113](https://doi.org/10.1103/PhysRevE.69.026113).
- Maoguo Gong, Qing Cai, Xiaowei Chen and Lijia Ma. "Complex Network
  Clustering by Multiobjective Discrete Particle Swarm Optimization Based on
  Decomposition." *IEEE Transactions on Evolutionary Computation* 18(1):
  82–97, 2014.
  [doi:10.1109/TEVC.2013.2260862](https://doi.org/10.1109/TEVC.2013.2260862).
- Park and Song. "A genetic algorithm for clustering problems." *Proceedings
  of the 3rd Annual Conference on Genetic Programming*, pp. 568–575, 1998
  (as cited by the module README).

## See also

- [NSGA-III-CCM](ccm.md): the same paper's first variant, with Community
  Score, Community Fitness and modularity.
- [MMCoMO](mmcomo.md): the KKM and RC pair under a macro-micro co-evolutionary
  NSGA-II.
- API reference: [`krm`][pymocd.krm], [`krm_fronts`][pymocd.krm_fronts],
  [`nmi`][pymocd.nmi].
- Guide: [Core concepts](../guide/concepts.md),
  [Working with Pareto fronts](../guide/fronts.md),
  [Evaluating partitions](../guide/evaluation.md),
  [Choosing a detector](../guide/choosing.md).
- Gallery: [Every detector on one graph](../examples/basics/plot_detector_tour.md),
  which runs `krm` alongside the other detectors, and
  [Choosing a partition from a front](../examples/fronts/plot_choosing_from_front.md).
- [Algorithms overview](index.md) and the
  [module README](https://github.com/csilab-ufop/pymocd/blob/master/src/core/algorithms/krm/README.md),
  which carries the full derivation and file map.
