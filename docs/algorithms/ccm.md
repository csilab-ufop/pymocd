---
description: NSGA-III-CCM in pymocd, a re-implementation of Shaik, Ravi and Deb (2021) that runs NSGA-III over Community Score, Community Fitness and modularity and returns the max-modularity member of its Pareto front.
---

# NSGA-III-CCM

<span class="pm-pill pm-pill--baseline">Re-implemented baseline</span>
<span class="pm-pill pm-pill--neutral">Multi-objective</span>

NSGA-III-CCM, from Shaik, Ravi and Deb (2021)
([doi:10.1007/s42979-020-00382-x](https://doi.org/10.1007/s42979-020-00382-x)),
maximizes three objectives together: Community Score and Community Fitness,
the objective pair of [MOGA-Net](moga-net.md), and Newman's modularity $Q$.
Its three letters name those three objectives. The search is the
reference-point-based NSGA-III of Deb and Jain (2014)
([doi:10.1109/TEVC.2013.2281535](https://doi.org/10.1109/TEVC.2013.2281535)),
run over a locus-based genome in which every genome decodes to a valid
partition, so no repair step is needed. [`ccm`][pymocd.ccm] returns the member of the final rank-1 Pareto
front with the highest modularity, and [`ccm_fronts`][pymocd.ccm_fronts]
returns the whole front. The same paper also defines
[NSGA-III-KRM](krm.md), which keeps the engine and the modularity objective but
replaces Community Score and Community Fitness with kernel k-means and ratio
cut.

<div class="pm-facts" markdown>

| | |
|---|---|
| **Entry points** | [`ccm`][pymocd.ccm], [`ccm_fronts`][pymocd.ccm_fronts] |
| **Objectives** | Community Score $\mathrm{CS}$, Community Fitness $\mathrm{CF}$ and modularity $Q$, all three maximized |
| **Search** | NSGA-III (Deb and Jain, 2014) over a locus-based genome, in a private, sequential engine |
| **Selection rule** | Highest modularity $Q$ on the rank-1 front |
| **Front accessor** | [`ccm_fronts`][pymocd.ccm_fronts], a `list[dict[int, int]]` |
| **Deterministic** | No. The engine draws from an OS-seeded generator and the API has no seed parameter (module README) |
| **Original code** | None released (module README) |
| **Paper** | Shaik, Ravi and Deb (2021), *SN Computer Science* 2(1):1 — [doi:10.1007/s42979-020-00382-x](https://doi.org/10.1007/s42979-020-00382-x) |

</div>

## Objectives

Community Score, Community Fitness and the genome come from Pizzuti (2012),
[doi:10.1109/TEVC.2011.2161090](https://doi.org/10.1109/TEVC.2011.2161090), the
paper behind [MOGA-Net](moga-net.md). Let $S$ range over the communities of a
partition, $|S|$ be the number of nodes in $S$, $k_{\mathrm{in}}(i)$ the number
of neighbors of node $i$ that lie in $i$'s own community, $\deg(i)$ the degree
of $i$, $m$ the number of edges, $l_c$ the number of edges inside community $c$
(each counted once) and $d_c$ the sum of the degrees of the nodes of $c$.

$$
\mu_i = \frac{k_{\mathrm{in}}(i)}{|S|}, \qquad
M(S) = \frac{1}{|S|} \sum_{i \in S} \mu_i^{\,r}, \qquad
v_S = \sum_{i \in S} k_{\mathrm{in}}(i)
$$

$$
\mathrm{CS} = \sum_{S} M(S)\, v_S, \qquad
\mathrm{CF} = \sum_{S} \sum_{i \in S} \frac{k_{\mathrm{in}}(i)}{\deg(i)^{\alpha}}, \qquad
Q = \sum_{c} \left[ \frac{l_c}{m} - \left( \frac{d_c}{2m} \right)^{2} \right]
$$

- $v_S$ is twice the number of edges inside $S$.
- $M(S)$ is the plain mean of $\mu_i^{\,r}$, with no outer $1/r$ root. This
  follows Pizzuti's definition. At the default $r = 1$ the two readings
  coincide.
- A node of degree 0 contributes 0 to $\mathrm{CF}$.
- $\mathrm{CF}$ is largest when no edge leaves a community, which on a
  connected graph happens only for the single-community partition, so
  maximizing it favors partitions with few edges between communities.
- $Q$ is the modularity of Newman and Girvan (2004),
  [doi:10.1103/PhysRevE.69.026113](https://doi.org/10.1103/PhysRevE.69.026113).
- The exponents $r$ and $\alpha$ are the parameters `r` and `alpha`. The module
  README states that the paper fixes $r = 1$, which weights all nodes equally,
  and $\alpha = 1$, following MOGA-Net.

The paper maximizes all three objectives. The engine minimizes, so the vector
it receives is $(-\mathrm{CS}, -\mathrm{CF}, -Q)$, and every rank, dominance
test and reference-point distance is computed in that sign-flipped space.

## Search

**Representation.** The genome is Pizzuti's locus-based encoding. It has one
cell per node, in the graph's stable node order. Cell $u$ holds either $u$
itself or the position of one of $u$'s neighbors. Decoding is a union-find
over the cells: each connected component is one community. Every genome that
an operator can produce is valid, so the module has no repair step.
Initialization draws each cell uniformly from $\{u\} \cup N(u)$, where $N(u)$
is the set of neighbors of $u$. Isolated nodes decode to singletons and are
reported as community `-1`.

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
       is their uniform crossover; otherwise it is a clone of one of the two,
       chosen with probability 1/2 each. Every gene of the child is then
       independently resampled from $\{u\} \cup N(u)$ with probability
       `mut_rate`.
    2. Run NSGA-III environmental selection over the $2N$ parents and
       offspring (Deb and Jain, 2014, Algorithms 1 to 4). Whole non-dominated
       fronts are kept while they fit. The remaining places go to members of
       the first front that does not fit: the kept members and that front
       together ($S_t$) are normalized onto a hyperplane, each is associated
       with its nearest reference line, and niching picks the survivors from
       that front. The least crowded reference point wins, with ties broken
       uniformly at random; an empty niche takes its closest member and a
       non-empty one takes a random member.
    3. Apply the paper's two customizations (Sec. 4), in this order: any
       individual whose canonical labeling repeats one already seen is
       replaced by a fresh random genome, and then so is any individual that
       puts every node in one community. Replacements are evaluated at once.
    4. Re-rank the population, because the customizations replaced members
       after selection had ranked them.
4. Sort the final population once more. Its rank-1 members are the result.

The engine is that of Deb and Jain (2014),
[doi:10.1109/TEVC.2013.2281535](https://doi.org/10.1109/TEVC.2013.2281535),
with the fast non-dominated sort of Deb et al. (2002),
[doi:10.1109/4235.996017](https://doi.org/10.1109/4235.996017).

### Randomness

NSGA-III-CCM is not deterministic: two consecutive calls can return different
results. The engine draws from the thread-local generator, which is seeded
from the operating system, and no seed parameter is exposed, so the initial
population, both parent picks, every crossover, clone and mutation draw, the
replacement genomes of the two customizations and the random tie-breaks inside
niching change from run to run. The search is sequential, so the thread count
is not a further source of variation.

The module README reports that 30 runs on Zachary's karate club returned 4
communities every time, attributes the absence of spread to that graph being
degenerate, and reports a rank-1 front of roughly 110 ± 30 members on a
250-node LFR graph.
<!-- TODO(docs): the ccm README also reports "AMI 0.2398 ± 0.0000" over those 30 karate runs without naming the reference partition; against the networkx `club` attribute the current package gives AMI ≈ 0.567 (≈ 0.671 with node 8 in the other faction). Name the ground truth and re-measure before quoting an AMI on this page. -->
See [Performance and reproducibility](../guide/performance.md#nondeterministic-detectors)
for how to handle nondeterministic detectors in experiments.

## Selecting a partition

[`ccm`][pymocd.ccm] returns the rank-1 member with the highest modularity $Q$.
The docstring calls this the paper's recommended decision rule when no ground
truth is available.

[`ccm_fronts`][pymocd.ccm_fronts] returns every rank-1 member as a list of
partitions. The paper reports both the best-NMI and the best-modularity member
of the front (its Tables 1 and 2), and reproducing that protocol needs the
whole candidate set; score each member against the ground truth with
[`nmi`][pymocd.nmi].

Each call runs its own search. Because the search is not deterministic, the
partition returned by `ccm` need not be a member of the front returned by a
separate `ccm_fronts` call. To obtain a front and its selected member from the
same run, call `ccm_fronts` and apply the max-modularity rule yourself, as the
[example](#example) below does. See
[Working with Pareto fronts](../guide/fronts.md) for more on fronts.

## Parameters

<!-- pymocd:defaults pymocd.ccm pymocd.ccm_fronts -->
| Parameter | Default | Description |
|---|---|---|
| `pop_size` | `200` | Population size $N$. The module README requires it to stay above the number of reference points, $H = 91$ at the default `divisions`; the code does not check this. The README notes that Sec. 8.1 of the paper says 20 for karate, while the paper's Table 3 grid and its Fig. S.22 use 200. |
| `num_gens` | `100` | Number of generations. The paper's best-reported combination on Zachary's karate club. |
| `cross_rate` | `0.8` | Probability that a child is the uniform crossover of its two parents rather than a clone of one of them. The paper's best combination. |
| `mut_rate` | `1/68` | Per-gene probability of resampling a cell from $\{u\} \cup N(u)$. The paper's best combination on the 34-node karate graph. |
| `r` | `1.0` | Exponent $r$ applied to each $\mu_i$ in $M(S)$ (there is no outer $1/r$ root, so $M(S)$ is not a power mean). Fixed by the paper at $r = 1$. |
| `alpha` | `1.0` | Exponent $\alpha$ of the degree in Community Fitness. Fixed by the paper at $\alpha = 1$, following MOGA-Net. |
| `divisions` | `12` | Das and Dennis subdivisions $p$; $H = \binom{p + 2}{2}$ reference points for the three objectives, 91 at $p = 12$. The paper does not state a value; the module README takes 12 from pymoo's convention for three objectives. |

Every parameter after `graph` is optional, and [`ccm`][pymocd.ccm] and
[`ccm_fronts`][pymocd.ccm_fronts] take the same ones with the same defaults.

!!! note "The `mut_rate` default is tuned to a 34-node graph"
    The default is the constant $1/68 \approx 0.0147$ (the stub prints it as
    `0.014705882352941176`). The module README gives it only as the paper's
    best combination on the 34-node karate graph, and does not express it as a
    function of the number of nodes $n$. It is not rescaled for other graphs.
    Because mutation acts on each gene independently, a child of an $n$-node
    graph has on average $n \cdot$ `mut_rate` resampled genes, which is $0.5$
    on karate. To use another rate, pass `mut_rate` explicitly.
    <!-- TODO(docs): state whether Shaik, Ravi and Deb (2021) give mut_rate for CCM as a rule in n, and whether its Table 3 has per-dataset rows for CCM like the D1-D4 rows the krm module README lists for KRM; source: the paper (doi:10.1007/s42979-020-00382-x), Table 3. -->

??? info "Internal constants"
    These are fixed in `src/core/algorithms/ccm/nsga3/normalize.rs` and are not
    exposed: `NEAR_ZERO = 1e-10` is the threshold below which a pivot, an
    intercept or a fallback scale counts as zero; `MIN_INTERCEPT = 1e-6` is the
    smallest accepted hyperplane intercept, below which a fallback scaling is
    used; and `ASF_OFF_AXIS_WEIGHT = 1e-6` is the off-axis weight of the
    achievement scalarizing function that picks the extreme points. In
    `nsga3/niching.rs`, a reference direction whose squared norm is below
    `1e-30` is treated as the origin.

### Degenerate inputs

The budget arguments are not validated. Every row was checked against the
current package.

| Input | Result |
|---|---|
| `pop_size=0` | `ccm` panics with "empty Pareto front"; `ccm_fronts` returns `[]`. |
| `cross_rate` or `mut_rate` outside $[0, 1]$ | Both functions panic with `p=… is outside range [0.0, 1.0]`. |
| An empty graph | `ccm` returns `{}`; `ccm_fronts` returns `pop_size` empty dicts. |
| A graph with no edges | `ccm` maps every node to `-1`; `ccm_fronts` returns `pop_size` copies of that partition. |

A panic reaches Python as `pyo3_runtime.PanicException`, which derives from
`BaseException`, so `except Exception` does not catch it.

## Complexity

The module README does not state a time or memory complexity bound.
<!-- TODO(docs): add a complexity bound for NSGA-III-CCM per generation and per run; sources: the paper (doi:10.1007/s42979-020-00382-x) or a derivation by the maintainer from src/core/algorithms/ccm/. -->
What the code does per run: it evaluates the $N$ initial genomes, then in each
of `num_gens` generations it evaluates $N$ offspring plus every replacement
made by the two customizations, and runs environmental selection over $2N$
individuals. Each evaluation decodes a genome by union-find and computes the
three objectives in passes over every node's neighbor list. The search runs
on one thread.

## When to use it

- **Reproducing or comparing with Shaik, Ravi and Deb (2021).** The
  `num_gens`, `cross_rate` and `mut_rate` defaults are the paper's
  best-reported combination on Zachary's karate club, and
  [`ccm_fronts`][pymocd.ccm_fronts] exposes the front for the paper's
  best-NMI and best-modularity protocol.
- **Comparing with MOGA-Net.** CCM keeps MOGA-Net's two objectives, adds
  modularity as a third, and replaces MOGA-Net's NSGA-II-ranked genetic
  algorithm with NSGA-III. Both detectors
  return the max-modularity member of their rank-1 front. Note that
  [`moga_net`][pymocd.moga_net] defaults to `r=2.0` while `ccm` defaults to
  `r=1.0`.
- **As a baseline, not as a default.** NSGA-III-CCM is one of the library's
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
  problems, which the paper adapted; there is no binary tournament.
- $M(S)$ is the plain mean of $\mu_i^{\,r}$ with no outer $1/r$ root, following
  Pizzuti's definition. This is not a dropped operation.
- `pop_size = 200` resolves a contradiction in the paper and `divisions = 12`
  an omission; see [Parameters](#parameters).
- Modularity is computed inside the module rather than by the library's shared
  modularity function, so the inner loop never builds a partition object. A
  test asserts that the two agree to within $10^{-12}$.
- The public functions return normalized partitions in which isolated nodes
  are community `-1`; the paper leaves them as singletons.
- The duplicate filter compares canonical labelings, so it also removes
  partitions that differ only by a permutation of community labels.

The code has one more difference that the CCM README does not list but the
KRM README records for the same engine design (see
[NSGA-III-KRM](krm.md#divergences-from-the-paper)): the ideal point and the
extreme points are recomputed from $S_t$ at every normalization
(`nsga3/normalize.rs`), whereas pymoo carries them across generations.

## Example

The block below runs [`ccm`][pymocd.ccm] and [`ccm_fronts`][pymocd.ccm_fronts]
on Zachary's karate club, then applies the max-modularity rule to the front by
hand. The two calls are independent runs. The modularity is computed with
`weight=None` because `nx.karate_club_graph()` carries edge weights, while the
modularity inside CCM counts edges.

```python exec="on" source="above" session="ccm" result="text"
import networkx as nx
import pymocd

G = nx.karate_club_graph()

part = pymocd.ccm(G)        # max-Q member of one run's rank-1 front
print("ccm:", len(set(part.values())), "communities")

front = pymocd.ccm_fronts(G)  # the rank-1 front of a second, independent run
print("ccm_fronts:", len(front), "partitions on the front")


def communities(partition):
    groups = {}
    for node, label in partition.items():
        groups.setdefault(label, set()).add(node)
    return list(groups.values())


# ccm's decision rule, applied to this front
scored = [
    (nx.community.modularity(G, communities(p), weight=None), p) for p in front
]
q_best, best = max(scored, key=lambda item: item[0])
sizes = sorted(len(set(p.values())) for p in front)
print(f"communities per front member: {sizes[0]} to {sizes[-1]}")
print(f"max-Q member: {len(set(best.values()))} communities, Q = {q_best:.4f}")
```

The first line is the number of communities `ccm` selected. The last line is
the member the same rule picks from the second run's front; the range above it
shows how far the front spreads in the number of communities.

## Reference

Tanveerulhuq Shaik, Vadlamani Ravi and Kalyanmoy Deb. "Evolutionary
Multi-Objective Optimization Algorithm for Community Detection in Complex
Social Networks." *SN Computer Science* 2(1), article 1, 2021 (published online
7 November 2020).
[doi:10.1007/s42979-020-00382-x](https://doi.org/10.1007/s42979-020-00382-x).

The volume, issue, article number and dates follow the Crossref record of the
DOI.
<!-- TODO(docs): the module READMEs src/core/algorithms/ccm/README.md, src/core/algorithms/krm/README.md and src/core/algorithms/README.md cite "SN Computer Science 2:13"; Crossref gives article number 1 for this DOI, and article 13 of SN Comput. Sci. 2(1) is doi:10.1007/s42979-020-00407-5. Confirm on the publisher's page and fix the READMEs. -->

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

## See also

- [NSGA-III-KRM](krm.md): the same paper's second variant, with kernel
  k-means, ratio cut and modularity.
- [MOGA-Net](moga-net.md): Community Score and Community Fitness under
  NSGA-II ranking with elitist generational replacement.
- API reference: [`ccm`][pymocd.ccm], [`ccm_fronts`][pymocd.ccm_fronts],
  [`nmi`][pymocd.nmi].
- Guide: [Core concepts](../guide/concepts.md),
  [Working with Pareto fronts](../guide/fronts.md),
  [Evaluating partitions](../guide/evaluation.md),
  [Choosing a detector](../guide/choosing.md).
- Gallery: [Every detector on one graph](../examples/basics/plot_detector_tour.md), which
  runs `ccm` alongside the other detectors, and
  [Choosing a partition from a front](../examples/fronts/plot_choosing_from_front.md).
- [Algorithms overview](index.md) and the
  [module README](https://github.com/csilab-ufop/pymocd/blob/master/src/core/algorithms/ccm/README.md),
  which carries the full derivation and file map.
