---
description: MOGA-Net in pymocd, a re-implementation of Pizzuti (2009, 2012) that maximizes community score and community fitness over a locus-based genome and returns the max-modularity member of its Pareto front.
---

# MOGA-Net

<span class="pm-pill pm-pill--baseline">Re-implemented baseline</span>
<span class="pm-pill pm-pill--neutral">Multi-objective</span>

MOGA-Net, from Pizzuti (2012)
([doi:10.1109/TEVC.2011.2161090](https://doi.org/10.1109/TEVC.2011.2161090)),
first presented at ICTAI 2009
([doi:10.1109/ICTAI.2009.58](https://doi.org/10.1109/ICTAI.2009.58)),
maximizes two objectives together: the community score, which favours dense
communities, and the community fitness, which grows as more of each node's
edges stay inside its own community. The search ranks a population of
locus-based genomes with the non-dominated sorting and crowding distance of
NSGA-II; every gene links a node to one of its neighbours by construction, so
no repair step is needed. [`moga_net`][pymocd.moga_net] returns the member of
the final rank-1 Pareto front with the highest modularity, the decision rule of
Pizzuti (2012, Sec. V-E), and [`moga_net_fronts`][pymocd.moga_net_fronts]
returns the whole front. The port's objective values match a result file
shipped with Pizzuti's code to the 7 significant digits it prints; its
generational loop is the part that cannot be checked against her code, which is
distributed as obfuscated MATLAB (see
[Divergences from the paper](#divergences-from-the-paper)).

<div class="pm-facts" markdown>

| | |
|---|---|
| **Entry points** | [`moga_net`][pymocd.moga_net], [`moga_net_fronts`][pymocd.moga_net_fronts] |
| **Objectives** | Community score $\mathrm{CS}$ and community fitness $\mathrm{CF}$, both maximized |
| **Search** | Non-dominated sorting and crowding distance of NSGA-II, with elitist generational replacement and roulette selection, over a locus-based genome; private, single-threaded engine |
| **Selection rule** | Highest modularity $Q$ on the rank-1 front (Pizzuti 2012, Sec. V-E) |
| **Front accessor** | [`moga_net_fronts`][pymocd.moga_net_fronts], a `list[dict[int, int]]` |
| **Deterministic** | No. The engine draws from the thread-local, OS-seeded generator and the API has no seed parameter (module README) |
| **Original code** | [Moganet2016.zip](https://staff.icar.cnr.it/pizzuti/codice/Moganet2016.zip), MATLAB shipped as obfuscated P-code |
| **Paper** | Pizzuti (2012), *IEEE Transactions on Evolutionary Computation* 16(3):418–430 — [doi:10.1109/TEVC.2011.2161090](https://doi.org/10.1109/TEVC.2011.2161090); first presented at ICTAI 2009 — [doi:10.1109/ICTAI.2009.58](https://doi.org/10.1109/ICTAI.2009.58) |

</div>

The module README treats the 2012 journal paper as the authority for
parameters and the 2009 conference paper as the one that defines the genetic
operators.

## Objectives

Let $S$ range over the communities of a partition, $|S|$ be the number of
nodes in $S$, $k_{\mathrm{in}}(i)$ the number of neighbours of node $i$ that
lie in $i$'s own community, and $\deg(i)$ the degree of $i$. The exponents $r$
and $\alpha$ are the parameters `r` and `alpha`. Following Pizzuti (2012,
Sec. V-A) and the module README:

$$
\mu_i = \frac{k_{\mathrm{in}}(i)}{|S|}, \qquad
M(S) = \frac{1}{|S|} \sum_{i \in S} \mu_i^{\,r}, \qquad
v_S = \sum_{i \in S} k_{\mathrm{in}}(i)
$$

$$
\mathrm{CS} = \sum_{S} M(S)\, v_S, \qquad
\mathrm{CF} = \sum_{S} \sum_{i \in S} \frac{k_{\mathrm{in}}(i)}{\deg(i)^{\alpha}}
$$

- $\mu_i$ lies in $[0, 1)$, because a node has at most $|S| - 1$ neighbours
  inside $S$.
- $v_S$ is twice the number of edges inside $S$.
- $M(S)$ is the plain mean of $\mu_i^{\,r}$, with no outer $1/r$ root.
- A node of degree 0 contributes 0 to $\mathrm{CF}$.
- The community fitness goes back to Lancichinetti, Fortunato and Kertész
  (2009),
  [doi:10.1088/1367-2630/11/3/033015](https://doi.org/10.1088/1367-2630/11/3/033015);
  MOGA-Net uses the per-node form written above.

Both objectives are maximized. The engine's dominance test minimizes, so each
individual stores $(-\mathrm{CS}, -\mathrm{CF})$. The module README warns that
treating $\mathrm{CF}$ as a quantity to minimize is a sign error that
over-fragments.

The two objectives pull in opposite directions. The README describes
$\mathrm{CS}$ as rising with cohesive but not-too-large communities, and
$\mathrm{CF}$ as rising monotonically with coarseness: merging two communities
never lowers any $k_{\mathrm{in}}(i)$, so a single community maximizes
$\mathrm{CF}$ and all-singletons gives 0. The module test
`cs_maximal_for_two_community_split` checks this on two triangles joined by
one edge: the two-triangle split has a larger $\mathrm{CS}$ than both the
single community and the all-singletons partition, and the single community
has a larger $\mathrm{CF}$ than the split.

The README reports that the port reproduces all three numeric columns of
`karate.pairs.ris`, a result file shipped with Pizzuti's Moganet2016 package,
on all 7 of her solutions, to 7 significant digits. Read against the
solutions in the companion file `karate.pairs.sol`, the three columns are
$-\mathrm{CS}$ at $r = 2$, $34 - \mathrm{CF}$ at $\alpha = 1$ (34 is the
number of nodes) and the modularity $Q$. The objectives are therefore not
among the port's divergences.

The definitions translate directly into a few lines of Python. The block below
evaluates them at the defaults $r = 2$ and $\alpha = 1$ for three partitions of
Zachary's karate club: one community, the two factions stored in networkx's
`club` node attribute, and all singletons.

```python exec="on" source="above" session="moga-net" result="text"
import networkx as nx


def community_score_fitness(G, partition, r=2.0, alpha=1.0):
    """CS and CF of a {node: community} dict, as defined above."""
    groups = {}
    for node, label in partition.items():
        groups.setdefault(label, set()).add(node)
    cs = cf = 0.0
    for S in groups.values():
        k_in = {i: sum(1 for j in G[i] if j in S) for i in S}
        M = sum((k_in[i] / len(S)) ** r for i in S) / len(S)
        cs += M * sum(k_in.values())
        cf += sum(k_in[i] / G.degree(i) ** alpha for i in S if G.degree(i) > 0)
    return cs, cf


G = nx.karate_club_graph()
candidates = {
    "one community": {v: 0 for v in G},
    "two factions": {v: int(G.nodes[v]["club"] == "Officer") for v in G},
    "singletons": {v: v for v in G},
}
for name, p in candidates.items():
    cs, cf = community_score_fitness(G, p)
    print(f"{name:14s} CS = {cs:7.3f}   CF = {cf:7.3f}")
print("CF bound at alpha = 1 (non-isolated nodes):", sum(1 for v in G if G.degree(v) > 0))
```

The single community reaches the largest $\mathrm{CF}$, equal to the bound
explained under [Parameters](#parameters); the singletons score 0 on both
objectives; the two-faction split gives up some $\mathrm{CF}$ for a larger
$\mathrm{CS}$.

## Search

**Representation.** The genome is the locus-based adjacency encoding of Park
and Song (1998), "A genetic algorithm for clustering problems", *Proc. 3rd
Annual Conference on Genetic Programming*, pp. 568–575 (Crossref lists no
DOI for it). For a graph of $n$ nodes it has $n$ genes, one per node in the
graph's sorted node order, and gene $p$ holds a position $q$, read as "link
node $p$ to node $q$". Decoding is a union-find over the $n$ links: each group
of linked nodes is one community. The decoder labels the groups in order of
first appearance by position, so two genomes that encode the same partition
decode to the same internal label array.

**Initialization.** The port uses the paper's safe variant: every gene is drawn
uniformly from the neighbours of its own node. Both operators keep that
property (crossover copies genes from safe parents, and mutation redraws from
neighbours), so the module has no repair pass. Because every link follows an
edge of the network, every community induces a connected subgraph. An isolated
node has no neighbour to point at and keeps the self-link $p \to p$; it becomes
a singleton and is reported as community `-1`.

**Loop.** With `pop_size` $= N$:

1. Draw $N$ random safe genomes, decode them and evaluate each into
   $(-\mathrm{CS}, -\mathrm{CF})$.
2. For each of `num_gens` generations:
    1. Run the fast non-dominated sort and the crowding distance over the whole
       population. The population size stays exactly $N$; there is never a
       combined pool of $2N$ parents and offspring.
    2. Order the members by rank (ascending), then crowding distance
       (descending).
    3. Copy the first $\max(1, \operatorname{round}(0.1\,N))$ members of that
       order unchanged into the next generation (30 at the default $N = 300$).
    4. Fill the remaining places with children. Two parents are drawn by
       roulette over the order, the member at position $j$ ($j = 1$ for the
       best) having weight $1/\sqrt{j}$. With probability `cross_rate` the
       child is their uniform crossover, which takes each gene from either
       parent with probability 1/2; otherwise it is a clone of one of the two,
       chosen by a fair coin. Every child is then mutated: each gene, with
       probability `mut_rate`, is redrawn uniformly from its node's
       neighbours. The child is evaluated.
    5. The elite copies and the children replace the population.
3. Rank the final population once more. Its rank-1 members are the result.

The non-dominated sort and crowding distance follow Deb et al. (2002),
[doi:10.1109/4235.996017](https://doi.org/10.1109/4235.996017). They are
re-implemented inside the module, single-threaded, rather than taken from the
library's shared NSGA-II machinery, so that this baseline's runtime tracks the
published single-threaded method (module README).

### Randomness

MOGA-Net is not deterministic: two consecutive calls can return different
results. The engine draws from the thread-local generator, which is seeded
from the operating system and cannot be re-seeded through the public API;
there is no seed parameter. Everything downstream of the random draws is
order-fixed, so the module README identifies the random stream as the only
source of run-to-run variation. The search is sequential, so the thread count
plays no part.

Calls can return different partitions even on a graph as small as Zachary's
karate club.
<!-- TODO(docs): the module README's "Measured reference on karate over 30 runs: AMI 0.2554 ± 0.0000, k = 4.00" was dropped because it does not reproduce with pymocd 6.0.1 (jackal, 30 calls of moga_net on nx.karate_club_graph()): 2 distinct partitions, k = 4 every time, AMI 0.474 ± 0.018 against networkx's 'club' labels and 0.584 ± 0.018 with node 8 moved to the Officers; the README does not name its ground truth. Maintainer: re-measure and update the README. -->
See [Performance and reproducibility](../guide/performance.md) for how to
handle nondeterministic detectors in experiments.

## Selecting a partition

[`moga_net`][pymocd.moga_net] returns the rank-1 member with the highest
modularity
$Q = \sum_c \left[ l_c/m - \left( d_c / 2m \right)^{2} \right]$
of Newman and Girvan (2004),
[doi:10.1103/PhysRevE.69.026113](https://doi.org/10.1103/PhysRevE.69.026113),
where $m$ is the number of edges, $l_c$ the number of edges inside community
$c$ and $d_c$ the sum of the degrees of its nodes. This is the reduction to one
partition described in Pizzuti (2012, Sec. V-E). $Q$ is computed by counting
edges; edge weights are not used.

[`moga_net_fronts`][pymocd.moga_net_fronts] returns every rank-1 member as a
list of partitions. Pizzuti's Table 1 reports the best-NMI member of the front
instead, and the docstring notes that reproducing that protocol needs the whole
candidate set; score each member against the ground truth with
[`nmi`][pymocd.nmi]. The list can contain the same partition more than once:
the locus-based encoding maps many genomes to one partition, and nothing in the
loop removes duplicates. Count distinct members if that matters, as the
[example](#example) does.

Each call runs its own search. Because the search is not deterministic, the
partition returned by `moga_net` need not be a member of the front returned by
a separate `moga_net_fronts` call. To obtain a front and its selected member
from the same run, call `moga_net_fronts` and apply the max-modularity rule
yourself. See [Working with Pareto fronts](../guide/fronts.md) for more on
fronts.

## Parameters

<!-- pymocd:defaults pymocd.moga_net pymocd.moga_net_fronts -->
| Parameter | Default | Description |
|---|---|---|
| `pop_size` | `300` | Population size $N$, constant across generations. Pizzuti's setting. |
| `num_gens` | `30` | Number of generations. Pizzuti's setting. |
| `cross_rate` | `0.8` | Probability that a child is the uniform crossover of its two parents rather than a clone of one of them. Pizzuti's setting. |
| `mut_rate` | `0.2` | Per-gene probability of redrawing the gene uniformly from the node's neighbours; applied to every child. Pizzuti's setting. |
| `r` | `2.0` | Exponent $r$ in $M(S)$, inside the community score. Fixed at 2 by Pizzuti (2012, Sec. VI-C). |
| `alpha` | `1.0` | Exponent $\alpha$ on $\deg(i)$ in the community fitness. Pizzuti's default. |

[`moga_net`][pymocd.moga_net] and
[`moga_net_fronts`][pymocd.moga_net_fronts] take the same parameters with the
same defaults. The module README takes `pop_size`, `num_gens`, `cross_rate`,
`mut_rate` and the 10% elite fraction from Pizzuti's own single-run settings,
which it quotes as: *"The population size was 300, the number of generations
30, the crossover rate 0.8, the mutation rate 0.2, elite reproduction 10% of
the population size, roulette selection function."* The driver script
`mainMOGANet.m` in Moganet2016, which is readable MATLAB, keeps
`mutationRate=0.2` and `CrossoverFraction=0.8` but sets `PopSize=30` and
`maxGen=50`, and its argument `alfa` is the community-score exponent $r$, not
the community-fitness exponent that pymocd calls `alpha`. Because `mut_rate` acts
on each gene independently, a child of an $n$-node graph has on average
$0.2\,n$ genes redrawn at the default; the rate is not rescaled with $n$.

### The exponents `r` and `alpha`

`r`
:   The exponent applied to each node's internal-degree fraction in the
    Community Score; $M(S)$ takes no outer $1/r$ root. The docstring of
    [`moga_net`][pymocd.moga_net] describes it as a resolution knob, adding
    that higher values help at high mixing. The module README states that a
    higher `r` penalizes nodes weakly attached to their community. TEVC 2012
    Sec. VI-C fixes `r` at 2, which is the default.

`alpha`
:   The Community Fitness exponent. The docstrings stress that it does **not**
    set a community size. In the per-node form used here,

    $$
    \mathrm{CF} \le \sum_{i:\ \deg(i) > 0} \deg(i)^{1-\alpha}
    \quad \text{for every } \alpha,
    $$

    because $k_{\mathrm{in}}(i) \le \deg(i)$ for every node; the bound is
    reached exactly when no edge leaves a community, which on a connected
    graph means the single-community partition. What `alpha` changes is who
    counts: $\alpha > 1$ discounts high-degree nodes, so the internal edges of
    low-degree nodes matter relatively more. At $\alpha = 1$, Pizzuti's
    default, each term is the fraction of node $i$'s edges that stay inside
    its community, and the bound is the number of non-isolated nodes.

!!! note "How the default `r = 2` was settled"
    The value comes from TEVC 2012 Sec. VI-C, which the README quotes as
    *"the parameter r ... has been set to 2"*. The README adds a numerical
    check against the Moganet2016 package: fitting `r` to the numeric columns
    of Pizzuti's `karate.pairs.ris` gives a unique minimum at `r = 2`. It
    warns against re-deriving `r` from the karate NMI reported in the paper,
    because Pizzuti's own binary does not reproduce that number. The README
    also records that the default was moved from an earlier 1.5 to 2, and
    that 2 is kept "because the point of this module is to reproduce the
    published method, not to win".

??? info "Internal constant"
    The elite fraction, `ELITE_FRACTION = 0.10` in `search/engine.rs`, is
    Pizzuti's "elite reproduction 10% of the population size". It is not a
    public parameter.

### Degenerate inputs

The budget arguments are not validated. Every row was checked against the
current package.

| Input | Result |
|---|---|
| `pop_size=0` | `moga_net` panics with "empty Pareto front"; `moga_net_fronts` returns `[]`. |
| `cross_rate` or `mut_rate` outside $[0, 1]$ | Both functions panic with `p=… is outside range [0.0, 1.0]`. |
| An empty graph | `moga_net` returns `{}`; `moga_net_fronts` returns `pop_size` empty dicts. |
| A graph with no edges | `moga_net` maps every node to `-1`; `moga_net_fronts` returns `pop_size` copies of that partition. |

A panic reaches Python as `pyo3_runtime.PanicException`, which derives from
`BaseException`, so `except Exception` does not catch it.

## Complexity

Neither the module README nor the docstrings state a time or memory complexity
bound.
<!-- TODO(docs): add a complexity bound for MOGA-Net per generation and per run; sources: Pizzuti (2012, doi:10.1109/TEVC.2011.2161090), Pizzuti (2009, doi:10.1109/ICTAI.2009.58), or a derivation by the maintainer from src/core/algorithms/moganet/. -->
What the code does per run: it evaluates the $N$ initial genomes, then in each
of `num_gens` generations it runs one non-dominated sort that compares every
pair of the $N$ members, computes the crowding distance by sorting each front
on each of the two objectives, and evaluates the $N - \max(1,
\operatorname{round}(0.1\,N))$ new children. Each evaluation decodes a genome
by union-find over $n$ links and computes both objectives in one pass over
every node's neighbour list. After the last generation the population is
ranked once more, and `moga_net` computes $Q$ for each rank-1 member. The
search runs on one thread.

## When to use it

- **Reproducing or comparing with Pizzuti (2012).** The defaults are the
  settings the module README quotes from Pizzuti, `r` follows Sec. VI-C, the
  objectives are checked against a
  result file shipped with her code, and
  [`moga_net_fronts`][pymocd.moga_net_fronts] exposes the
  front for the best-NMI protocol of her Table 1.
- **Comparing with NSGA-III-CCM.** [NSGA-III-CCM](ccm.md) keeps MOGA-Net's two
  objectives, adds modularity as a third, and searches with NSGA-III. Both
  detectors return the max-modularity member of their rank-1 front. Note that
  [`moga_net`][pymocd.moga_net] defaults to `r=2.0` while
  [`ccm`][pymocd.ccm] defaults to `r=1.0`.
- **As a baseline, not as a default.** MOGA-Net is one of the library's
  re-implemented baselines, which exist for papers and benchmarks; the module
  README states that its purpose is to reproduce the published method. Its
  single-threaded engine tracks the cost of the published method rather than
  the library's optimizations. For general use, start from
  [Choosing a detector](../guide/choosing.md); [RIMPSO](rimpso.md) is the
  library's recommended default.

## Divergences from the paper

The module README lists these. Pizzuti's search loop is distributed as
obfuscated P-code (`MOGANet.p`, `create_population.p`, `crossover_cluster.p`,
`mutate_cluster.p`), which the README notes cannot be read, so the first three
items, which concern that loop, cannot be settled against her source. The
README marks the first two as unverifiable, and they are also recorded in the
header of `search/engine.rs`.

- **Replacement scheme (unverifiable).** The paper describes NSGA-II, whose
  replacement merges parents and offspring into a pool of $2N$ and truncates
  by rank and crowding distance. This loop instead copies the top 10% and
  fills the rest with fresh children, so a child never competes against the
  parent generation. The README identifies this as the MATLAB
  `gamultiobj`-style elitist and roulette generational replacement of the
  paper's own tooling.
- **Mutation applied to every child (unverifiable).** The paper applies
  mutation to crossover offspring; here every child is mutated, including one
  that was cloned rather than crossed. The README notes that
  `crossoverFraction` in `gamultiobj` can be read either way and keeps the
  composed reading.
- **Selection fitness.** The paper only says "roulette selection function".
  The weights used are `gamultiobj`'s default rank scaling, $1/\sqrt{j}$ over
  the rank-then-crowding order.
- **No shared engine, no thread parallelism.** The non-dominated sort and the
  crowding distance are re-implemented in the module, single-threaded, rather
  than reusing the library's shared NSGA-II and NSGA-III machinery.
- **No repair operator.** The paper describes repairing an invalid gene to one
  of the node's neighbours. Safe initialization and neighbour-respecting
  operators make every gene valid by construction, so the repair pass is
  absent rather than implemented and never triggered.
- **Decision rule.** Reducing the rank-1 front to one partition by maximum
  modularity follows TEVC 2012 Sec. V-E. Pizzuti's Table 1 protocol instead
  reports the best-NMI member of the front, which is what
  [`moga_net_fronts`][pymocd.moga_net_fronts] exposes.

The README states that the objectives are not a divergence: they reproduce
the values in Pizzuti's `karate.pairs.ris` result file to 7 significant
digits (see [Objectives](#objectives)).

## Example

The block below continues the session of the block under
[Objectives](#objectives), reusing `G`, the two-faction ground truth and
`community_score_fitness`. It runs [`moga_net`][pymocd.moga_net] and
[`moga_net_fronts`][pymocd.moga_net_fronts] on Zachary's karate club, then
applies both selection rules to the front by hand: the max-modularity rule of
`moga_net`, and the best-NMI rule of Pizzuti's Table 1, which needs the ground
truth. The two calls are independent runs. The modularity is computed with
`weight=None` because `nx.karate_club_graph()` carries edge weights, while
MOGA-Net counts edges.

```python exec="on" source="above" session="moga-net" result="text"
import pymocd

part = pymocd.moga_net(G)          # max-Q member of one run's rank-1 front
print("moga_net:", len(set(part.values())), "communities")

front = pymocd.moga_net_fronts(G)  # rank-1 front of a second, independent run
distinct = {tuple(sorted(p.items())) for p in front}
print(f"moga_net_fronts: {len(front)} members, {len(distinct)} distinct partitions")


def modularity(partition):
    groups = {}
    for node, label in partition.items():
        groups.setdefault(label, set()).add(node)
    return nx.community.modularity(G, list(groups.values()), weight=None)


factions = candidates["two factions"]
max_q = max(front, key=modularity)                            # Pizzuti Sec. V-E
best_nmi = max(front, key=lambda p: pymocd.nmi(p, factions))  # Pizzuti Table 1
for name, p in [("max-Q member", max_q), ("best-NMI member", best_nmi)]:
    cs, cf = community_score_fitness(G, p)
    print(f"{name:15s} k = {len(set(p.values()))}  Q = {modularity(p):.4f}  "
          f"NMI = {pymocd.nmi(p, factions):.4f}  CS = {cs:.3f}  CF = {cf:.3f}")
```

The first line is the number of communities `moga_net` selected. The last two
lines compare the two rules on the same front: the member `moga_net`'s rule
would pick, and the member Pizzuti's Table 1 protocol would report. Both lie
on the rank-1 front, so neither dominates the other in
$(\mathrm{CS}, \mathrm{CF})$: when the rules pick different partitions, one
has the larger $\mathrm{CS}$ and the other the larger $\mathrm{CF}$, unless
their scores tie.

## Reference

Clara Pizzuti. "A Multiobjective Genetic Algorithm to Find Communities in
Complex Networks." *IEEE Transactions on Evolutionary Computation* 16(3):
418–430, June 2012.
[doi:10.1109/TEVC.2011.2161090](https://doi.org/10.1109/TEVC.2011.2161090).

```bibtex
@article{Pizzuti_2012,
  title     = {A Multiobjective Genetic Algorithm to Find Communities in Complex Networks},
  volume    = {16},
  ISSN      = {1941-0026},
  url       = {http://dx.doi.org/10.1109/TEVC.2011.2161090},
  DOI       = {10.1109/tevc.2011.2161090},
  number    = {3},
  journal   = {IEEE Transactions on Evolutionary Computation},
  publisher = {Institute of Electrical and Electronics Engineers (IEEE)},
  author    = {Pizzuti, Clara},
  year      = {2012},
  month     = jun,
  pages     = {418–430}
}
```

Clara Pizzuti. "A Multi-objective Genetic Algorithm for Community Detection in
Networks." In *2009 21st IEEE International Conference on Tools with
Artificial Intelligence (ICTAI)*, Newark, New Jersey, USA, pp. 379–386,
November 2009.
[doi:10.1109/ICTAI.2009.58](https://doi.org/10.1109/ICTAI.2009.58).

```bibtex
@inproceedings{Pizzuti_2009,
  title     = {A Multi-objective Genetic Algorithm for Community Detection in Networks},
  url       = {http://dx.doi.org/10.1109/ICTAI.2009.58},
  DOI       = {10.1109/ictai.2009.58},
  booktitle = {2009 21st IEEE International Conference on Tools with Artificial Intelligence},
  publisher = {IEEE},
  author    = {Pizzuti, Clara},
  year      = {2009},
  month     = Nov,
  pages     = {379–386}
}
```

Both entries are the publisher's records from doi.org. The only change is the
month of the journal entry: doi.org writes `month=June`, which is not a BibTeX
month macro, so it is given here as the standard macro `jun`.

Works the implementation builds on, as cited in the module README:

- Y. Park and M. Song. "A genetic algorithm for clustering problems." In
  *Proceedings of the 3rd Annual Conference on Genetic Programming*,
  pp. 568–575, 1998 (the locus-based representation; Crossref lists no DOI).
- K. Deb, A. Pratap, S. Agarwal and T. Meyarivan. "A fast and
  elitist multiobjective genetic algorithm: NSGA-II." *IEEE Transactions on
  Evolutionary Computation* 6(2): 182–197, 2002.
  [doi:10.1109/4235.996017](https://doi.org/10.1109/4235.996017) (non-dominated
  sort and crowding distance).
- Andrea Lancichinetti, Santo Fortunato and János Kertész. "Detecting the
  overlapping and hierarchical community structure in complex networks."
  *New Journal of Physics* 11(3): 033015, 2009.
  [doi:10.1088/1367-2630/11/3/033015](https://doi.org/10.1088/1367-2630/11/3/033015)
  (the community fitness).

The modularity of the decision rule is that of M. E. J. Newman and M. Girvan,
"Finding and evaluating community structure in networks," *Physical Review E*
69(2): 026113, 2004,
[doi:10.1103/PhysRevE.69.026113](https://doi.org/10.1103/PhysRevE.69.026113).

## See also

- [NSGA-III-CCM](ccm.md): community score and community fitness, plus
  modularity, under NSGA-III.
- API reference: [`moga_net`][pymocd.moga_net],
  [`moga_net_fronts`][pymocd.moga_net_fronts], [`nmi`][pymocd.nmi].
- Guide: [Core concepts](../guide/concepts.md),
  [Working with Pareto fronts](../guide/fronts.md),
  [Evaluating partitions](../guide/evaluation.md),
  [Choosing a detector](../guide/choosing.md).
- Gallery: [Choosing a partition from a front](../examples/fronts/plot_choosing_from_front.md).
- [Algorithms overview](index.md) and the
  [module README](https://github.com/csilab-ufop/pymocd/blob/master/src/core/algorithms/moganet/README.md),
  which carries the full derivation and file map.
