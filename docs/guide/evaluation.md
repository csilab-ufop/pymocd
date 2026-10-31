---
description: How pymocd scores a partition against ground truth with NMI, AMI, ARI and pairwise F1, the exact definitions implemented, and how missing and isolated nodes are handled.
---

# Evaluating partitions

pymocd ships four agreement scores between a detected partition and a reference
partition (usually ground truth): normalized mutual information (NMI), adjusted
mutual information (AMI), the adjusted Rand index (ARI) and a pairwise F1 score.
They are computed in Rust from the contingency table of the two labelings
and need no other package. This page gives the definition of each score exactly as implemented,
explains which nodes are scored and how the label `-1` of isolated nodes is
treated, and lists the conventions applied to degenerate inputs. The scores
compare two labelings only; the graph plays no part in them.

## The functions

| Function | Returns |
|---|---|
| [`gt_metrics(partition, gt)`][pymocd.gt_metrics] | the tuple `(NMI, AMI, ARI, F1)` |
| [`nmi(partition, gt)`][pymocd.nmi] | NMI as a `float` |
| [`ami(partition, gt)`][pymocd.ami] | AMI as a `float` |
| [`ari(partition, gt)`][pymocd.ari] | ARI as a `float` |
| [`f1(partition, gt)`][pymocd.f1] | pairwise F1 as a `float` |

Both arguments are `{node: community}` dicts, the shape the detectors return.
Node ids and community labels must be integers (Python `int` or NumPy integer
scalars) that fit in a signed 32-bit integer. A string or `float` raises
`TypeError` and a larger integer raises `OverflowError`, so ground truth stored
as names has to be mapped to integers first. Labels are identifiers only:
renaming the communities of either argument changes no score beyond
floating-point rounding. All four formulas below are symmetric in their two
arguments.

A ground-truth example on Zachary's karate club, whose NetworkX copy stores
each member's faction in the `club` node attribute:

```python exec="on" source="above" session="evaluation" result="text"
import networkx as nx
import pymocd

G = nx.karate_club_graph()
gt = {v: 0 if G.nodes[v]["club"] == "Mr. Hi" else 1 for v in G}

partition = pymocd.rimpso(G)
nmi, ami, ari, f1 = pymocd.gt_metrics(partition, gt)
print(f"communities found: {len(set(partition.values()))}")
print(f"gt_metrics:  NMI={nmi:.4f}  AMI={ami:.4f}  ARI={ari:.4f}  F1={f1:.4f}")
print(f"one by one:  NMI={pymocd.nmi(partition, gt):.4f}  AMI={pymocd.ami(partition, gt):.4f}"
      f"  ARI={pymocd.ari(partition, gt):.4f}  F1={pymocd.f1(partition, gt):.4f}")
```

## Which nodes are scored

Only nodes that appear as keys in **both** dicts are scored. A node missing from
either dict is skipped without a warning, so a partition that omits nodes is
compared on the shared nodes and is not penalized for the omission. Check
coverage yourself when it matters:

```python exec="on" source="above" session="evaluation" result="text"
first_half = {v: c for v, c in partition.items() if v < 17}
print(f"scored nodes: {len(first_half)} of {len(gt)}")
print(f"AMI on the shared nodes only: {pymocd.ami(first_half, gt):.4f}")
print("same node set:", first_half.keys() == gt.keys())
```

## Definitions

Write $U$ for the first argument, $V$ for the second, and $n$ for the number of
scored nodes. Let $n_{ij}$ be the number of scored nodes with label $i$ in $U$
and label $j$ in $V$, and let $a_i = \sum_j n_{ij}$ and $b_j = \sum_i n_{ij}$ be
the community sizes. Sums over cells run over the non-empty cells. Logarithms
are natural (nats); NMI and AMI are ratios, so the base cancels.

**Entropy and mutual information.**

$$
H(U) = -\sum_i \frac{a_i}{n}\ln\frac{a_i}{n},
\qquad
H(V) = -\sum_j \frac{b_j}{n}\ln\frac{b_j}{n},
\qquad
I(U,V) = \sum_{i,j} \frac{n_{ij}}{n}\ln\frac{n\,n_{ij}}{a_i\,b_j}.
$$

**NMI** normalizes by the arithmetic mean of the two entropies, the variant
Vinh et al. (2010) call $\mathrm{NMI}_{\mathrm{sum}}$
([JMLR 11:2837–2854](https://jmlr.org/papers/v11/vinh10a.html)). It is not
corrected for chance.

$$
\mathrm{NMI}(U,V) = \frac{I(U,V)}{\tfrac12\bigl(H(U)+H(V)\bigr)}
$$

**AMI** is the arithmetic-mean variant of the adjusted mutual information of
Vinh et al. (2010), their $\mathrm{AMI}_{\mathrm{sum}}$. The expected mutual
information $\mathbb{E}[I]$ is taken under the permutation model (random
labelings with the same community sizes $a_i$ and $b_j$):

$$
\mathrm{AMI}(U,V) = \frac{I(U,V) - \mathbb{E}[I]}{\tfrac12\bigl(H(U)+H(V)\bigr) - \mathbb{E}[I]},
$$

$$
\mathbb{E}[I] = \sum_{i}\sum_{j}\sum_{k=\max(1,\,a_i+b_j-n)}^{\min(a_i,\,b_j)}
\frac{k}{n}\ln\frac{n\,k}{a_i\,b_j}\cdot
\frac{a_i!\,b_j!\,(n-a_i)!\,(n-b_j)!}{n!\,k!\,(a_i-k)!\,(b_j-k)!\,(n-a_i-b_j+k)!}.
$$

The expectation is computed as this full triple sum over every pair of
communities and every feasible cell count $k$, with no truncation or asymptotic
approximation of the sum. Each hypergeometric weight is accumulated in log
space and then exponentiated, because computed directly it overflows for
moderate $n$. The log-factorials
come from a Lanczos log-gamma approximation that the source
(`src/core/utils/special.rs`) states is accurate to about $2\times10^{-10}$.
The paper's sum (their eq. 2) starts at $\max(0, a_i+b_j-n)$; the $k = 0$
term is zero, so starting at 1 changes nothing. The source describes AMI as 0 for labelings
that agree only by chance and negative for agreement worse than chance.

**ARI** is the adjusted Rand index of Hubert and Arabie (1985)
([doi:10.1007/BF01908075](https://doi.org/10.1007/BF01908075)): the number of
node pairs placed together in both partitions, corrected for its expectation.

$$
\mathrm{ARI}(U,V) = \frac{\sum_{i,j}\binom{n_{ij}}{2} - E}{\tfrac12\Bigl[\sum_i\binom{a_i}{2} + \sum_j\binom{b_j}{2}\Bigr] - E},
\qquad
E = \frac{\sum_i\binom{a_i}{2}\,\sum_j\binom{b_j}{2}}{\binom{n}{2}}.
$$

**F1** is a pair-counting score. With $T = \sum_{i,j}\binom{n_{ij}}{2}$ the
pairs placed together in both partitions, $P_U = \sum_i\binom{a_i}{2}$ and
$P_V = \sum_j\binom{b_j}{2}$, pair precision is $T/P_U$, pair recall is
$T/P_V$, and their harmonic mean is

$$
F_1(U,V) = \frac{2\,T}{P_U + P_V}.
$$

The source notes that this is not the community-matching F1 that some papers
report under the same name, so F1 values from those papers are not comparable
with this one. F1 is not corrected for chance.

??? info "Reproducing the four scores from the definitions"
    The block below recomputes the four scores in plain Python from the
    formulas above, using `math.lgamma` for the log-factorials, and prints the
    absolute difference from `gt_metrics` on the karate example. It omits the
    degenerate-input conventions listed further down.

    ```python exec="on" source="above" session="evaluation" result="text"
    import math
    from collections import Counter

    def from_definitions(U, V):
        nodes = [v for v in U if v in V]
        n = len(nodes)
        cells = Counter((U[v], V[v]) for v in nodes)
        a = Counter(U[v] for v in nodes)
        b = Counter(V[v] for v in nodes)
        h_u = -sum(c / n * math.log(c / n) for c in a.values())
        h_v = -sum(c / n * math.log(c / n) for c in b.values())
        mi = sum(c / n * math.log(n * c / (a[i] * b[j])) for (i, j), c in cells.items())
        lf = lambda x: math.lgamma(x + 1)
        emi = 0.0
        for ai in a.values():
            for bj in b.values():
                for k in range(max(1, ai + bj - n), min(ai, bj) + 1):
                    log_w = (lf(ai) + lf(bj) + lf(n - ai) + lf(n - bj) - lf(n)
                             - lf(k) - lf(ai - k) - lf(bj - k) - lf(n - ai - bj + k))
                    emi += k / n * math.log(n * k / (ai * bj)) * math.exp(log_w)
        pairs = lambda x: x * (x - 1) / 2
        t = sum(pairs(c) for c in cells.values())
        p_u = sum(pairs(c) for c in a.values())
        p_v = sum(pairs(c) for c in b.values())
        expected = p_u * p_v / pairs(n)
        return (mi / ((h_u + h_v) / 2),
                (mi - emi) / ((h_u + h_v) / 2 - emi),
                (t - expected) / ((p_u + p_v) / 2 - expected),
                2 * t / (p_u + p_v))

    for name, mine, theirs in zip(("NMI", "AMI", "ARI", "F1"),
                                  from_definitions(partition, gt),
                                  pymocd.gt_metrics(partition, gt)):
        print(f"{name}: |difference| = {abs(mine - theirs):.1e}")
    ```

## Chance correction in practice

Vinh et al. (2010, Sec. 4) observe that, between random clusterings, the
average NMI grows with the number of clusters and moves away from zero as the
number of points per cluster shrinks, while the average ARI stays close to
zero, and conclude that an adjusted-for-chance version of mutual information
is desirable. The block below scores random
labelings of the karate club with $k$ equally likely labels against the
two-faction ground truth, averaged over 200 labelings per $k$:

```python exec="on" source="above" session="evaluation" result="text"
import random
import statistics

rng = random.Random(0)
for k in (2, 4, 8, 16):
    scores = [pymocd.gt_metrics({v: rng.randrange(k) for v in G}, gt) for _ in range(200)]
    nmi_m, ami_m, ari_m, f1_m = (statistics.mean(col) for col in zip(*scores))
    print(f"k={k:>2}  NMI={nmi_m:.3f}  AMI={ami_m:+.3f}  ARI={ari_m:+.3f}  F1={f1_m:.3f}")
```

Read the columns against each other: a labeling that carries no information
about the ground truth should score near zero, which is what the
chance-corrected AMI and ARI are built to do. NMI and F1 have no such
correction, so their value for an uninformative labeling depends on $k$. The
source code warns that NMI and AMI are different quantities and that results
depend on which score a table reports; state the score, and for NMI the
normalization, whenever you report numbers.

## Isolated nodes and the label `-1`

Every detector returns the label `-1` for nodes of degree zero. The metrics
have no special case for `-1`: it is an ordinary label, so **all nodes labeled
`-1` in a dict form one community when that dict is scored**. A graph with
several isolated nodes is therefore scored as if the detector had put all of
them into a single extra community. The same rule applies to the second
argument, so a ground truth that marks unlabeled nodes with `-1` also groups
them together.

If that is not the comparison you want, change the dict before scoring:
give each `-1` node its own label to score isolated nodes as singletons, or
drop them, which excludes them from scoring because only shared nodes count.
The block below adds three isolated nodes to the karate club and, in a toy
ground truth, makes each of them its own community:

```python exec="on" source="above" session="evaluation" result="text"
H = G.copy()
H.add_nodes_from([34, 35, 36])  # three isolated nodes
part = pymocd.rimpso(H)
print("labels of the isolated nodes:", [part[v] for v in (34, 35, 36)])

gt_h = gt | {34: 2, 35: 3, 36: 4}

isolated = [v for v, c in part.items() if c == -1]
as_singletons = dict(part)
next_label = max(part.values()) + 1
for offset, v in enumerate(isolated):
    as_singletons[v] = next_label + offset
excluded = {v: c for v, c in part.items() if c != -1}

for name, p in [("as returned (-1 shared)", part),
                ("-1 nodes as singletons", as_singletons),
                ("-1 nodes dropped", excluded)]:
    print(f"{name:24s} AMI={pymocd.ami(p, gt_h):.4f}  scored nodes={sum(v in gt_h for v in p)}")
```

## Degenerate inputs

The implementation applies these conventions, read from
`src/core/metrics/` and `src/api/metrics.rs`:

- **No shared nodes.** `gt_metrics` returns `(0.0, 0.0, 0.0, 0.0)`. The single
  functions have no such check: `nmi`, `ami` and `f1` return `1.0` and `ari`
  returns `nan`.
- **One shared node.** $\binom{n}{2} = 0$, so the expected index $E$ of ARI
  is $0/0$ and ARI is `nan`, from `gt_metrics` and from `ari` alike.
- **NMI** is `1.0` when both entropies are zero, that is, when each argument
  puts all scored nodes in one community.
- **AMI** is `1.0` when both entropies are zero, and `0.0` when the
  denominator $\tfrac12(H(U)+H(V)) - \mathbb{E}[I]$ is below $10^{-15}$ in
  absolute value. When both arguments put every scored node in its own
  community, $\mathbb{E}[I] = I(U,V) = H(U) = H(V) = \ln n$, so numerator and
  denominator are both zero in exact arithmetic; the computed AMI is then
  `0.0` or `1.0` depending on rounding, as the last line printed by the block
  below shows for a few values of $n$, while NMI, ARI and F1 are `1.0`.
- **ARI** is `1.0` when the denominator is below $10^{-15}$ in absolute value,
  which the source describes as both partitions being all singletons or both
  being a single community.
- **F1** is `1.0` when neither partition has a pair of nodes in the same
  community ($P_U + P_V = 0$).

```python exec="on" source="above" session="evaluation" result="text"
cases = {
    "no shared nodes": ({0: 0, 1: 1}, {5: 0, 6: 1}),
    "one shared node": ({0: 0}, {0: 7}),
    "both one community": ({0: 0, 1: 0, 2: 0}, {0: 4, 1: 4, 2: 4}),
    "one community vs singletons": ({v: 0 for v in range(4)}, {v: v for v in range(4)}),
}
fmt = lambda xs: "(" + ", ".join(f"{x:.3g}" for x in xs) + ")"
for name, (p, g) in cases.items():
    single = (pymocd.nmi(p, g), pymocd.ami(p, g), pymocd.ari(p, g), pymocd.f1(p, g))
    print(f"{name:28s} gt_metrics={fmt(pymocd.gt_metrics(p, g))}  nmi/ami/ari/f1={fmt(single)}")

singletons = {n: {v: v for v in range(n)} for n in range(2, 7)}
print("identical, all singletons:  ",
      "  ".join(f"n={n}: AMI={pymocd.ami(s, s):.3g}" for n, s in singletons.items()))
```

## The karate ground-truth convention

A comment in `src/core/graph/karate.rs` notes that NetworkX's `club` attribute
puts node 8 with Mr. Hi, while much of the literature puts node 8 with the
Officer faction, and that a karate score which disagrees with a published one
often reflects only this difference. The two conventions differ in one node:

```python exec="on" source="above" session="evaluation" result="text"
print("node 8 in NetworkX:", G.nodes[8]["club"])
gt_alt = gt | {8: 1}  # node 8 with the Officer faction
print(f"NMI between the two conventions: {pymocd.nmi(gt, gt_alt):.4f}")
print(f"AMI of the RIMPSO partition, NetworkX convention:  {pymocd.ami(partition, gt):.4f}")
print(f"AMI of the RIMPSO partition, node 8 with Officers: {pymocd.ami(partition, gt_alt):.4f}")
```

State which convention you used when you report karate results.

## Scoring a whole front

The metrics take one partition at a time, so scoring a Pareto front is a loop.
The best score over a front uses the ground truth to choose, so it is an upper
bound for any rule that selects from that front without labels; report it
separately from the score of the partition the detector selects. See
[Working with Pareto fronts](fronts.md) for the front accessors.

```python exec="on" source="above" session="evaluation" result="text"
partitions, points, selected = pymocd.rimpso_fronts(G)
amis = [pymocd.ami(p, gt) for p in partitions]
print(f"front size: {len(partitions)}")
print(f"AMI of the selected member: {amis[selected]:.4f}")
print(f"best AMI in the front:      {max(amis):.4f}")
```

## Scoring without ground truth

The four functions need a reference partition. pymocd does not export a
quality function such as modularity; NetworkX provides one. It takes the
communities as a collection of node sets. Pass `weight=None` so that it scores
the unweighted graph the detectors see: `nx.karate_club_graph()` carries a
`weight` attribute on its edges, which the detectors ignore (see
[Input graphs](input-graphs.md#edge-weights)) but NetworkX would otherwise use:

```python exec="on" source="above" session="evaluation" result="text"
from collections import defaultdict

communities = defaultdict(set)
for v, c in partition.items():
    communities[c].add(v)
print(f"modularity: {nx.community.modularity(G, communities.values(), weight=None):.4f}")
```

## See also

- API reference: [`gt_metrics`][pymocd.gt_metrics], [`nmi`][pymocd.nmi],
  [`ami`][pymocd.ami], [`ari`][pymocd.ari], [`f1`][pymocd.f1].
- [Input graphs](input-graphs.md): node ids and isolated nodes on the input side.
- [Core concepts](concepts.md): Pareto fronts and selecting one partition.
- [Performance and reproducibility](performance.md): repeated runs of
  nondeterministic detectors.
- Gallery: [Comparing detectors on a planted partition](../examples/benchmarks/plot_compare_detectors.md)
  and [Recovery as mixing increases](../examples/benchmarks/plot_mixing_sweep.md).
