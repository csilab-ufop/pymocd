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
Node ids and community labels can be any hashable values, such as integers of
any size, strings or tuples. They are matched the way dict keys are, so `1`,
`1.0` and `True` name the same node or community, and ground truth stored as
names can be scored directly. Labels are identifiers only: renaming the
communities of either argument changes no score beyond floating-point
rounding. All four formulas below are symmetric in their two arguments.

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
compared on the shared nodes and is not penalized for the omission. If the two
dicts share no node at all and are not both empty, every function raises
`ValueError`: nothing could be scored, and the usual cause is node ids of
different types, such as `str` ids read from a file against the `int` ids a
detector returns. A partial mismatch goes unnoticed, so check coverage
yourself when it matters:

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

The fraction in the sum is the hypergeometric probability that a community of
size $a_i$ and a random community of size $b_j$ share $k$ nodes. The paper's
sum (their eq. 2) starts at $\max(0, a_i+b_j-n)$; the $k = 0$ term is zero, so
starting at 1 changes nothing. The source describes AMI as 0 for labelings
that agree only by chance and negative for agreement worse than chance.

`src/core/metrics/ami.rs` does not evaluate the AMI formula as written. Write
$\bar H = \tfrac12\bigl(H(U)+H(V)\bigr)$. $\bar H$, $I(U,V)$ and
$\mathbb{E}[I]$ all contain the term $\ln n$, so when both labelings are nearly
all singletons the numerator and the denominator are differences of values
close to $\ln n$, and most of their digits cancel. Under the permutation model,
the expected overlap of community $i$ with all the communities of $V$ together
is $a_i$, and likewise $b_j$ for community $j$; with that, the shared terms
cancel exactly and

$$
n\bigl(I(U,V) - \bar H\bigr) = Y = \sum_{i,j} n_{ij}\ln\frac{n_{ij}}{\sqrt{a_i\,b_j}},
\qquad
n\bigl(\mathbb{E}[I] - \bar H\bigr) = Z = \mathbb{E}[Y],
\qquad
\mathrm{AMI}(U,V) = \frac{Z - Y}{Z},
$$

where $\mathbb{E}[Y]$ is the expectation of the same sum under the permutation
model. Every term of $Y$ and $Z$ is at most zero, since
$n_{ij} \le \min(a_i, b_j) \le \sqrt{a_i b_j}$, so neither sum cancels.

The terms of $Z$ depend on a pair of communities only through their sizes, so
the code sums over pairs of *distinct* sizes and multiplies each by the number
of community pairs with those sizes. A labeling of $n$ nodes has fewer than
$\sqrt{2n}$ distinct community sizes, so this visits at most $2n$ size pairs,
where the sum as written visits every pair of communities ($n^2$ when both
sides are all singletons). For each size pair, the probabilities are obtained
from the exact ratio of consecutive terms, starting at the most likely $k$ and
normalized by their sum, rather than from log-factorials: a log-factorial
difference cancels terms of size $n \ln n$ and leaves a relative error of order
$10^{-9}$ in each weight at $n = 10^6$. Terms smaller than $10^{-40}$ times the
largest one for that size pair are skipped.

In NMI and AMI, the logarithm of a ratio of counts is computed as
$\ln(1+x)$ with $x$ formed from the exact integer difference whenever the
ratio is at least $\tfrac12$, so that a ratio close to 1, such as that of a
community holding almost every node, keeps its low bits. The floating-point
sums use compensated (Neumaier) summation.

**ARI** is the adjusted Rand index of Hubert and Arabie (1985)
([doi:10.1007/BF01908075](https://doi.org/10.1007/BF01908075)): the number of
node pairs placed together in both partitions, corrected for its expectation.
The pair counts are exact integers, and the formula below is evaluated as one
integer ratio with its numerator and denominator rounded once each, as
scikit-learn does; in floating point, $E$ rounds and the subtraction then
cancels when both partitions are nearly a single community.

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
`src/core/metrics/` and `src/api/metrics.rs`. Except where a bullet says
otherwise, scikit-learn returns the same values up to floating-point rounding.

- **Same partition.** When the two arguments are the same partition up to
  renaming the communities, all four scores are exactly `1.0`. This includes
  both arguments putting every scored node in one community, and both putting
  every scored node in its own community, where AMI is $0/0$ in exact
  arithmetic because $I(U,V) = \mathbb{E}[I] = H(U) = H(V) = \ln n$.
- **No shared nodes.** All five functions raise `ValueError` when the two
  dicts share no node and are not both empty, as described in
  [Which nodes are scored](#which-nodes-are-scored); scikit-learn takes two
  aligned arrays and has no such case. Two empty dicts are the same (empty)
  partition and score `1.0`, as scikit-learn scores empty labelings.
- **One shared node.** One-node labelings are the same partition, so all four
  scores are `1.0`, as scikit-learn returns for one-element labelings.
- **One argument is a single community** and the other is not: NMI and AMI
  are `0.0`, since $I(U,V) = 0$.
- **One argument is all singletons** and the other is not: AMI is exactly
  `0.0`, because then $\mathbb{E}[I] = I(U,V)$, both equal to the entropy of the
  other argument. scikit-learn does not special-case this; it evaluates the
  formula and returns rounding noise, up to $2.8\times10^{-6}$ in the
  comparison in the next section.
- **AMI near zero.** scikit-learn moves the numerator and the denominator of
  AMI at least machine epsilon ($2.2\times10^{-16}$) away from zero, keeping
  their signs, to avoid $0/0$. pymocd does not: outside the cases above, its
  denominator $-Z/n$ is a sum of positive terms and is never zero, and moving
  the numerator would add up to $2.2\times10^{-16}$ divided by the denominator
  to the result, $3.2\times10^{-10}$ when each argument has one pair and
  otherwise singletons on $2\times10^6$ nodes.
- **ARI** is `1.0` when no pair of scored nodes is together in one argument and
  apart in the other, which covers the same-partition cases above.
- **F1** is `1.0` when neither partition has a pair of nodes in the same
  community ($P_U + P_V = 0$), where its ratio is $0/0$. scikit-learn has no
  pairwise F1.

```python exec="on" source="above" session="evaluation" result="text"
cases = {
    "two empty dicts": ({}, {}),
    "one shared node": ({0: 0}, {0: 7}),
    "both one community": ({0: 0, 1: 0, 2: 0}, {0: 4, 1: 4, 2: 4}),
    "one community vs singletons": ({v: 0 for v in range(4)}, {v: v for v in range(4)}),
    "singletons vs a pair": ({v: v for v in range(4)}, {0: 0, 1: 0, 2: 1, 3: 2}),
}
fmt = lambda xs: "(" + ", ".join(f"{x:.3g}" for x in xs) + ")"
for name, (p, g) in cases.items():
    single = (pymocd.nmi(p, g), pymocd.ami(p, g), pymocd.ari(p, g), pymocd.f1(p, g))
    print(f"{name:28s} gt_metrics={fmt(pymocd.gt_metrics(p, g))}  nmi/ami/ari/f1={fmt(single)}")

singletons = {n: {v: v for v in range(n)} for n in (2, 4, 14, 309)}
print("identical, all singletons:  ",
      "  ".join(f"n={n}: AMI={pymocd.ami(s, s):.3g}" for n, s in singletons.items()))

try:
    pymocd.gt_metrics({0: 0, 1: 0, 2: 1}, {"0": 0, "1": 0, "2": 1})
except ValueError as err:
    print("int ids against str ids:     ValueError:", err)
```

## Agreement with scikit-learn

NMI, AMI and ARI follow scikit-learn's `normalized_mutual_info_score`,
`adjusted_mutual_info_score` and `adjusted_rand_score` with their default
`average_method="arithmetic"`, including the conventions above (with the
exception of AMI near zero); the comparisons below used scikit-learn 1.7.2.
`tests/test_metrics.py` compares the two on random labelings and on the
degenerate cases when scikit-learn is installed. On 30,000 random pairs of
labelings with 1 to $10^6$ nodes:

- **ARI** returned the same float as scikit-learn in every pair.
- **NMI** differed from scikit-learn by at most $8.8\times10^{-14}$.
- **AMI** agreed within $10^{-12}$ in 26,498 of the 27,962 pairs where it was
  computed. The other 1,464 were compared with the exact value, computed with
  40-digit arithmetic, in the 1,447 where that was affordable: pymocd was
  within $2.9\times10^{-16}$ of it in all of them, scikit-learn up to
  $2.8\times10^{-6}$ away. The remaining 17 differ from scikit-learn by at
  most $4.3\times10^{-11}$.

scikit-learn's AMI error is largest when one argument is all or almost all
singletons, which makes the denominator small, and with thousands of small
communities, where its expected mutual information adds up millions of terms
one at a time and computes each probability from log-factorials. Against the
exact value, for example, two independent random labelings of $10^5$ nodes
into $10^5$ labels score $-1.2097898\times10^{-5}$: pymocd is $3\times10^{-17}$
away and scikit-learn $4.2\times10^{-7}$. The worst conditioned inputs are
labelings that are nearly all singletons on both sides. On 50 of them, with 1
to 100 small communities per side among $10^3$ to $2\times10^6$ nodes, pymocd
was within $2.9\times10^{-16}$ of the 60-digit exact value; scikit-learn, run
up to $10^4$ nodes, was up to $1.0\times10^{-4}$ away. On twelve labelings of
$10^4$ to $10^6$ nodes timed on one machine, AMI took under 0.25 s each and
scikit-learn 2 to 34 s.

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
- [Benchmark graphs](benchmarks.md): LFR graphs with a planted ground truth to
  score against.
- Gallery: [Comparing detectors on a planted partition](../examples/benchmarks/plot_compare_detectors.md)
  and [Recovery as mixing increases](../examples/benchmarks/plot_mixing_sweep.md).
