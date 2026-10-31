---
description: The ideas behind pymocd — partitions, modularity and the Constant Potts Model, Pareto dominance and fronts, and the rules that select one partition from a front.
---

# Core concepts

pymocd treats community detection as an optimization problem over partitions of
a graph. Most of its detectors optimize two or more objectives at once, so their
search does not end with one answer: it ends with a set of trade-off partitions,
a Pareto front, and a selection rule then picks one member to return (two
re-implemented baselines, GDPSO and CDRME, optimize a single objective). This
page defines the pieces of that pipeline: what a partition is and how pymocd returns it, the
two classical quality functions the library's own methods are built on
(modularity and the Constant Potts Model), Pareto dominance and fronts, the
selection rules the detectors apply, and one caveat for evaluation: the
ground-truth partition need not be on the front at all.

## Partitions and what pymocd returns

Let $G = (V, E)$ be an undirected graph with $n = \lvert V \rvert$ nodes and
$m = \lvert E \rvert$ edges. A **partition** $\mathcal{C}$ divides $V$ into
disjoint communities: every node belongs to exactly one community
$c \in \mathcal{C}$.

Every detector returns a partition as a `dict[int, int]` that maps each node id
to a community id. The ids are labels: only which nodes share an id carries
information, so partitions are compared with label-invariant measures (see
[Evaluating partitions](evaluation.md)) rather than by their ids. Nodes without
edges (isolated nodes) are given the community id `-1` by every detector. The
detectors read the graph's nodes and edges and nothing else, so edge attributes
such as `weight` play no role; [Input graphs](input-graphs.md) covers the
accepted graph types and node ids.

```python exec="on" source="above" session="concepts" result="text"
from collections import Counter

import networkx as nx
import pymocd

G = nx.karate_club_graph()          # Zachary's karate club
partition = pymocd.rimpso(G)        # dict[node, community]

print("nodes:", G.number_of_nodes(), "| edges:", G.number_of_edges())
print("first five entries:", dict(list(partition.items())[:5]))
sizes = Counter(partition.values())
print("communities:", len(sizes), "| sizes:", sorted(sizes.values(), reverse=True))

# NetworkX's community functions expect a list of node sets instead.
groups = {}
for node, c in partition.items():
    groups.setdefault(c, set()).add(node)
print("as node sets:", [sorted(g) for g in groups.values()][:2], "...")
```

## Community detection as optimization

A quality function scores a partition, and a detector searches for a partition
with a good score. For a community $c$, write

- $n_c$ for its number of nodes,
- $e_c$ for the number of edges with both endpoints in $c$, each counted once
  (the HP-MOCD and Shi-MOCD module READMEs write it $l_c$),
- $d_c = \sum_{v \in c} k_v$ for its degree sum, where $k_v$ is the degree of
  node $v$; every internal edge contributes twice to $d_c$.

### Modularity

The modularity of [Newman and Girvan (2004)](https://doi.org/10.1103/PhysRevE.69.026113)
is

$$
Q(\mathcal{C}) = \sum_{c \in \mathcal{C}} \left[ \frac{e_c}{m} - \left( \frac{d_c}{2m} \right)^{2} \right].
$$

The first term is the fraction of edges that fall inside communities. The
second is a null term: $m\,(d_c / 2m)^2 = d_c^2 / 4m$ is the number of internal
edges that the degree sequence alone would produce (the RIMPSO module README
calls it $E_c$). A high $Q$ means that communities hold more internal edges
than the degrees explain. Maximizing $Q$ has a known blind spot:
[Fortunato and Barthélemy (2007)](https://doi.org/10.1073/pnas.0605965104) show
that it may fail to identify modules smaller than a scale that depends on the
size of the network and on how interconnected the modules are, even when the
modules are unambiguously defined (the resolution limit).

### The Constant Potts Model

The Constant Potts Model (CPM) of
[Traag, Van Dooren and Nesterov (2011)](https://doi.org/10.1103/PhysRevE.84.016114)
replaces the degree-based null term with a constant density $\gamma$, the
resolution. Written, as in the RIMPSO docstring, as a quantity to maximize:

$$
H_\gamma(\mathcal{C}) = \sum_{c \in \mathcal{C}} \left[ e_c - \gamma \binom{n_c}{2} \right].
$$

Each community is rewarded for its internal edges and charged $\gamma$ for
every pair of its nodes, so a community pays off only when it is internally
denser than $\gamma$. The RIMPSO module README draws the consequence for the
range of useful resolutions: above $\gamma = 1$ everything shatters, and on a
connected graph nothing splits below $1/n^2$, with $n$ the number of nodes
(separate components cut no edge between them, so they split apart at any
$\gamma > 0$). The catch is that $\gamma$ is a parameter, and different values
give different partitions.

The block below evaluates both functions directly from these definitions on the
two factions of the karate club (the `club` node attribute in NetworkX). NetworkX's
karate graph carries edge weights; pymocd ignores them, so the comparison with
NetworkX passes `weight=None`.

```python exec="on" source="above" session="concepts" result="text"
def stats(G, partition):
    """Per community: node count n_c, internal edges e_c, degree sum d_c."""
    n_c, e_c, d_c = Counter(), Counter(), Counter()
    for v, c in partition.items():
        n_c[c] += 1
        d_c[c] += G.degree(v)                 # unweighted degree
    for u, v in G.edges():
        if partition[u] == partition[v]:
            e_c[partition[u]] += 1
    return n_c, e_c, d_c

def modularity(G, partition):
    m = G.number_of_edges()
    n_c, e_c, d_c = stats(G, partition)
    return sum(e_c[c] / m - (d_c[c] / (2 * m)) ** 2 for c in n_c)

def cpm(G, partition, gamma):
    n_c, e_c, _ = stats(G, partition)
    return sum(e_c[c] - gamma * n_c[c] * (n_c[c] - 1) / 2 for c in n_c)

club = {v: int(G.nodes[v]["club"] != "Mr. Hi") for v in G}   # the two factions
factions = [{v for v in G if club[v] == k} for k in (0, 1)]

print("Q from the formula:", round(modularity(G, club), 6))
print("Q from NetworkX:   ", round(nx.community.modularity(G, factions, weight=None), 6))
for gamma in (0.05, 0.1, 0.5):
    print(f"H at gamma = {gamma}: {cpm(G, club, gamma):.2f}")
```

## Multi-objective formulation

Instead of one score, a multi-objective detector evaluates a partition with a
vector of $M$ objectives,
$f(\mathcal{C}) = (f_1(\mathcal{C}), \dots, f_M(\mathcal{C}))$. Here every
objective is written as minimized; an objective that a paper maximizes enters
with its sign flipped, which leaves the definitions below unchanged.

**Pareto dominance.** A partition $\mathcal{A}$ dominates a partition
$\mathcal{B}$, written $\mathcal{A} \prec \mathcal{B}$, when it is no worse on
every objective and strictly better on at least one:

$$
\mathcal{A} \prec \mathcal{B}
\iff
f_i(\mathcal{A}) \le f_i(\mathcal{B}) \;\; \text{for all } i
\quad \text{and} \quad
f_j(\mathcal{A}) < f_j(\mathcal{B}) \;\; \text{for some } j .
$$

This is the test the library's code applies (for example
`pareto/dominance.rs` in the RIMPSO module).

**Non-dominated set and Pareto front.** Among a set of partitions, the
non-dominated ones are those that no other member dominates. Over all possible
partitions, the non-dominated set is the Pareto-optimal set, and its image in
objective space is the Pareto front. No member of a front is better than
another on every objective: each is a different trade-off. The searches in
pymocd are heuristics, so what a detector holds at the end is a set of mutually
non-dominated partitions that approximates the Pareto front: the rank-1 members of the final population for the
NSGA-II and NSGA-III detectors (for MMCoMO, of its two co-evolving populations
merged), or the external archive for RIMPSO's swarm and Shi-MOCD's PESA-II. The documentation calls that set the detector's *front*,
and the `*_fronts` functions return it as a list of partitions.

Two objective pairs matter most in pymocd, because the library's own methods
are built on them. The other detectors' objectives are listed on the
[Algorithms](../algorithms/index.md) pages.

### Decomposed modularity: HP-MOCD and Shi-MOCD

[HP-MOCD](../algorithms/hpmocd.md) and [Shi-MOCD](../algorithms/shi-mocd.md)
minimize the pair of [Shi et al. (2012)](https://doi.org/10.1016/j.asoc.2011.10.005),
Eqs. 3.5 and 3.6:

$$
\operatorname{intra}(\mathcal{C}) = 1 - \frac{1}{m} \sum_{c} e_c,
\qquad
\operatorname{inter}(\mathcal{C}) = \sum_{c} \left( \frac{d_c}{2m} \right)^{2}.
$$

Despite the names, $\operatorname{intra}$ is the fraction of edges that run
*between* communities, and $\operatorname{inter}$ is modularity's null term. The
pair is an exact decomposition of modularity:

$$
Q(\mathcal{C}) = 1 - \operatorname{intra}(\mathcal{C}) - \operatorname{inter}(\mathcal{C}).
$$

The two objectives conflict. A single community scores
$(\operatorname{intra}, \operatorname{inter}) = (0, 1)$; the Shi-MOCD module
README summarizes the trade-off as "fragmenting lowers inter and raises intra;
coarsening does the reverse", which is why the front spans the range from
coarse to fine partitions.

### The CPM split: RIMPSO

[RIMPSO](../algorithms/rimpso.md) splits the Constant Potts Model the same way.
It minimizes

$$
\operatorname{cut}(\mathcal{C}) = 1 - \frac{1}{m} \sum_{c} e_c,
\qquad
\operatorname{pair}(\mathcal{C}) = \binom{n}{2}^{-1} \sum_{c} \binom{n_c}{2},
$$

where $n$ counts every node of the graph, isolated ones included, and each
isolated node is a singleton that adds nothing to the sum. $\operatorname{cut}$
is the fraction of edges leaving their community (the same quantity as
$\operatorname{intra}$), and $\operatorname{pair}$ is the fraction of node pairs
that share a community. A single community scores
$(\operatorname{cut}, \operatorname{pair}) = (0, 1)$ and all singletons score
$(1, 0)$. Where $\operatorname{inter}$ weighs a community by its degree,
$\operatorname{pair}$ weighs it by its node count: it is the expected fraction
of edges inside communities if the $m$ edges were placed uniformly at random
among node pairs.

Let $\gamma_d = m / \binom{n}{2}$ be the edge density of the graph. For every
partition and every resolution,

$$
H_\gamma(\mathcal{C}) = m \left[ 1 - \operatorname{cut}(\mathcal{C}) - \frac{\gamma}{\gamma_d}\, \operatorname{pair}(\mathcal{C}) \right].
$$

Maximizing CPM at resolution $\gamma$ is therefore the same as minimizing the
weighted sum $\operatorname{cut} + (\gamma / \gamma_d)\,\operatorname{pair}$,
and minimizing any positive weighting $a\,\operatorname{cut} + b\,\operatorname{pair}$
is the same as maximizing CPM at $\gamma = \gamma_d\, b / a$. The resolution is
not a search parameter but the exchange rate between the two objectives, and a
CPM optimum at any $\gamma > 0$ lies on the
$(\operatorname{cut}, \operatorname{pair})$ Pareto front. That is the sense in
which the RIMPSO docstring and module README call the front the graph's
**resolution profile**: one front, computed once, answers every resolution. The
unit test `the_split_reproduces_cpm_at_every_resolution` in
`objectives/cpm.rs` checks the identity.

The block below computes $(\operatorname{cut}, \operatorname{pair})$ by hand for
every member of RIMPSO's front on the karate club and compares the result with
the points that [`rimpso_fronts`][pymocd.rimpso_fronts] reports. It then reads
the front at several resolutions: for each $\gamma$, the member that minimizes
$\operatorname{cut} + (\gamma/\gamma_d)\,\operatorname{pair}$ is the best
partition the front holds for CPM at that $\gamma$.

```python exec="on" source="above" session="concepts" result="text"
def cut_pair(G, partition):
    n, m = G.number_of_nodes(), G.number_of_edges()
    n_c, e_c, _ = stats(G, partition)
    cut = 1 - sum(e_c.values()) / m
    # Label -1 marks isolated nodes; each one is a singleton with no pairs.
    pairs = sum(s * (s - 1) / 2 for c, s in n_c.items() if c != -1)
    return cut, pairs / (n * (n - 1) / 2)

partitions, points, selected = pymocd.rimpso_fronts(G)
gap = max(
    max(abs(a - b) for a, b in zip(cut_pair(G, p), point))
    for p, point in zip(partitions, points)
)
print("front members:", len(partitions))
print("largest difference, hand-computed vs reported (cut, pair):", gap)

n, m = G.number_of_nodes(), G.number_of_edges()
gamma_d = m / (n * (n - 1) / 2)
print(f"edge density gamma_d = {gamma_d:.4f}")
for gamma in (0.01, 0.05, 0.1, 0.2, 0.5, 1.0):
    i = min(range(len(points)), key=lambda i: points[i][0] + gamma / gamma_d * points[i][1])
    k = len(set(partitions[i].values()))
    print(f"gamma = {gamma:<4}  best member {i:>2}: {k:>2} communities, H = {cpm(G, partitions[i], gamma):7.2f}")
```

As $\gamma$ grows, the weight on $\operatorname{pair}$ grows, so the minimizer's
$\operatorname{pair}$ can only stay the same or decrease: the best member moves
toward the fine end of the front. The last column is $H_\gamma$ of the chosen
member, evaluated directly from the definition of CPM.

## Selecting one partition from a front

A front is a set of candidates, and a detector must return one of them. Each
detector applies one fixed rule, stated in its docstring and module README:

| Rule | Uses ground truth | Applied by |
|---|---|---|
| Maximum modularity $Q$ over the front | no | [`hpmocd`][pymocd.hpmocd], [`HpMocd.run`][pymocd.HpMocd.run] (built-in objectives), [`mocd_q`][pymocd.mocd_q], [`moga_net`][pymocd.moga_net], [`ccm`][pymocd.ccm], [`krm`][pymocd.krm], [`mmcomo`][pymocd.mmcomo] |
| Best fit of a degree-corrected assortative block model, Schwarz-penalized | no | [`rimpso`][pymocd.rimpso]; the index `selected` of [`rimpso_fronts`][pymocd.rimpso_fronts]; [`rimpso_select`][pymocd.rimpso_select] for any candidate set |
| Max–min distance to random-graph control fronts | no | [`mocd_d`][pymocd.mocd_d] |
| Best match to a reference partition (the front oracle) | **yes** | none; evaluation only, applied by hand to a `*_fronts` output |

Two detectors have no front to select from. [`gdpso`][pymocd.gdpso] optimizes
modularity as a single objective and returns the best position its swarm held;
[`cdrme`][pymocd.cdrme] maximizes one scalar, the linkage sum of its paper's
Eq. (12), over a population of merge chains and selects its final partition by
maximum modularity.

### Maximum modularity

The most common rule returns the front member of highest $Q$. With the
decomposed-modularity pair this is

$$
\mathcal{C}^\star = \operatorname*{arg\,max}_{\mathcal{C} \in F} Q(\mathcal{C})
= \operatorname*{arg\,min}_{\mathcal{C} \in F} \bigl( \operatorname{intra}(\mathcal{C}) + \operatorname{inter}(\mathcal{C}) \bigr),
$$

a fixed, equal-weight scalarization of the two objectives. It needs no labels
and no parameter, but it inherits the resolution limit of modularity whatever
the front contains; the RIMPSO module README gives this as the reason maximum
modularity loses as a selector on its front. A minimizer of a positive weighted sum
also always lies on the lower convex hull of the front, so this rule can never
return a front member that lies off that hull. Applied to RIMPSO's front, the
equal-weight rule would maximize $(1 - \operatorname{cut}) - \operatorname{pair}$,
which the RIMPSO module README calls Erdős–Rényi modularity and lists among the
selection criteria it measured and rejected.

### RIMPSO: a degree-corrected assortative block model

RIMPSO's front spans every resolution, so the rule that picks from it must
decide which resolution suits the graph, without labels. RIMPSO scores each
member $\mathcal{C}$ with a degree-corrected block model in which every
community $c$ has its own internal edge density and all pairs between
communities share one density. Profiling out the densities leaves a
log-likelihood ratio against the configuration model:

$$
L(\mathcal{C}) = \sum_c e_c \ln\frac{e_c}{E_c} + e_{\text{out}} \ln\frac{e_{\text{out}}}{E_{\text{out}}},
\qquad
E_c = \frac{d_c^2}{4m}, \quad
E_{\text{out}} = m - \sum_c E_c, \quad
e_{\text{out}} = m - \sum_c e_c .
$$

Each community brings one free density, so the $B + 1$ densities are charged a
Schwarz-type penalty of half a nat per parameter per log observation, over the
$2m$ edge endpoints
([Schwarz 1978](https://doi.org/10.1214/aos/1176344136)):

$$
\operatorname{score}(\mathcal{C}) = L(\mathcal{C}) - \frac{B + 1}{2} \ln(2m),
$$

where $B$ counts the communities that hold at least one non-isolated node. The
member with the highest score is selected, ties going to the lower index. A
partition whose communities are not denser inside than between, that is with
$e_{\text{out}} > 0$ and
$e_{\text{in}} / E_{\text{in}} \le e_{\text{out}} / E_{\text{out}}$ for
$e_{\text{in}} = \sum_c e_c$ and $E_{\text{in}} = \sum_c E_c$, is rejected
outright. Neither degenerate partition needs a special case: a single community
has $L = 0$ and still pays the penalty, and all singletons put no edge inside a
community and are rejected. The [RIMPSO page](../algorithms/rimpso.md#selecting-a-partition)
gives the model behind $L$.

The block below implements the score from these equations and checks which
member it picks against the index that `rimpso_fronts` reports.

```python exec="on" source="above" session="concepts" result="text"
import math

def rimpso_score(G, partition):
    m = G.number_of_edges()
    n_c, e_c, d_c = stats(G, partition)
    live = [c for c in n_c if d_c[c] > 0]        # communities with a non-isolated node
    E = {c: d_c[c] ** 2 / (4 * m) for c in live}
    e_in, E_in = sum(e_c[c] for c in live), sum(E.values())
    e_out, E_out = m - e_in, m - E_in
    if e_out > 0 and e_in * E_out <= e_out * E_in:
        return -math.inf                          # not assortative: rejected
    L = sum(e_c[c] * math.log(e_c[c] / E[c]) for c in live if e_c[c] > 0)
    if e_out > 0 and E_out > 0:
        L += e_out * math.log(e_out / E_out)
    return L - (len(live) + 1) / 2 * math.log(2 * m)

scores = [rimpso_score(G, p) for p in partitions]
best = max(range(len(scores)), key=lambda i: scores[i])   # first maximum: lowest index
print("highest score in Python: member", best)
print("selected by rimpso_fronts: member", selected)
print("communities in the selected member:", len(set(partitions[selected].values())))
```

[`rimpso_select`][pymocd.rimpso_select] applies the same rule to candidate
partitions produced elsewhere, which separates the contribution of the search
from that of the selector; [Working with Pareto fronts](fronts.md#with-rimpsos-selection-rule)
shows it on another detector's front.

### Shi's max–min distance (MOCD-D)

Shi et al. (2012) propose a second rule, implemented by
[`mocd_d`][pymocd.mocd_d] (Shi Eqs. 3.9–3.11). It evolves full fronts
$R_1, \dots, R_r$ on $r$ = `rand_networks` Erdős–Rényi $G(n, m)$ control graphs
with the same node set and number of edges, and returns the member of the real
front $F$ whose $(\operatorname{intra}, \operatorname{inter})$ point is furthest
from all of them:

$$
D(\mathcal{C}) = \min_{1 \le j \le r} \; \min_{\mathcal{C}' \in R_j} \bigl\lVert f(\mathcal{C}) - f(\mathcal{C}') \bigr\rVert_2,
\qquad
\mathcal{C}^\star = \operatorname*{arg\,max}_{\mathcal{C} \in F} D(\mathcal{C}).
$$

The module README states the idea: a random graph's front shows what "no
community structure" looks like in objective space, so the most deviant real
point is taken as the most meaningful one. Each control front costs a full
search, so `mocd_d` runs `rand_networks + 1` searches, and the control graphs are
drawn afresh on every call. The rule is available only inside `mocd_d`; there
is no function that applies it to another front.

### Ground truth: the front oracle (evaluation only)

When a reference partition is known, as for benchmark graphs, one can score
every front member against it and keep the best one. The RIMPSO module README,
which scores members by AMI, calls that member the **front oracle**; the
`hpmocd_fronts` docstring calls the comparison *best-in-front*, or
selector-free. The oracle uses the labels, so it is not a detector: it measures
what the search found, and the gap between the oracle and a detector's own
choice measures the selection rule. Since the oracle is the maximum over the
front, no rule applied to the same front can score higher against the same
reference on the same measure. The block below scores by ARI. Several of the re-implemented papers
report best-in-front numbers; [Working with Pareto fronts](fronts.md#reproducing-best-in-front-protocols)
shows how to reproduce them.

```python exec="on" source="above" session="concepts" result="text"
ari = [pymocd.ari(p, club) for p in partitions]
oracle = max(range(len(ari)), key=lambda i: ari[i])
max_q = max(range(len(partitions)), key=lambda i: modularity(G, partitions[i]))

for name, i in [("RIMPSO rule", selected),
                ("maximum modularity", max_q),
                ("front oracle (uses labels)", oracle)]:
    k = len(set(partitions[i].values()))
    print(f"{name:<27} member {i:>2}: {k:>2} communities, ARI = {ari[i]:.3f}")
```

The oracle can also be out of reach of a label-free rule. The RIMPSO module
README reports that on the karate club the best-matching member of its
front does not compress the graph, so no description-length criterion of the
family it tested can return it.

## The ground truth may not be on the front

The front is defined by the search objectives, not by agreement with any
reference. A reference partition can be *dominated* under those objectives: some
other partition is at least as good on every objective and better on one. A
dominated partition is not Pareto-optimal, and a detector's front cannot hold
it together with a partition that dominates it, so whenever the front holds
such a partition the reference is absent, and even the front oracle then falls
short of a perfect match. The front bounds what any selection rule can recover.
A reference partition that is not dominated can still be missing, because the
search is a heuristic that may not find it.

The check is direct for RIMPSO, whose objectives can be computed for any
partition: compute the reference partition's
$(\operatorname{cut}, \operatorname{pair})$ point and look for front members that
dominate it.

```python exec="on" source="above" session="concepts" result="text"
def dominates(a, b):
    return all(x <= y for x, y in zip(a, b)) and any(x < y for x, y in zip(a, b))

def as_sets(p):
    return frozenset(frozenset(v for v in p if p[v] == c) for c in set(p.values()))

gt_point = cut_pair(G, club)
print(f"ground truth: (cut, pair) = ({gt_point[0]:.4f}, {gt_point[1]:.4f})")
print("ground truth is a front member:", as_sets(club) in {as_sets(p) for p in partitions})
for i, point in enumerate(points):
    if dominates(point, gt_point):
        k = len(set(partitions[i].values()))
        print(f"member {i:>2} dominates it: (cut, pair) = ({point[0]:.4f}, {point[1]:.4f}), {k} communities")
```

The block prints every front member that dominates the reference point; if it
prints none, the reference is non-dominated with respect to this front. For the
other detectors, the same test needs the reference's objective values under
that detector's objectives. For HP-MOCD these are $\operatorname{intra}$ and
$\operatorname{inter}$ above, and
[`HpMocd.generate_pareto_front`][pymocd.HpMocd.generate_pareto_front] returns
them for every front member.

## See also

- [Working with Pareto fronts](fronts.md): the `*_fronts` functions, their
  return shapes, and selection by hand.
- [Choosing a detector](choosing.md) and the [Algorithms](../algorithms/index.md)
  overview: every detector's objectives and selection rule.
- [Evaluating partitions](evaluation.md): the metrics behind the front oracle.
- Gallery: [the resolution profile](../examples/fronts/plot_resolution_profile.md)
  and [choosing from a front](../examples/fronts/plot_choosing_from_front.md).
- API: [`rimpso_fronts`][pymocd.rimpso_fronts],
  [`rimpso_select`][pymocd.rimpso_select],
  [`hpmocd_fronts`][pymocd.hpmocd_fronts],
  [`HpMocd`][pymocd.HpMocd].
