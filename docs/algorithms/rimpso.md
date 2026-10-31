---
description: RIMPSO, pymocd's own detector, is a memetic particle swarm over an exact two-objective split of the Constant Potts Model whose Pareto front is the graph's resolution profile, with a parameter-free block-model selector.
---

# RIMPSO

<span class="pm-pill pm-pill--own">pymocd method</span>
<span class="pm-pill pm-pill--neutral">Multi-objective</span>
<span class="pm-pill pm-pill--neutral">Deterministic</span>

RIMPSO (Resolution-Indexed Memetic Particle Swarm Optimisation) minimizes two
objectives: the fraction of edges cut, and the fraction of node pairs placed in
the same community. Together they form an exact split of the Constant Potts
Model (CPM), so the Pareto front over them holds the CPM optimum at every
resolution. A particle swarm approximates that front. Each particle is pinned to one
rung of a geometric resolution ladder, and the swarm periodically repairs
particles and merges communities. A selector then picks one member of the front
without labels or parameters: the partition that a degree-corrected assortative
block model fits best once its free densities are penalized.
[`rimpso`][pymocd.rimpso] returns that partition,
[`rimpso_fronts`][pymocd.rimpso_fronts] returns the whole front, and
[`rimpso_select`][pymocd.rimpso_select] applies the selector to candidate
partitions produced elsewhere.

<div class="pm-facts" markdown>

| | |
|---|---|
| **Entry points** | [`rimpso`][pymocd.rimpso], [`rimpso_fronts`][pymocd.rimpso_fronts], [`rimpso_select`][pymocd.rimpso_select] |
| **Objectives** | `cut`, the fraction of edges between communities, and `pair`, the fraction of node pairs within a community. Both are minimized, and together they are an exact split of CPM. |
| **Search** | Memetic particle swarm with one particle per rung of a geometric resolution ladder over $[1/n^2, 1]$. Every `ls_period` iterations it runs a local repair and a community-merge sweep. The archive is truncated by rung. |
| **Selection rule** | The front member with the highest Schwarz-penalized fit under a degree-corrected assortative block model. It uses no labels and has no parameter. |
| **Front accessor** | [`rimpso_fronts`][pymocd.rimpso_fronts], plus [`rimpso_select`][pymocd.rimpso_select] for external candidates |
| **Deterministic** | Yes, for a fixed input. The same graph (same nodes and edges, listed in the same order) and the same parameters, `seed` included, give the same partition on any number of threads (function docstring; module README, section 8). Listing the nodes or edges in another order can change the partition (see [Input order](../guide/performance.md#input-order)). |
| **Original code** | This library ([`src/core/algorithms/rimpso`](https://github.com/csilab-ufop/pymocd/tree/master/src/core/algorithms/rimpso)) |
| **Paper** | Santos, in preparation (2026). It has no DOI yet. |

</div>

## Objectives

Let the graph have $n$ nodes and $m$ edges. Isolated nodes count toward $n$. For
a community $c$ of a partition $\mathcal{C}$, let $n_c$ be its number of nodes
and $e_c$ its number of internal edges. The Constant Potts Model of
[Traag, Van Dooren and Nesterov (2011)](https://doi.org/10.1103/PhysRevE.84.016114)
at resolution $\gamma$ is, in the sign convention the module uses (a quantity to
maximize),

$$
H_\gamma(\mathcal{C}) = \sum_{c \in \mathcal{C}} \left[ e_c - \gamma \binom{n_c}{2} \right].
$$

RIMPSO never scores $H_\gamma$ directly. It minimizes two terms that do not
contain $\gamma$:

$$
\operatorname{cut}(\mathcal{C}) = 1 - \frac{1}{m} \sum_{c} e_c,
\qquad
\operatorname{pair}(\mathcal{C}) = \binom{n}{2}^{-1} \sum_{c} \binom{n_c}{2}.
$$

`cut` is the fraction of edges that leave their community. The module README
calls it the partition's own realized mixing parameter. `pair` is the fraction
of node pairs that share a community. A single community scores
$(\operatorname{cut}, \operatorname{pair}) = (0, 1)$, and all singletons score
$(1, 0)$. The function docstring notes that this mirrors how
[HP-MOCD](hpmocd.md) splits modularity.

Let $\gamma_d = m / \binom{n}{2}$ be the edge density of the graph. Then for
every partition and every resolution

$$
H_\gamma(\mathcal{C}) = m \left[ 1 - \operatorname{cut}(\mathcal{C}) - \frac{\gamma}{\gamma_d}\, \operatorname{pair}(\mathcal{C}) \right].
$$

Maximizing CPM at resolution $\gamma$ is therefore the same as minimizing
$\operatorname{cut} + (\gamma/\gamma_d)\,\operatorname{pair}$. Minimizing any
positive weighting $a\,\operatorname{cut} + b\,\operatorname{pair}$ is the same
as maximizing CPM at $\gamma = \gamma_d\, b / a$. The resolution is not a search
parameter. It is the exchange rate between the two objectives. A CPM optimum at
any $\gamma > 0$ minimizes such a weighting, so it lies on the
$(\operatorname{cut}, \operatorname{pair})$ Pareto front, and one front holds
the optimum at every resolution. The module README calls this front the graph's
**resolution profile**. The archive that RIMPSO returns is a heuristic
approximation of it: at each rung it holds the best partition the search found,
not a certified optimum. The unit test
`the_split_reproduces_cpm_at_every_resolution` in `objectives/cpm.rs` checks the
identity.

Both raw counts, the number of intra-community edges and
$\sum_c \binom{n_c}{2}$, are integers. They are updated incrementally on every
node move, so they cannot drift from a full recount. On a graph with no edges
`cut` is defined as 0, and with fewer than two nodes `pair` is defined as 0.

## Search

**Resolution ladder.** The swarm has $P$ particles, where $P$ is `pop_size`.
Particle $k \in \{0, \dots, P-1\}$ is pinned to the resolution

$$
\gamma_k = n^{\,2k/(P-1) - 2},
$$

a geometric ladder from $\gamma_0 = 1/n^2$ to $\gamma_{P-1} = 1$. The module
README gives the reasons for this range. A CPM community must be internally
denser than $\gamma$ to survive, so above 1 every community shatters and, on a
connected graph, nothing splits below $1/n^2$ (separate components cut no edge
between them, so they split apart at any $\gamma > 0$). The spacing is
geometric because resolution is a multiplicative quantity. If $n < 2$ or the graph has no edges, every rung is 1.

**Particles.** A particle holds a partition (one label per node), a velocity
$v_j \in [0, 1]$ per node, and a personal best: the partition with the highest
$H_{\gamma_k}$ the particle has seen. The velocity $v_j$ is the probability
that node $j$ is unstable.

**Seeding.** Every non-isolated node takes the label of a random neighbor. The start is
locally coherent but unoptimized. RIMPSO has no seeding local search and no
`seed_rounds` parameter.

**Flight.** Each iteration updates every particle in parallel. A leader is drawn
from the archive by a binary tournament on crowding distance
([Deb et al. 2002](https://doi.org/10.1109/4235.996017)): of two uniformly drawn
members, the one with the larger crowding distance wins, so leaders tend to come
from sparsely populated parts of the front. Two pulls are drawn once per particle per
iteration: $r_1 = U(0,1) \cdot$ `cognitive` and $r_2 = U(0,1) \cdot$ `social`.
For each non-isolated node $j$, with current label $x_j$, personal-best label
$b_j$ and leader label $\ell_j$:

$$
v_j \leftarrow \operatorname{clamp}\!\left( \omega\, v_j + r_1\,[\,b_j \neq x_j\,] + r_2\,[\,\ell_j \neq x_j\,],\ 0,\ 1 \right),
$$

where $\omega$ is `inertia` and $[\cdot]$ is 1 when the condition holds and 0
otherwise. A single uniform draw $u$ then decides between two mutually
exclusive branches:

- If $u < v_j$, the node is unstable. It adopts $b_j$ or $\ell_j$ in proportion
  to the two pulls. If neither attractor pulls, it takes a random neighbor's
  label.
- If $u \ge v_j$ and $u > 1 -$ `local_rate`, the node is stable and gets one
  resolution-directed local move, the best CPM move at the particle's own
  $\gamma_k$. This branch runs only on iterations where the full local search
  does not.

**Local search.** Every `ls_period`-th iteration (0 disables it), the flight
ends with two steps:

1. **Repair.** Local moves drive the particle back to a local optimum of
   $H_{\gamma_k}$. Repair starts from a queue of only the vertices the flight
   disturbed and their neighbors, so its work scales with what the
   perturbation touched rather than with $n$.
2. **Merge sweep.** The sweep merges pairs of communities $a, b$ whose CPM gain
   $e_{ab} - \gamma_k\, n_a n_b$ is positive, where $e_{ab}$ is the number of
   edges between them. Candidate merges are taken in decreasing order of gain,
   and the sweep accepts a matching, so each community takes part in at most one
   merge per sweep. Single-node moves cannot merge two communities
   that are each stable. The repair refines, and the merge coarsens.

The module README reports ablations of both steps. Without the repair, the
flight's perturbations are never undone and its net effect on the README's LFR
grid is negative. Without the merge, the community count inflates on the
README's real-world graphs (dblp, amazon and youtube). After the flight, the particle's labels are canonicalized. Its
personal best is replaced when $H_{\gamma_k}$ has improved.

**Archive.** After seeding and after every iteration, each particle's partition
is offered to the external archive in slot order:

- **Offer.** A candidate is rejected if an archive member dominates it or
  already occupies its point. Otherwise it evicts every member it dominates and
  joins.
- **Prune.** When the archive holds more than `archive` members, it keeps, for
  each rung $k$, the member that minimizes
  $\operatorname{cut} + w_k\,\operatorname{pair}$ with $w_k = \gamma_k / \gamma_d$.
  Ties go to the lower index. By the identity above, this keeps the best
  partition the archive holds at every resolution the swarm searches.
  Crowding distance is recomputed after each prune and is used only to select
  leaders.

The prune replaced truncation by crowding distance. The module README records
that crowding evicted the archive's best member while nothing dominated it. The
archive is also not refilled up to `archive` after a prune; the README records
that refilling hurt the selector. The front returned after the last iteration is
the archive.

## Selecting a partition

The front covers every resolution on the ladder, so the profile alone cannot
say which $\gamma$ suits the graph. RIMPSO answers with one
criterion that uses no ground truth, no target number of communities and no
parameter. In a degree-corrected assortative block model, every community $c$
has its own internal density $\omega_c$, and all pairs between communities share
one density $\omega_{\text{out}}$:

$$
\begin{aligned}
A_{ij} &\sim \operatorname{Poisson}(k_i k_j\, \omega_c) && \text{if } i, j \in c,\\
A_{ij} &\sim \operatorname{Poisson}(k_i k_j\, \omega_{\text{out}}) && \text{otherwise},
\end{aligned}
$$

where $k_i$ is the degree of node $i$. Let $d_c$ be the degree sum of community
$c$. Profiling out the densities leaves a log-likelihood ratio against the
configuration model:

$$
L(\mathcal{C}) = \sum_c e_c \ln\frac{e_c}{E_c} + e_{\text{out}} \ln\frac{e_{\text{out}}}{E_{\text{out}}},
$$

$$
E_c = \frac{d_c^2}{4m}, \qquad
E_{\text{out}} = m - \sum_c E_c, \qquad
e_{\text{out}} = m - \sum_c e_c .
$$

$E_c$ is the number of intra-community edges that the degrees alone would
produce. Each community adds a free density, so the $B + 1$ densities pay the
Schwarz penalty ([Schwarz 1978](https://doi.org/10.1214/aos/1176344136)) over
the $2m$ edge endpoints:

$$
\operatorname{score}(\mathcal{C}) = L(\mathcal{C}) - \frac{B + 1}{2} \ln(2m),
$$

where $B$ counts the communities that hold at least one non-isolated node.
RIMPSO returns the member with the highest score, and ties go to the lower
index.

A candidate whose intra-community edges are not over-represented relative to
the edges between communities is rejected outright. With
$e_{\text{in}} = \sum_c e_c$ and $E_{\text{in}} = \sum_c E_c$, that is the case
when at least one edge runs between communities and
$e_{\text{in}}\, E_{\text{out}} \le e_{\text{out}}\, E_{\text{in}}$, that is,
$e_{\text{in}} / E_{\text{in}} \le e_{\text{out}} / E_{\text{out}}$. Such a
partition is not a community structure. The module README reports that this
guard changed no choice on any of the 362 archives it measured; it matters for
`rimpso_select`, which accepts candidates from anywhere. Neither degenerate
partition needs a special case. A single community
has $L = 0$ and scores $-\ln(2m)$. All singletons place no edge inside a
community, so the assortativity guard rejects them. There is no degeneracy
filter and no fallback stage.

### The whole front: `rimpso_fronts`

[`rimpso_fronts`][pymocd.rimpso_fronts] takes the same arguments as
[`rimpso`][pymocd.rimpso] and runs the same search. It returns a tuple
`(partitions, points, selected)`:

- `partitions` is a list of `dict[node, community]`.
- `points` holds the matching `(cut, pair)` tuples.
- `selected` is the index of the member that `rimpso` would return.

```python exec="on" source="above" session="rimpso" result="text"
import networkx as nx
import pymocd

G = nx.karate_club_graph()
partitions, points, selected = pymocd.rimpso_fronts(G)

print("front size:", len(partitions), "| selected index:", selected)
print("rimpso returns the selected member:", pymocd.rimpso(G) == partitions[selected])

# Walk the resolution profile from coarse (low cut) to fine (high cut).
for i in sorted(range(len(points)), key=lambda i: points[i][0])[::8]:
    cut, pair = points[i]
    k = len(set(partitions[i].values()))
    print(f"cut={cut:.3f}  pair={pair:.3f}  communities={k}")
```

The identity in [Objectives](#objectives) turns any front point into its CPM
value at any resolution, without recounting edges:

```python exec="on" source="above" session="rimpso" result="text"
from collections import Counter

n, m = G.number_of_nodes(), G.number_of_edges()
gamma_d = m / (n * (n - 1) / 2)
part = partitions[selected]
cut, pair = points[selected]

gamma = 0.3
internal = sum(part[u] == part[v] for u, v in G.edges())
penalty = sum(s * (s - 1) / 2 for s in Counter(part.values()).values())
print("H_gamma counted directly:", round(internal - gamma * penalty, 6))
print("H_gamma from (cut, pair): ", round(m * (1 - cut - gamma / gamma_d * pair), 6))
```

### The selector alone: `rimpso_select`

[`rimpso_select`][pymocd.rimpso_select] applies the same rule to a list of
`dict[node, community]` candidates produced elsewhere. One use is a front pooled
from an external CPM solver. The module README describes this as the control
that separates the search's contribution from the selector's. The function
returns a two-element list `[selected_index, points]`, where `points` holds each
candidate's `[cut, pair]`.

```python exec="on" source="above" session="rimpso" result="text"
resolutions = [0.25, 0.5, 1.0, 2.0, 4.0]
candidates = []
for r in resolutions:
    comms = nx.community.louvain_communities(G, resolution=r, seed=0)
    candidates.append({v: i for i, c in enumerate(comms) for v in c})

index, cand_points = pymocd.rimpso_select(G, candidates)
print("communities per candidate:", [len(set(c.values())) for c in candidates])
print(f"selected: candidate {index} (Louvain resolution {resolutions[index]})")
```

If every candidate is rejected, the lowest index is returned. The same happens
when the graph has no edges. The returned index does not signal either case, so
check the candidates if it matters.

!!! note "Isolated and missing nodes in candidates"
    `rimpso` and `rimpso_fronts` label every isolated node `-1`. `rimpso_select`
    treats all nodes that share a label as one community, so on a graph with
    more than one isolated node the `pair` it reports for such a candidate is
    larger than the `pair` that `rimpso_fronts` reports for the same partition.
    The selection score skips isolated nodes, so the selected index is not
    affected.

    Each candidate should label every node of the graph. A node missing from a
    candidate is given its internal position index as a label (`api.rs`), and
    that value can coincide with a community id the candidate uses.

## Parameters

[`rimpso`][pymocd.rimpso] and [`rimpso_fronts`][pymocd.rimpso_fronts] take the
same parameters with the same defaults. The `graph` argument is described in
[Input graphs](../guide/input-graphs.md).

<!-- pymocd:defaults pymocd.rimpso pymocd.rimpso_fronts -->
| Parameter | Default | Description |
|---|---|---|
| `pop_size` | `100` | Number of particles, and the number of rungs on the resolution ladder. |
| `num_gens` | `100` | Number of iterations. The search always runs all of them. |
| `inertia` | `0.4` | Fraction of a node's instability carried to the next iteration ($\omega$). |
| `cognitive` | `0.7` | Pull toward the particle's own best partition. |
| `social` | `0.7` | Pull toward a leader drawn from the archive by binary tournament on crowding distance. |
| `local_rate` | `0.35` | Per-node rate of the resolution-directed CPM local move, applied on iterations without the full local search. |
| `archive` | `100` | Archive size above which the rung prune runs; not a cap on the front (see below). |
| `ls_period` | `10` | Run the full local search (repair, then merge sweep) every `ls_period` iterations. 0 turns it off. |
| `seed` | `0` | Run seed. The default of 0 contributes nothing to the random stream. Any other value runs an independent trajectory. |

Out-of-range values that the binding accepts are clamped rather than rejected.
`pop_size` and `archive` are raised to at least 2. `inertia`, `cognitive`, `social` and `local_rate` are
clamped to $[0, 1]$, and NaN becomes 0. The integer parameters are unsigned, so
a negative integer is not clamped: the Python binding raises `OverflowError`.

`archive` sets when the prune runs, not a hard upper bound on the front. A prune
keeps one member per rung, so it leaves at most `pop_size` members. With
`archive < pop_size`, the returned front can therefore hold more than `archive`
members. At the defaults the two values are equal, so the front never holds
more than `pop_size` members.

`seed` is a random seed, not a seeding budget.

## Complexity

The module README describes the cost of individual steps:

- Node moves update the objective counts incrementally. A full recount happens
  only at seeding and after a merge sweep.
- Repair work is proportional to the vertices the perturbation disturbed, not to
  $n$.
- Selection costs one pass over the labels and one over the edges for each
  front member, and members are scored in parallel.
- The flight's `n`-sized scratch buffers are allocated once per worker thread
  for the whole run, not once per particle per iteration.

The RIMPSO paper is not published, so the bounds below are derived from this
implementation's source, not taken from the paper. Write $P$ for `pop_size`,
$T$ for `num_gens`, $L$ for `ls_period`, $A$ for `archive` and $W$ for the
number of Rayon threads.

- **Flight, per particle and iteration: $O(n + m)$.** `advance`
  (`swarm/motion.rs`) visits every node once. It reads a node's neighbors only
  when the node moves (`link_pair`) or takes the local move (`best_move`,
  `swarm/local.rs`).
- **Local search, per particle every $L$-th iteration: $O(n + m \log m)$ plus
  the repair.** `merge_sweep` (`swarm/merge.rs`) counts the edges between
  communities in $O(n + m)$, sorts at most $m$ candidate merges, and recounts
  both objectives in $O(n + m)$. The repair (`repair`, `swarm/motion.rs`) costs
  $O(\deg v)$ for each queued vertex $v$, and a vertex that moves queues itself
  and its neighbors again. The code sets no cap on the number of moves, so the
  repair's cost depends on the data.
- **Archive, per iteration, on one thread: $O(P(A + P + n))$.** Each of the $P$
  offers scans the archive, which holds at most $\max(A, P) + P$ members before
  a prune, and copies the particle's $n$ labels if it is accepted. The rung
  prune scans the archive once per rung (`pareto/archive.rs`).
- **Seeding and selection, once each: $O((P + F)(n + m))$**, for a returned
  front of $F \le \max(A, P)$ members (`swarm/init.rs`, `front/assortative.rs`).

A run therefore does $O(T P (n + m + A + P) + (T/L)\, P\, (m \log m + R))$
work, where $R$ is the repair work of one particle at one local-search
iteration. The flights, the seeding and the selection are split across the $W$
threads by particle or by front member; the archive update is not. Memory is
$O((P + A)\, n + W (n + m))$: each particle holds three length-$n$ arrays
(labels, velocities and personal best), the archive holds up to
$\max(A, P) + P$ label arrays, and each worker thread holds one set of
length-$n$ scratch arrays and a list of at most $m$ candidate merges
(`ScratchPool`, `swarm/particle.rs`).

**Measured scaling.** These runs used jackal, a shared server with an AMD Ryzen
Threadripper 7960X (48 hardware threads), with Python 3.10.12 and pymocd 6.1.0
from PyPI. The graphs came from `pymocd.lfr_graph(n, 2.5, 1.5, 0.3,
average_degree=10, max_degree=50, min_community=20, max_community=100,
seed=1)`. Each run used a fresh interpreter with the pool set by
[`max_cores`][pymocd.max_cores] to 8 threads, except where the thread count is
varied, and every figure is the median of 3 runs. CPU time is the process's
total over all threads. Peak
extra memory is the process's peak resident memory minus its resident memory
just before the call. The `num_gens=0` column times the setup alone: graph
conversion, seeding, one archive update and the selection.

| $n$ | $m$ | Default run, wall | Default run, CPU | `num_gens=0`, wall | Peak extra memory |
|---|---|---|---|---|---|
| 10,000 | 49,851 | 0.86 s | 6.43 s | 0.014 s | 18.5 MiB |
| 20,000 | 100,227 | 1.79 s | 13.38 s | 0.031 s | 38.8 MiB |
| 40,000 | 200,564 | 3.77 s | 28.18 s | 0.096 s | 82.5 MiB |
| 80,000 | 400,885 | 8.58 s | 63.43 s | 0.13 s | 163.6 MiB |

The growth is close to linear in $n$, as the bound predicts for a sparse graph:
each doubling multiplied the wall time by 2.09, 2.10 and 2.27, the CPU time by
2.08, 2.11 and 2.25, and the peak extra memory by 2.10, 2.13 and 1.98. At
$n$ = 20,000 the other terms behaved as derived:

- **`num_gens`.** CPU time was 3.39, 7.38, 13.54 and 24.88 s at 25, 50, 100 and
  200 iterations, close to linear.
- **`pop_size`.** With `archive` set equal to `pop_size`, CPU time was 3.79,
  7.82, 13.44 and 28.37 s, and peak extra memory 14.2, 23.3, 38.6 and
  70.1 MiB, at 25, 50, 100 and 200 particles.
- **`ls_period`.** CPU time was 10.52 s with the local search off (0), 14.31 s at
  the default 10 and 41.00 s at 1. An iteration that runs the repair and the
  merge sweep therefore costs about four times a flight-only iteration.
- **Threads.** At the defaults, wall time was 12.28, 6.06, 3.33, 1.81 and
  1.10 s on 1, 2, 4, 8 and 16 threads (11.2 times faster on 16), while total CPU
  time rose from 12.27 to 15.15 s. Every thread count returned a partition with
  the same 451 communities.

## When to use it

- **As the default detector.** The library recommends `rimpso` or `hpmocd` (the
  module docstring in `src/lib.rs`), and earlier versions of these docs called
  `rimpso` the recommended default. One run covers the whole resolution ladder,
  and you do not set a resolution or a number of communities.
- **When the right resolution is unknown.** The front holds the best partition
  the search found at every rung, so use `rimpso_fronts` to inspect partitions
  at several scales or to apply your own selection rule. See
  [Working with Pareto fronts](../guide/fronts.md).
- **When results must be reproducible.** The same graph (same nodes and edges,
  listed in the same order) and the same parameters, `seed` included, give the
  same partition on any number of threads; listing the nodes or edges in
  another order can change it. See
  [Performance and reproducibility](../guide/performance.md#input-order).
- **When you compare searches or selectors.** `rimpso_select` scores candidate
  sets from any source with the same rule.

Modularity maximization can fail to resolve communities smaller than a scale
that depends on the total number of edges. This is the resolution limit of
[Fortunato and Barthélemy (2007)](https://doi.org/10.1073/pnas.0605965104), who
illustrate it with a ring of cliques joined by single edges. The
module's test suite uses a ring of thirty five-node cliques as a
resolution-limit case (`beats_the_resolution_limit_on_a_ring_of_cliques` in
`api.rs`). The block below runs that graph through `rimpso` and through
NetworkX's Louvain method, a modularity-maximization heuristic. Compare both
counts with the 30 planted cliques.

```python exec="on" source="above" session="rimpso" result="text"
ring = nx.ring_of_cliques(30, 5)
print("rimpso communities:   ", len(set(pymocd.rimpso(ring).values())))
print("Louvain communities:  ", len(nx.community.louvain_communities(ring, seed=0)))
```

**Limits of the selector.** The selector cannot see ground truth, so it does not
always return the front member that best matches the reference labels (the
"front oracle", the member with the highest AMI against ground truth). The
module README reports a mean AMI gap of 0.021 to that oracle over 46 benchmark
cells. It traces the remaining gap to cases where the best-matching member does
not compress the graph, such as LFR graphs at mixing $\mu \ge 0.7$ and
Zachary's karate club, and states that no description-length criterion in the
family it tested can return such a member. If you have a selection criterion of your own,
apply it to the output of `rimpso_fronts`.

## Implementation notes

- **Determinism.** Each (iteration, particle) pair has its own random stream,
  derived from `seed`. Every tie is broken by the lower index or the lower
  label. The archive is offered candidates in slot order, never in completion
  order. As a result, the output does not depend on how Rayon schedules work
  (module README, section 8). It does depend on the input order: nodes are
  indexed in `graph.nodes()` order and each neighbor list is built in
  `graph.edges()` order (`CsrGraph::from_edges` in `src/core/graph/csr.rs`).
- **Isolated nodes** never move and are returned as community `-1`. They count
  in the `pair` denominator but add nothing to its numerator, and the selection
  score ignores them.
- **Removed parameters.** `ls_start` had no effect at its shipped value of 0 and
  was removed from the Rust API, the Python signature and the docs. There is
  no `seed_rounds` parameter, because RIMPSO has no seeding local search.
- **Deprecated aliases.** `pymocd.mr_mocd`, `pymocd.mr_mocd_fronts` and
  `pymocd.mr_mocd_select` are the names the detector had before it was renamed
  RIMPSO, and `pymocd.scale` and `pymocd.scale_fronts` are older still. All five
  are the same function objects as `rimpso`, `rimpso_fronts` and
  `rimpso_select`. They emit no warning and are not in the type stubs. New code
  should use the `rimpso*` names.

## Example

Detect communities in Zachary's karate club and check that a second call
returns the same partition:

```python exec="on" source="above" session="rimpso" result="text"
from collections import Counter

import networkx as nx
import pymocd

G = nx.karate_club_graph()
communities = pymocd.rimpso(G)  # dict[node, community]

sizes = sorted(Counter(communities.values()).values(), reverse=True)
print("communities:", len(sizes), "| sizes:", sizes)
print("same partition on a second call:", pymocd.rimpso(G) == communities)
```

## Reference

RIMPSO has not been published yet. pymocd cites it as *Santos, in preparation
(2026)*. No DOI or BibTeX record exists yet, so none is given here. To cite the
library itself, see [Citing pymocd](../citing.md).
<!-- TODO(docs): add the full RIMPSO citation, DOI and doi.org BibTeX (as a bibtex code block) once the paper is published; currently only "Santos, in prep. (2026)" appears in README.md and docs/algorithms/index.md. -->

Works cited on this page:

- V. A. Traag, P. Van Dooren and Y. Nesterov. Narrow scope for
  resolution-limit-free community detection. *Physical Review E* 84(1):016114,
  2011. [doi:10.1103/PhysRevE.84.016114](https://doi.org/10.1103/PhysRevE.84.016114)
- K. Deb, A. Pratap, S. Agarwal and T. Meyarivan. A fast and elitist
  multiobjective genetic algorithm: NSGA-II. *IEEE Transactions on Evolutionary
  Computation* 6(2):182–197, 2002.
  [doi:10.1109/4235.996017](https://doi.org/10.1109/4235.996017)
- G. Schwarz. Estimating the dimension of a model. *The Annals of Statistics*
  6(2):461–464, 1978.
  [doi:10.1214/aos/1176344136](https://doi.org/10.1214/aos/1176344136)
- S. Fortunato and M. Barthélemy. Resolution limit in community detection.
  *Proceedings of the National Academy of Sciences* 104(1):36–41, 2007.
  [doi:10.1073/pnas.0605965104](https://doi.org/10.1073/pnas.0605965104)

## See also

- API: [`rimpso`][pymocd.rimpso], [`rimpso_fronts`][pymocd.rimpso_fronts],
  [`rimpso_select`][pymocd.rimpso_select]
- Guide: [Core concepts](../guide/concepts.md),
  [Choosing a detector](../guide/choosing.md),
  [Working with Pareto fronts](../guide/fronts.md),
  [Evaluating partitions](../guide/evaluation.md)
- Gallery: [RIMPSO's resolution profile](../examples/fronts/plot_resolution_profile.md),
  [Choosing a partition from a front](../examples/fronts/plot_choosing_from_front.md),
  [Communities in Zachary's karate club](../examples/basics/plot_karate_club.md)
- The other library method: [HP-MOCD](hpmocd.md);
  all detectors: [Algorithms overview](index.md)
- Full derivation and measurements: the
  [module README](https://github.com/csilab-ufop/pymocd/blob/master/src/core/algorithms/rimpso/README.md)
