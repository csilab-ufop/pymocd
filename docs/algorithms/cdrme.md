---
description: CDRME in pymocd — random-walk community composition, stochastic merge chains under the Eq. (12) linkage objective, and similarity-driven mutation (Dabaghi-Zarandi, Afkhami & Ashoori, JNCA 2025).
---

# CDRME

<span class="pm-pill pm-pill--baseline">Re-implemented baseline</span>
<span class="pm-pill pm-pill--neutral">Single-objective</span>
<span class="pm-pill pm-pill--neutral">Deterministic</span>

CDRME (Community Detection based on Random walk and Multi-objective
Evolutionary algorithm) is a port of
[Dabaghi-Zarandi, Afkhami & Ashoori (2025)](https://doi.org/10.1016/j.jnca.2024.104070),
written from the paper. Softmax-weighted random walks started at
degree-weighted centres compose a primary set of communities; a population of
stochastic agglomerative merge chains diversifies that set under the paper's
linkage objective (Eq. 12); a similarity-driven mutation reassigns weakly
attached nodes; and the candidate with the highest modularity is returned as
`dict[int, int]`. Despite the "multi-objective" in its name, the method
maximizes one scalar: Eq. (12) is the sum of two linkage measures, and, as
the module README notes, the paper has no dominance test or non-dominated
sorting. There is therefore no Pareto front and no `cdrme_fronts`.

<div class="pm-facts" markdown>

| | |
|---|---|
| **Entry points** | [`cdrme`][pymocd.cdrme] |
| **Objectives** | One, maximized: $\mathrm{ObjFunc} = \mathrm{innerLinkage} + \mathrm{outerLinkage}$ (Eq. 12) |
| **Search** | Softmax random walks from degree-weighted centres build a primary partition; `pop_size` independent agglomerative merge chains; similarity-driven mutation of the `elite_size` best chains |
| **Selection rule** | Maximum modularity ($\gamma = 1$); ties broken by Eq. (12), then fewer communities, then the label vector |
| **Front accessor** | none (single objective) |
| **Deterministic** | Yes: bit-identical output on repeated calls with the same graph and for any thread count (module README; tests in `api.rs`). There is no `seed` argument. The result can change if the same graph lists its nodes or edges in a different order (see [Parameters](#parameters)). |
| **Original code** | A private Python notebook supplied by the authors; the module README states that no public repository exists. It is vendored in this repository at [`res/original_algs/cdrme`](https://github.com/csilab-ufop/pymocd/tree/master/res/original_algs/cdrme). |
| **Paper** | Dabaghi-Zarandi et al. (2025), *Journal of Network and Computer Applications* 234:104070 — [DOI](https://doi.org/10.1016/j.jnca.2024.104070) |

</div>

## Objectives

Let $C$ be a partition into $k$ communities. For a community $c$,
$\mathrm{innerLinks}(c)$ counts the links with both ends in $c$,
$\mathrm{totalLinks}(c)$ counts the links with at least one end in $c$, and
$\mathrm{outerLinks}(c, c')$ counts the links between $c$ and $c'$. The single
objective, maximized, is Eq. (12) of the paper:

$$
\mathrm{ObjFunc}(C) = \mathrm{innerLinkage}(C) + \mathrm{outerLinkage}(C)
\tag{12}
$$

with

$$
\mathrm{innerLinkage}(C) = \sum_{c \in C} \frac{\mathrm{innerLinks}(c)}{\mathrm{totalLinks}(c)}
\tag{9}
$$

$$
\mathrm{outerLinkage}(c) = \max_{c' \neq c} \frac{\mathrm{outerLinks}(c, c')}{\mathrm{totalLinks}(c)}
\tag{11}
$$

$$
\mathrm{outerLinkage}(C) = \sum_{c \in C} \Bigl[\, \max_{b \in C} \mathrm{outerLinkage}(b) - \mathrm{outerLinkage}(c) \Bigr]
\tag{10}
$$

The implementation (`objective.rs`) evaluates Eq. (10) in the algebraically
identical one-pass form
$k \cdot \max_b \mathrm{outerLinkage}(b) - \sum_c \mathrm{outerLinkage}(c)$.

!!! note "Both terms grow with the number of communities"
    The module README points out that on the paper's own Fig. 6 example
    $\mathrm{ObjFunc} = 3.97$ at $k = 5$ and exactly $1.0$ at $k = 1$, so a
    merge almost always lowers Eq. (12). This is why the merge acceptance rule
    (below) carries the search, and why the default `elite_size` lets every
    chain reach mutation instead of keeping only the best by Eq. (12).

Degenerate cases are guarded rather than divided by zero:

- a community with $\mathrm{totalLinks}(c) = 0$ contributes $0$ to Eqs. (9)
  and (11);
- at $k = 1$ the maximum in Eq. (11) is over an empty set, so
  $\mathrm{outerLinkage}(C)$ is defined as $0$ and $\mathrm{ObjFunc}$ is
  $\mathrm{innerLinkage}$ alone;
- an empty $\mathrm{primWalk}$ ($\sum_r \mathrm{freq}(r) = 0$ in Eq. 8)
  leaves the similarity column zero except for the centre's own gene;
- a graph with $m = 0$ edges has modularity $0$.

## Search

The search has four stages. Throughout, $N(v)$ is the neighbour set of $v$,
$V_{\mathrm{active}}$ the set of non-isolated nodes and $m$ the number of
edges of the simple graph. "Node order" is the order in which
`graph.nodes()` lists the nodes; the detector indexes nodes in that order.

**1. Pre-processing (Sec. 4.1).** Self-loops and repeated edges are dropped,
so every degree and $m$ refer to the simple graph. Then

$$
\mathrm{AvgDegree}(G) = \frac{2m}{\lvert V_{\mathrm{active}} \rvert}
\tag{5}
$$

$$
\mathrm{ENC} = \operatorname{round}\!\left( \frac{\lvert V_{\mathrm{active}} \rvert}{\mathrm{AvgDegree}(G)} \right)
\tag{6}
$$

clamped to $[1, \lvert V_{\mathrm{active}} \rvert]$. $\mathrm{ENC}$ is the
number of walk rounds, each started from one centre. The primary partition
has one community per distinct centre, so at most $\mathrm{ENC}$ communities
(a node can be drawn as a centre more than once).

**2. Primary composition (Sec. 4.2), $\mathrm{ENC}$ sequential rounds.** Each
round:

1. Draws a centre $x$ with probability proportional to its draw weight,
   which starts at its degree (Sec. 4.2.1) and decays in step 6.

2. Sets the walk length (Eq. 7). Because Eq. (6) sets
   $\mathrm{ENC}$ to $\lvert V \rvert / \mathrm{AvgDegree}(G)$ (before the
   rounding above), the paper's $\lvert V \rvert / \mathrm{ENC}$ term is
   $\mathrm{AvgDegree}(G)$, and the code computes the collapsed form

    $$
    \mathrm{walkLength}(x) = \operatorname{round}\bigl( \mathrm{Degree}(x) + \alpha_{\mathrm{walk}} \cdot \mathrm{AvgDegree}(G) \bigr) \ge 1
    \tag{7}
    $$

3. Runs Algorithm 1: `n_walk` walks of that length from $x$, each step moving
   from $u$ to a neighbour $v$ with the softmax probability

    $$
    P(u \to v) = \frac{\exp \lvert N(u) \cap N(v) \rvert}{\sum_{w \in N(u)} \exp \lvert N(u) \cap N(w) \rvert}, \qquad v \in N(u).
    $$

    The exponent is shifted by the row maximum before exponentiation, the
    accurate softmax form of
    [Blanchard, Higham & Higham (2021)](https://doi.org/10.1093/imanum/draa038),
    which the paper cites for this rule. The $\mathrm{walkLength}(x)$ most
    frequent nodes form $\mathrm{primWalk}$, where $\mathrm{freq}(r)$ is the
    number of walks that visit $r$.

4. Scores every node against the walk (Eq. 8), a frequency-weighted mean of
   the pairwise similarity of Eq. (4). As the code evaluates it,
   Eq. (4) is $\mathrm{sim}(v, u) = \mathrm{connection}(v, u) + \lvert N(v) \cap N(u) \rvert$,
   with $\mathrm{connection}(v, u) = 1$ when $v$ and $u$ are adjacent and $0$
   otherwise:

    $$
    \bar{s}_x(v) = \frac{\sum_{r \in \mathrm{primWalk}} \mathrm{freq}(r)\, \mathrm{sim}(v, r)}{\sum_{r \in \mathrm{primWalk}} \mathrm{freq}(r)}
    \tag{8}
    $$

    The column is then divided by its maximum, so it lies in $[0, 1]$, and
    the centre's own value is set to $1$ (see
    [Divergences](#divergences-from-the-paper)).

5. Folds the column into a single chromosome with one gene
   $(\text{centre}, \text{similarity})$ per node: each gene keeps the centre
   with the highest similarity seen so far (ties go to the centre that comes
   first in node order). This running element-wise arg-max is the paper's
   Sec. 4.2.5 pairwise merge. After the last round, every chosen centre's
   own gene is reset to that centre with similarity $1$.

6. Divides the draw weight of every $\mathrm{primWalk}$ node by its frequency
   (Sec. 4.2.3). Weights decay but are never set to zero, so a node can be
   drawn as a centre twice.

**3. Merge chains (Sec. 4.3).** The primary partition, with $k$ communities,
is converted to a community relation graph (one gene per community: its
internal links, its links to each neighbouring community and its
$\mathrm{totalLinks}$). `pop_size` chains start from it and run in parallel;
chain $s \in \{0, \dots, \texttt{pop\_size} - 1\}$ receives
$\lfloor 2 k s / \texttt{pop\_size} \rfloor$ merge attempts. The chains
therefore spread along the merge ladder, from the primary partition itself
($s = 0$, no attempts) to the deepest chains, with almost $2k$ attempts; the
module README notes that this budget lets the deepest chains run the ladder
down to a single community. A chain stops early once one community is left.
One attempt draws a live community $i$ uniformly, draws a
neighbouring community $j$ of $i$ uniformly, evaluates Eq. (12) for the merged
partition without modifying it, and applies the Sec. 4.3.2 acceptance rule
with $\delta = \mathrm{ObjFunc}(C') - \mathrm{ObjFunc}(C)$:

| $\delta$ | Paper's $1 - 1/\delta$ | pymocd |
|---|---|---|
| $\delta < 0$ | greater than 1 | accept |
| $\delta = 0$ or NaN | undefined | accept |
| $0 < \delta < 1$ | negative | reject |
| $\delta \ge 1$ | in $[0, 1)$ | accept with probability $1 - 1/\delta$ |

The shipped rule is the literal formula clamped into $[0, 1]$, with a merge
that leaves Eq. (12) unchanged ($\delta = 0$) accepted.

**4. Mutation (Sec. 4.4.1–4.4.3).** The chains are ranked by Eq. (12)
(descending, then fewer communities, then the label vector) and the best
`elite_size` are decoded back to the one-gene-per-node form, each keeping
the per-node similarities $s(v)$ of the primary composition. Each of these
chains then sweeps the active nodes in node order, rewriting genes in place.
A node $v$ whose similarity $s(v)$ is below
`alpha_mut` sums its neighbours' similarities per neighbouring community,
$S_c = \sum_{u \in N(v),\, u \in c} s(u)$, out of the total
$T = \sum_{u \in N(v)} s(u)$, and takes the arg-max community $c^\ast$ (ties to
the lower community index). It moves to $c^\ast$, with new similarity
$S_{c^\ast} / T$, only if $T > 0$, $c^\ast$ differs from its current
community and $S_{c^\ast} / T > s(v)$.
Sweeps stop at the first sweep that moves no node, or after `mut_sweeps`
sweeps.

## Selecting a partition

The paper's Sec. 4.4.4 selects the final chromosome by "the quality of
evaluation measures described in Section 2.3.3", which names modularity,
NMI and density. NMI needs ground truth, and density is maximized by the
single-community partition, so the shipped rule is the remaining label-free
measure, Newman–Girvan modularity at resolution $\gamma = 1$ on the simple,
unweighted graph:

$$
Q(C) = \sum_{c \in C} \left[ \frac{L_c}{m} - \left( \frac{d_c}{2m} \right)^2 \right],
$$

where $L_c$ is the number of edges inside $c$ and $d_c$ the sum of the degrees
of its nodes. This is also the measure the authors' notebook sorts on. Among
the `elite_size` mutated candidates the one with the highest $Q$ is returned;
ties go to the higher Eq. (12), then to fewer communities, then to the
lexicographically smaller label vector, so the choice is reproducible.
Community ids in the returned `dict` are dense and assigned in node order
(the first non-isolated node listed by `graph.nodes()` is in community `0`);
isolated nodes are reported as community `-1`.

## Parameters

<!-- pymocd:defaults pymocd.cdrme -->
| Parameter | Default | Description |
|---|---|---|
| `alpha_walk` | `1.0` | Eq. (7) walk-length coefficient. The paper gives "typically 1 to 2"; its worked example (Fig. 3-b) is the $\alpha = 1$ end. Clamped into $[0, 2]$; a non-finite value falls back to the default. |
| `n_walk` | `50` | Walks per centre (Algorithm 1). The paper gives no value. At least 1. |
| `pop_size` | `300` | $N_p$, the number of merge chains. All chains start from the same primary partition, so this sets how many points along the merge ladder are sampled, not the size of a breeding pool. The docstring and module README describe the cost as linear in it. At least 1. |
| `elite_size` | `300` | $N_{sp} \le N_p$, the chains that reach mutation (Sec. 4.4.1), the best by Eq. (12). The default equals the default `pop_size`, so at default settings every chain is kept; it is a fixed `300`, not tied to `pop_size`, so with `pop_size` above 300 only the best 300 chains reach mutation unless `elite_size` is raised. Clamped into $[1, \texttt{pop\_size}]$. |
| `alpha_mut` | `0.5` | Sec. 4.4.2 mutation threshold on the $[0, 1]$ similarity scale: a node below it is offered a new community. The paper gives no value or range. A non-finite value falls back to the default. |
| `mut_sweeps` | `10` | Cap on the Sec. 4.4.2–4.4.3 loop, which the paper leaves unbounded. The loop also stops at the first sweep that moves no node; `0` skips mutation. |

Edge weights are not used: the detector receives the edge list only. There is
no `seed` argument; every random stream is derived from fixed constants, one
independent stream per (stage, slot) pair. Because the draws are made over
the internal node index, listing the same nodes in a different order can give
a different partition; the [Example](#example) shows this on Zachary's karate
club. The edge order can change the node order: a NetworkX graph built from an
edge list lists its nodes in the order the edges first mention them, so the
same edges in another order give another node order (see
[Input order](../guide/performance.md#input-order)). One constant is
internal and not reachable from Python: `MERGE_ATTEMPTS_PER_COMMUNITY = 2`,
which sets the merge budget $2k$ of which chain $s$ receives the fraction
$s / \texttt{pop\_size}$, as described above.

## Complexity

The paper (Dabaghi-Zarandi et al., 2025) is behind the publisher's paywall and
was not read for this page, so whether it states an overall bound has not
been checked and none is quoted here. The module README does not state one
either. It states these per-step costs:

- one Eq. (4) similarity column is accumulated in
  $O\bigl(\sum_{w \in N(r)} \mathrm{deg}(w)\bigr)$ without materializing the
  pairwise similarity matrix;
- the merge phase runs on the community relation graph, so applying a merge
  costs $O(\mathrm{deg}(c))$ on that graph rather than $O(m)$ on the network,
  and a rejected trial merge costs nothing to undo;
- the cost is linear in `pop_size`; the README reports 6.0 s at
  `pop_size=100` and 14.3 s at `pop_size=300` on a 50k-node, 241k-edge
  stochastic block model with 48 cores.

The bounds below are derived from the Rust source of this implementation.
Write $n$ for the number of nodes, $m$ for the number of edges of the simple
graph, $\Delta$ for the largest degree, $\ell$ for the longest walk (Eq. 7),
$k \le \mathrm{ENC}$ for the number of primary communities, $E_C$ for the
number of community pairs joined by at least one edge in the primary
partition, $P$ for `pop_size`, $E = \min(\texttt{elite\_size}, P)$, $S$ for
`mut_sweeps`, $W$ for `n_walk` and $T$ for the number of Rayon threads.

- **Pre-processing** (`topology.rs`) sorts and deduplicates every adjacency
  row: $O(n + m \log \Delta)$.
- **Primary composition** (`primary.rs`), $\mathrm{ENC}$ rounds. Each round
  passes over all $n$ nodes a fixed number of times: to clear, normalize and
  rescale the Eq. (8) column (`avg_similarity` in `similarity/walk_sim.rs`)
  and to fold it into the chromosome (`Chromosome::absorb` in
  `chromosome/genes.rs`). Each round also runs $W$ walks of at most $\ell$
  steps, each step a binary search in a cached transition row, and ranks the
  visited nodes by frequency, $O(W \ell \log n)$; it then accumulates
  Eq. (8) over at most $\ell$ walk nodes, each a two-hop scan of
  $O(\Delta^2)$ entries. The softmax rows (`walk/procedure.rs`) are built at
  most once per node, $O\bigl(\sum_v \mathrm{deg}(v)^2\bigr)$ in all. The
  stage costs
  $O\bigl(\mathrm{ENC} \cdot (n + W \ell \log n + \ell \Delta^2) + \sum_v \mathrm{deg}(v)^2\bigr)$.
- **Merge chains** (`evolve/merge.rs`). Each chain builds its own relation
  graph and decodes its labels, $O(n + m)$. Chain $s$ makes
  $\lfloor 2 k s / P \rfloor$ attempts, about $k P$ over all chains, and every
  attempt scores the trial merge with `objective_if_merged`
  (`objective.rs`), which reads the relation-graph row of every live
  community, $O(k + E_C)$. Applying an accepted merge is cheaper, as the
  README states. The stage costs $O\bigl(P (n + m) + P k (k + E_C)\bigr)$.
- **Mutation and selection** (`evolve/mutate.rs`, `evolve/select.rs`,
  `search` in `api.rs`). Each of the $E$ candidates runs at most $S$ sweeps of
  $O(n + m)$, then a modularity evaluation, a relation-graph build and a sort
  of its labels: $O\bigl(E (S (n + m) + n \log n)\bigr)$.

The whole run is the sum of the four stages. For a graph of fixed mean degree
$\bar d = 2m/n$, $\mathrm{ENC} \approx n / \bar d$, so the
$\mathrm{ENC} \cdot n$ term of the composition grows as $n^2 / \bar d$; and
since $k$ and $E_C$ also grow in proportion to $n$, the merge stage grows as
$P n^2$. The run is therefore quadratic in $n$ and linear in `pop_size`.
Memory is $O\bigl(P n + T (n + k + E_C) + m\bigr)$: the label vectors of all
$P$ chains (4 bytes per node each) are held together until they are ranked;
each running chain or mutation task holds its own relation graph, and each
mutation task a copy of the similarity vector (8 bytes per node); the cached
transition rows hold one 8-byte value per adjacency entry, at most $16m$
bytes.

The bounds are total work. Pre-processing and the $\mathrm{ENC}$ composition
rounds are sequential (each round's decay feeds the next draw) and do not
depend on `pop_size`, so the $n^2 / \bar d$ term is not divided among
threads. The `pop_size` merge chains and the `elite_size` mutation tasks run
in parallel on the Rayon pool (`into_par_iter` in `api.rs`), so the wall time
of those stages is about their work divided by $\min(P, T)$ and
$\min(E, T)$.

The scaling was checked on jackal (48 cores, shared with other jobs; Python
3.10; pymocd 6.1.0 from PyPI) on LFR graphs from
`pymocd.lfr_graph(n, 2.5, 1.5, 0.3, average_degree=10, max_degree=50, min_community=20, max_community=100, seed=1)`,
which have $m \approx 5n$ and no isolated nodes. Each setting ran in a fresh
interpreter, on one thread (`pymocd.max_cores(1)`), and only the `cdrme` call
was timed. With `pop_size=1` and `mut_sweeps=0` the single chain makes no
merge attempt and mutation is skipped, so the returned partition is the
primary one and the time covers only pre-processing, composition and
linear-time bookkeeping:

| $n$ | 2,500 | 5,000 | 10,000 | 20,000 | 40,000 | 80,000 |
|---|---|---|---|---|---|---|
| $\mathrm{ENC}$ (Eq. 6) | 251 | 490 | 1,003 | 1,995 | 3,989 | 7,982 |
| $k$ (primary) | 251 | 489 | 1,002 | 1,993 | 3,983 | 7,969 |
| $E_C$ (primary) | 3,950 | 8,437 | 17,345 | 34,935 | 70,489 | 141,864 |
| Time, `pop_size=1`, `mut_sweeps=0` (s) | 0.014 | 0.035 | 0.11 | 0.33 | 1.18 | 4.96 |
| Time, `pop_size=30` (s) | 0.11 | 0.40–0.43 | 1.24–1.31 | 5.06–5.36 | 22.2 | — |

$k$ is within 0.2% of $\mathrm{ENC}$ and $E_C$ is 16 to 18 times $k$ on
these graphs, so both grow in proportion to $n$, as assumed above. The
composition time grows by 2.5, 3.1, 3.1, 3.6 and 4.2 times per doubling of
$n$, rising toward the factor of 4 of a quadratic term as the
$\mathrm{ENC} \cdot n$ passes overtake the linear terms. With 30 chains,
where the merge stage dominates, the time grows by 2.9 to 4.4 times per
doubling and by 4.1 to 4.4 times from 20,000 to 40,000 nodes (two runs per
setting, except a single run at $n = 40{,}000$). At $n = 10{,}000$,
`pop_size` = 1, 25, 50,
100 and 200 took 0.10, 1.11–1.19, 2.27–2.52, 4.30–4.36 and 9.01–9.22 s,
0.042 to 0.046 s per added chain: linear, as derived. At $n = 20{,}000$
(one run each), `pop_size` = 50, 200 and 400 took 8.6, 34.9 and 67.1 s,
0.17 s per added chain, about four times the per-chain cost at half the
size. At the defaults on the $n = 10{,}000$ graph, one thread took
12.2–13.0 s and 48 threads took 0.58–0.63 s of wall time (17.9–18.5 s of CPU
time), with other jobs keeping the machine's load average at about 18.

The memory terms were checked from the peak resident set size during the
call minus the resident size before it. On one thread at $n = 20{,}000$, the
call added 9.4, 20.2 and 36.4 MiB at `pop_size` = 50, 200 and 400, or 4.0
bytes per node per extra chain, against the 4 bytes of a chain's label
vector. On 48 threads the per-thread term dominates: at $n = 50{,}000$ the
call added 641 MiB at `pop_size=100` and 723 MiB at `pop_size=400` (one run
each). This measurement does not separate the relation graphs from the
memory allocator's own per-thread overhead.

## When to use it

- As a baseline for the CDRME method of Dabaghi-Zarandi et al. (2025), with
  the divergences listed below, in comparisons with the library's own methods
  ([RIMPSO](rimpso.md), [HP-MOCD](hpmocd.md)) and the other baselines.
- When you need one deterministic partition per graph: the module README notes
  that, because the detector is deterministic, one run per graph is exact.
  That holds for a fixed node and edge order; listing the nodes or edges in
  another order can change the result. Run-to-run variability cannot be
  sampled by changing a seed.
- Not when you need a Pareto front or a resolution profile: CDRME returns a
  single partition. In best-in-front (oracle) comparisons it can contribute
  only that partition; see [Working with Pareto fronts](../guide/fronts.md).
- Edge weights are ignored, so a weighted graph is treated as unweighted.

??? info "Results reported in the module README"
    One run at the shipped defaults on 48 cores. $k$ is the number of
    communities of non-isolated nodes; NMI and AMI are against each network's
    ground truth, $Q$ is modularity.

    | Network | $n$ | $m$ | $k$ | NMI | AMI | $Q$ | Time |
    |---|---|---|---|---|---|---|---|
    | karate | 34 | 78 | 2 | 0.837 | 0.833 | 0.371 | <0.01 s |
    | dolphins | 62 | 159 | 6 | 0.395 | 0.369 | 0.438 | <0.01 s |
    | polbooks | 105 | 441 | 3 | 0.406 | 0.394 | 0.404 | <0.01 s |
    | football | 115 | 613 | 11 | 0.857 | 0.814 | 0.560 | <0.01 s |
    | email-EU | 1,005 | 16,064 | 25 | 0.531 | 0.469 | 0.316 | 0.02 s |
    | SBM 5k | 5,000 | 24,094 | 49 | 0.445 | 0.410 | 0.346 | 0.15 s |
    | SBM 50k | 50,000 | 241,088 | 139 | 0.436 | 0.388 | 0.350 | 14.3 s |

    The README adds that on karate $k = 2$ is the exact two-club split: the NMI
    of 0.837 is the fixed value between the two ground-truth conventions in use
    for that network, not a partial recovery. The README does not record the
    node and edge order each network was loaded in; since the partition
    depends on them (see [Parameters](#parameters)), another loader can give
    other values.

## Divergences from the paper

CDRME was written from the paper. The authors' notebook was consulted only
where the text is silent and is not followed where it contradicts the text.
The module README records the following departures.

### From the paper

- **$\mathrm{ENC}$ cancels out of Eq. (7).** Since
  $\mathrm{ENC} = \lvert V \rvert / \mathrm{AvgDegree}(G)$, Eq. (7) collapses
  to $\mathrm{Degree}(v) + \alpha_{\mathrm{walk}} \cdot \mathrm{AvgDegree}(G)$;
  $\mathrm{ENC}$ survives only as the number of walk rounds.
- **Eq. (8) is rescaled onto $[0, 1]$.** Eq. (8) is an unbounded count whose
  scale depends on the density of the centre's neighbourhood, yet the paper's
  Figs. 4 and 5 show every gene in $[0, 1]$ with each centre at $1.0$, and the
  Sec. 4.4.2 threshold needs a bounded scale. Each walk's column is divided by
  its own maximum.
- **The acceptance probability is clamped into $[0, 1]$.** The literal
  $1 - 1/\delta$ exceeds 1 for a worsening merge and is negative for a gain
  below 1. Both branches the paper's stop condition needs are kept; the middle
  branch (a gain below 1 is rejected) is kept exactly as written rather than
  inverted.
- **The merge schedule is stratified, not drawn.** Sec. 4.3.2 hands each
  iteration to a uniformly random chromosome; pymocd gives chain $s$ a fixed
  $s / \texttt{pop\_size}$ share of the same total budget, so the whole merge
  ladder is sampled for every $k$ and the chains are independent (hence
  parallel and independent of the thread count).
- **Sec. 4.2.5's pairwise merge is a running arg-max fold**, the same operator
  applied incrementally.
- **Mutation order and iteration count are chosen.** The paper fixes neither;
  nodes are rewritten in place in node order, and the loop stops at a
  fixed point or after `mut_sweeps` sweeps.
- **A node whose neighbours are all in its own community is left alone**,
  since Fig. 7 describes a membership change; rewriting such a gene in place
  would pin its similarity at $1.0$ and make it immune to mutation.
- **`elite_size` defaults to the default `pop_size` (300)**, so at default
  settings the Sec. 4.4.1 filter is effectively off. Neither the size of the
  filtered set nor "acceptable"
  $\mathrm{ObjFunc}$ is defined in the paper, and ranking by Eq. (12) would
  drop the coarse chromosomes, because Eq. (12) grows with $k$.
- **Empty communities count in the $k$ of Eq. (10)** when Eq. (12) is
  re-evaluated after mutation (mutation can empty a community). This affects
  only the tie-break of the final selection, never the merge chains.
- **Isolated nodes are reported as `-1`** and excluded from $\lvert V \rvert$
  in Eqs. (5) and (6); the paper does not mention them.
- **Self-loops and repeated edges are dropped**, because Eq. (4) uses
  simple-graph quantities.
- **Of the paper's three Sec. 2.3.3 evaluation measures, only modularity is
  used inside the detector**; see
  [Selecting a partition](#selecting-a-partition).
- **Parameter values the paper does not give** (`n_walk`, `alpha_mut`,
  `mut_sweeps`, `elite_size`) were chosen as described in
  [Parameters](#parameters).

### From the authors' notebook

The notebook is vendored at
[`res/original_algs/cdrme`](https://github.com/csilab-ufop/pymocd/tree/master/res/original_algs/cdrme),
converted to plain Python with its code cells verbatim and its outputs
dropped. Four of its choices are deliberately not reproduced:

- **Fitness.** Its `kkm_rc()` fitness is not the paper's objective and is
  degenerate: both of its terms equal $0.5$ for every partition with at least
  one link inside a community and one link between communities, so on all
  such partitions the summed fitness is the constant $1.0$ and carries no
  search signal. Eqs. (9)–(12) are implemented instead.
- **Similarity.** It uses
  $\alpha \lvert N(v) \cap N(u) \rvert / \lvert N(v) \rvert + (1 - \alpha)\, \mathrm{connection}(v, u)$,
  which is not Eq. (4). Eq. (4) is implemented.
- **Walker.** Its walks are uniform over neighbours; Algorithm 1 line 10 is
  softmax-weighted by the common-neighbour count. Algorithm 1 is implemented.
- **Number of centres.** It uses 10% of the nodes; Eq. (6) is
  $\lvert V \rvert / \mathrm{AvgDegree}(G)$. Eq. (6) is implemented.

Where the notebook and the paper agree it is followed: its final selection by
modularity resolves the paper's unspecified Sec. 4.4.4 rule. For `n_walk`, the
README notes that the notebook runs 500 walks of fixed length 4, while Eq. (7)
makes each walk $\mathrm{Degree}(v) + \mathrm{AvgDegree}(G)$ steps long, so 50
walks cover a comparable number of steps per centre.

## Example

CDRME on Zachary's karate club. Modularity is computed with NetworkX on the
unweighted graph (`weight=None`), which is what the detector's selection rule
uses; NetworkX's karate graph carries edge weights that CDRME ignores.

```python exec="on" source="above" session="cdrme" result="text"
import networkx as nx
import pymocd

G = nx.karate_club_graph()
partition = pymocd.cdrme(G)

communities = {}
for node, community in partition.items():
    communities.setdefault(community, set()).add(node)

print("communities:", len(communities))
print("sizes:", sorted(len(c) for c in communities.values()))
print("modularity:", round(nx.community.modularity(G, communities.values(), weight=None), 4))
print("same partition on a second call:", partition == pymocd.cdrme(G))
```

The last line checks the determinism stated above: a second call with the same
arguments returns the identical partition.

The same edges with the nodes listed in reverse order form an equal graph,
but the detector indexes nodes in `graph.nodes()` order, so its random draws
fall on different nodes and the partition can differ:

```python exec="on" source="above" session="cdrme" result="text"
R = nx.Graph()
R.add_nodes_from(sorted(G.nodes(), reverse=True))
R.add_edges_from(G.edges())

reordered = pymocd.cdrme(R)
same_edges = {frozenset(e) for e in G.edges()} == {frozenset(e) for e in R.edges()}
print("same nodes and edges:", set(G) == set(R) and same_edges)
print("communities (reversed node order):", len(set(reordered.values())))
print("same partition as the default order:", reordered == partition)
```

When results must be reproduced, keep both the node order and the edge order
fixed, for example by adding the nodes in sorted order and then the edges in
sorted order.

## Reference

Fahimeh Dabaghi-Zarandi, Mohammad Mehdi Afkhami and Mohammad Hossein Ashoori.
"Community Detection method based on Random walk and Multi objective
Evolutionary algorithm in complex networks." *Journal of Network and Computer
Applications* 234:104070, February 2025.
[doi:10.1016/j.jnca.2024.104070](https://doi.org/10.1016/j.jnca.2024.104070)

```bibtex
@article{Dabaghi_Zarandi_2025,
  title={Community Detection method based on Random walk and Multi objective Evolutionary algorithm in complex networks},
  volume={234},
  ISSN={1084-8045},
  url={http://dx.doi.org/10.1016/j.jnca.2024.104070},
  DOI={10.1016/j.jnca.2024.104070},
  journal={Journal of Network and Computer Applications},
  publisher={Elsevier BV},
  author={Dabaghi-Zarandi, Fahimeh and Afkhami, Mohammad Mehdi and Ashoori, Mohammad Hossein},
  year={2025},
  month=Feb,
  pages={104070}
}
```

## See also

- API: [`pymocd.cdrme`][pymocd.cdrme]
- [Choosing a detector](../guide/choosing.md) and
  [Core concepts](../guide/concepts.md) (multi-objective detection and
  selection)
- [Input graphs](../guide/input-graphs.md) (node ids, isolated nodes,
  self-loops) and [Performance and reproducibility](../guide/performance.md)
- [Evaluating partitions](../guide/evaluation.md)
- [GDPSO](gdpso.md), the other single-objective baseline, and the
  [algorithms overview](index.md)
- Gallery: [Detector tour](../examples/basics/plot_detector_tour.md)
- The module README, with the full derivation:
  [`src/core/algorithms/cdrme/README.md`](https://github.com/csilab-ufop/pymocd/blob/master/src/core/algorithms/cdrme/README.md)
