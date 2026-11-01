---
description: Shi-MOCD in pymocd — a port of Shi et al.'s PESA-II community detector over the decomposed-modularity pair (intra, inter), with the MOCD-Q (max modularity) and MOCD-D (max–min distance to random-graph fronts) selection rules.
---

# Shi-MOCD (MOCD-Q / MOCD-D)

<span class="pm-pill pm-pill--baseline">Re-implemented baseline</span>
<span class="pm-pill pm-pill--neutral">Multi-objective</span>

Shi-MOCD is the multi-objective community detection method of
[Shi, Yan, Cai & Wu (2012)](https://doi.org/10.1016/j.asoc.2011.10.005). It
minimizes the two terms into which Newman's modularity $Q$ splits,
$\mathrm{intra}$ and $\mathrm{inter}$, with the PESA-II evolutionary algorithm
over a locus-based adjacency encoding, and keeps the non-dominated partitions
in an external archive. The paper proposes two rules for picking one partition
from that front, and pymocd exposes each as its own function:
[`mocd_q`][pymocd.mocd_q] returns the front member of highest modularity
(MOCD-Q), and [`mocd_d`][pymocd.mocd_d] returns the member that lies furthest
from the fronts of Erdős–Rényi control graphs (MOCD-D). Both return a single
partition; the front itself is not exposed.

<div class="pm-facts" markdown>

| | |
|---|---|
| **Entry points** | [`mocd_q`][pymocd.mocd_q] (MOCD-Q), [`mocd_d`][pymocd.mocd_d] (MOCD-D) |
| **Objectives** | $\mathrm{intra}$ and $\mathrm{inter}$ (Shi et al., 2012, Eqs. 3.5 and 3.6), both minimized; $Q = 1 - \mathrm{intra} - \mathrm{inter}$ (Eq. 3.7) |
| **Search** | PESA-II (Corne et al., 2001) with an adaptive hyper-grid, over the locus-based adjacency encoding (Park & Song, 1998) |
| **Selection rule** | MOCD-Q: maximum $Q$ (Eq. 3.8). MOCD-D: maximum of the minimum distance to Erdős–Rényi control fronts (Eqs. 3.9–3.11) |
| **Front accessor** | None; there is no `mocd_fronts` |
| **Deterministic** | No: there is no seed parameter, and the random generator is seeded from the operating system (module README, *Determinism*) |
| **Original code** | The first author's [publications page](http://www.shichuan.org/html/publications_ch.html) links a C++ archive, `MOGAForClustering28.zip`, as the code of this paper. The port was written from the paper text. Its variation operators, parent selection and grid resolution, which the paper leaves open, follow that code (see [Divergences from the paper](#divergences-from-the-paper)) |
| **Paper** | Shi et al. (2012), *Applied Soft Computing* 12(2):850–859 — [DOI](https://doi.org/10.1016/j.asoc.2011.10.005); the same publications page links an [article-in-press copy](http://www.shichuan.org/doc/11.pdf) |

</div>

!!! note "The paper's rates, a per-graph budget"
    [`mocd_q`][pymocd.mocd_q] and [`mocd_d`][pymocd.mocd_d] default to
    `cross_rate=0.6` and `mut_rate=0.4`, the crossover and mutation rates of
    Shi et al. (2012), which "are fixed for all experiments" (§4.3). The paper
    sets the internal population, the archive size and the number of
    generations per graph in its Table 1, so to reproduce one of its runs, set
    `pop_size` and `num_gens` to the values listed in
    [The paper's settings](#shi-table-1). `pop_size` sets both the internal
    population and, capped at 100, the archive size (see
    [Parameters](#parameters)); every column of Table 1 follows that rule.
    Up to 6.1.0 the defaults were `cross_rate=0.9` and `mut_rate=0.1`, the
    setting of the MOCD baseline in the HP-MOCD benchmark of
    [Santos et al. (2025)](https://doi.org/10.1007/s13278-025-01519-7); pass
    them explicitly for that setting.

## Search

**Representation.** Shi-MOCD uses the locus-based adjacency encoding of
Park & Song (1998), which Shi et al. adopt in their §3.1.3. A genome holds one
gene per node; gene $i$ names one neighbor of node $i$, read as "link node $i$
to that neighbor". The communities are the connected components of the graph
formed by these links. A node without neighbors points at itself and becomes a
singleton. Initialization and mutation draw every gene from the node's
neighbors, and crossover copies each gene from one parent at the same
position, so the module README notes that every genome is valid by
construction and no repair step exists. Decoding is a
union-find over the $n$ links followed by a relabeling of the communities in
order of first appearance.

**PESA-II.** The search is PESA-II
([Corne, Jerram, Knowles & Oates, 2001](#reference)), reimplemented inside the
module. It keeps an internal population of size `pop_size` that is
regenerated every generation and an external archive of capacity
`min(pop_size, 100)` that holds the non-dominated partitions found
so far. The archive is both the result and the mating pool.

1. **Seeding.** `pop_size` random genomes are decoded and evaluated; each is
   offered to the archive, which is then truncated to its capacity.
2. **Generations.** For each of the `num_gens` generations (stopping early if
   the archive is empty):
    1. The archive is binned into a hyper-grid with 10 divisions per
       objective (100 cells), each axis normalized to the range the archive
       currently spans. The occupancy of a cell is its *squeeze factor*.
    2. `pop_size` children are bred as the authors' implementation breeds
       them. Each child starts as a copy of an archive member chosen by
       region-based selection (below) and is mutated once. A mutation redraws
       a fixed number of genes, each uniformly from its node's neighbors
       (which may return the same neighbor); the number depends only on the
       node count $n$: $\lfloor n/2 \rfloor$ below 20 nodes,
       $\lfloor n/3 \rfloor$ below 50, $\lfloor n/5 \rfloor$ below 100,
       $\lfloor n/50 \rfloor$ below 1,000, $\lfloor n/200 \rfloor$ below
       10,000 and $\lfloor n/1000 \rfloor$ above (11 genes on karate, 2 on a
       128-node graph, 24 on the power grid). Then $\lceil P/2 \rceil$
       variation steps follow, with $P$ = `pop_size`: with probability
       `mut_rate` a step mutates one random child again; otherwise, with
       probability `cross_rate`, it applies uniform crossover to two distinct
       random children in place, swapping each gene between them with
       probability 1/2, so that both offspring are kept.
    3. Each child is evaluated and offered to the archive. It is rejected only
       if an archive member dominates it; archive members it dominates are
       removed.
    4. If the archive exceeds its capacity, it is truncated.
3. **Result.** The final archive is the Pareto front.

**Region-based selection.** Parents are chosen as in the authors'
implementation: a uniformly random occupied cell, then a uniformly random
member of it. Selection is over cells, not members: a cell with one member is
drawn as often as a cell with twenty, which spreads the parents along the
front. (Corne et al. (2001) instead hold a binary tournament between two
occupied cells on their squeeze factors.)

**Truncation.** While the archive is over capacity, a most crowded cell is
chosen (ties uniformly at random), a uniformly random member of it is removed,
and the grid is rebuilt before the next removal. Since only dominated
candidates are rejected, members with identical objective values can
accumulate; truncation is what bounds them.

## Selecting a partition

Both functions run the PESA-II search above on the input graph and apply one
of Shi's two model-selection rules to the resulting front $F$. Write
$f(C) = (\mathrm{intra}(C), \mathrm{inter}(C))$ for the objective point of a
partition.

**MOCD-Q** ([`mocd_q`][pymocd.mocd_q], Shi Eq. 3.8) returns the front member
of highest modularity. Since $Q = 1 - \mathrm{intra} - \mathrm{inter}$, this
is

$$
C^{\star} = \operatorname*{arg\,max}_{C \in F} Q(C)
          = \operatorname*{arg\,min}_{C \in F} \bigl( \mathrm{intra}(C) + \mathrm{inter}(C) \bigr).
$$

**MOCD-D** ([`mocd_d`][pymocd.mocd_d], Shi Eqs. 3.9–3.11) compares the front
with fronts evolved on random graphs. It generates $r$ = `rand_networks`
Erdős–Rényi $G(n, m)$ control graphs with the same node set and the same
number of edges as the input, evolves a full front $R_1, \dots, R_r$ on each,
and returns the real front member that is furthest from all of them:

$$
D(C) = \min_{1 \le j \le r} \; \min_{C' \in R_j} \bigl\lVert f(C) - f(C') \bigr\rVert_2,
\qquad
C^{\star} = \operatorname*{arg\,max}_{C \in F} D(C).
$$

The module README explains the idea: a random graph's front shows what "no
community structure" looks like in objective space, so the real point that
deviates most from it is taken as the most meaningful one. The equations above
describe the rule as implemented (`selection.rs`); the module README gives the
same rule in words and attributes it to Shi Eqs. 3.9–3.11.

The result is a `dict[int, int]` mapping each node to a community numbered
from 0. Isolated nodes are reported as community `-1`.

!!! note "Each call is a new search"
    Every call runs the search from scratch, and the search is not seeded. A
    call to `mocd_q` and a call to `mocd_d` on the same graph therefore select
    from the fronts of two independent searches, and `mocd_d` also draws new
    control graphs on every call. The module README reports that over 20 runs
    on the Florentine families graph, `mocd_d` returned between 2 and 4
    communities while `mocd_q` returned 3 every time, and attributes the
    spread to the control fronts rather than to the search.

!!! info "No front accessor"
    The Rust class behind both functions can produce the whole front, but the
    module README notes that it is not registered with the Python module, and
    there is no `mocd_fronts`. To study the $(\mathrm{intra}, \mathrm{inter})$
    front itself, [HP-MOCD](hpmocd.md) optimizes the same pair and exposes its
    front through [`hpmocd_fronts`][pymocd.hpmocd_fronts] and
    [`HpMocd.generate_pareto_front`][pymocd.HpMocd.generate_pareto_front]; see
    [Working with Pareto fronts](../guide/fronts.md).

## Parameters

One table covers both functions; `rand_networks` belongs to `mocd_d` only.

<!-- pymocd:defaults pymocd.mocd_q pymocd.mocd_d -->
| Parameter | Default | Description |
|---|---|---|
| `pop_size` | `100` | PESA-II internal population size: the number of random genomes in the seed population and of children bred per generation. Also sets the archive capacity to `min(pop_size, 100)`. Values below 1 are raised to 1. |
| `num_gens` | `100` | Number of generations. `0` returns a selection from the seeded archive. |
| `cross_rate` | `0.6` | Probability that a variation step that does not mutate applies uniform crossover to two random children. Must lie in $[0, 1]$. Shi et al.'s value. |
| `mut_rate` | `0.4` | Probability that a variation step mutates one random child, on top of the mutation every child receives when it is copied from the archive. A step does nothing with probability $\max(0, 1 - \texttt{cross\_rate} - \texttt{mut\_rate})$. Must lie in $[0, 1]$. Shi et al.'s value. |
| `rand_networks` | `3` | `mocd_d` only. Number of Erdős–Rényi control graphs; each one costs a full evolutionary run with the same budget. The module README attributes the value 3 to Shi §3.2, which uses three control fronts for the artificial networks (and one for the timed real-network runs). |

The archive capacity (`epsize` in PESA-II) is not a parameter: it is derived
as `min(pop_size, 100)` with a floor of 1. The module README states that the
cap of 100 encodes the largest archive in Shi's Table 1, and that the archive
of 200 used in the paper's §4.1.2 cannot be reached through the public API.
The number of grid divisions per objective is fixed at 10.

The `cross_rate` and `mut_rate` defaults are the paper's; `pop_size` and
`num_gens` are per graph in the paper (below).

### The paper's settings {#shi-table-1}

Shi et al. (2012) compare MOCD-Q and MOCD-D with other detectors on ten real
networks, with the per-graph settings of their Table 1 (§4.2.2), reproduced
below together with the node and edge counts the table gives. Every run uses
$p_c = 0.6$ and $p_m = 0.4$, and on these networks the Max–Min Distance rule
uses one random control front (§4.2.2), that is `rand_networks=1`. The archive
size `ep` equals `min(ip, 100)` in every column, the rule pymocd uses to
derive it, so `pop_size` = `ip` and `num_gens` = `gen` reproduce each column.

| Network | Nodes | Edges | `ip` (`pop_size`) | `ep` | `gen` (`num_gens`) |
|---|---:|---:|---:|---:|---:|
| Karate (P1) | 34 | 78 | 50 | 50 | 50 |
| Lesmis (P2) | 77 | 254 | 50 | 50 | 50 |
| Polbooks (P3) | 105 | 441 | 50 | 50 | 100 |
| Adjnoun (P4) | 112 | 425 | 50 | 50 | 100 |
| Football (P5) | 115 | 613 | 50 | 50 | 100 |
| Celegansneural (P6) | 297 | 2,345 | 100 | 100 | 100 |
| Celegansmetabolic (P7) | 453 | 2,025 | 100 | 100 | 100 |
| Netscience (P8) | 1,589 | 2,742 | 200 | 100 | 200 |
| Power (P9) | 4,941 | 6,594 | 300 | 100 | 300 |
| Hepth (P10) | 8,361 | 15,751 | 400 | 100 | 400 |

The table spells P6 "Celegansnearal". For example, the paper's MOCD-D run on
Netscience corresponds to
`pymocd.mocd_d(G, pop_size=200, num_gens=200, cross_rate=0.6, mut_rate=0.4, rand_networks=1)`.
The paper's other experiments give their budgets in the text: `ip` = `ep` =
100 with 800 generations for the hierarchical network (§4.1.1), 200 and 500
for the overlapping network (§4.1.2), and 100 and 200 for the four-group
artificial networks (§4.2.1).

### Degenerate inputs

The module README records the first two behaviors below as known issues that
are left unfixed on purpose. Every row was checked against the current
package.

| Input | Result |
|---|---|
| A graph with no edges (including an empty or single-node graph) | `mocd_q` raises `RuntimeError("MOCD-Q produced an empty front")`. `mocd_d` panics with "Real Pareto front is empty."; in Python this is `pyo3_runtime.PanicException`, which derives from `BaseException`, so `except Exception` does not catch it. |
| `rand_networks=0` | There are no control fronts, every member gets the same score, and `mocd_d` returns the first member of the archive: an arbitrary partition, not a selected one. |
| `cross_rate` or `mut_rate` outside $[0, 1]$ | Both functions panic (`pyo3_runtime.PanicException`, message `p=… is outside range [0.0, 1.0]`). |
| `num_gens=0` | The seeded archive is returned as the front, and the selection rule is applied to it. |

## Complexity

The module README states that both functions run single-threaded, with no
parallelism anywhere in the module, so the runtime of this baseline tracks the
published single-threaded method rather than the repository's optimizations.
It also states that `mocd_d` costs `rand_networks + 1` full searches: the real
one plus one per control graph, each with the same `pop_size` and `num_gens`.

Below, $n$ and $m$ are the numbers of nodes and edges, $P$ is `pop_size`,
$E = \min(P, 100)$ is the archive capacity, $G$ is `num_gens` and $r$ is
`rand_networks`. The hyper-grid has a fixed 100 cells. File names are relative
to `src/core/algorithms/mocd/` unless given in full.

**Time of this implementation.** One generation of `evolutionary_phase`
(`pesa2/engine.rs`) breeds and evaluates $P$ children. Breeding them
(`breed`) copies $P$ archive members, each found by one draw
(`select_member`, `pesa2/grid.rs`) that scans the archive, and runs at most
$\lceil P/2 \rceil$ crossovers over all $n$ genes; mutations
(`operators.rs`) redraw at most $n/2$ genes each. One child costs
$O(n + E)$. Evaluating it
(`evaluate`, `pesa2/solution.rs`) decodes the genome with a union-find over
the $n$ links (`decode`, `locus.rs`) and computes both objectives by visiting
every node and both directions of every edge (`calculate_objectives`,
`objectives.rs`), so it costs $O(n + m)$. Offering each child to the archive
(`insert_nondominated`, `pesa2/archive.rs`) scans at most $E + P$ members,
and `truncate` removes at most $P$ members, rebuilding the grid and rescanning
the archive after each removal; together these cost $O(P(E + P))$ per
generation and never read the graph. Seeding the archive costs about one
generation, and indexing the graph (`NodeIndex::build`, `locus.rs`) costs
$O(n + m)$ once. Since $E \le P$, one search costs

$$
O\bigl((G + 1)\,P\,(n + m + P)\bigr).
$$

`mocd_q` runs one search, converts every front member to a node-to-community
map (`generate_pareto_front`, `api.rs`), which costs $O(E\,n)$, and returns
the member of maximum $Q$ (`mocd_q_fn`, `src/api/detectors.rs`). `mocd_d` runs
$r + 1$ searches. Drawing each Erdős–Rényi control graph (`null_model.rs`)
takes an expected $O(m)$ edge draws on a sparse graph, and the selection rule
(`min_max_selection`, `selection.rs`) compares every real front member with
every control-front member, $O(r E^2)$, so `mocd_d` takes about $r + 1$ times
as long as `mocd_q`.

**Memory of this implementation.** Every child and every archive member holds
an $n$-gene genome (`Vec<usize>`, 8 bytes per gene) and, once evaluated, an
$n$-entry label vector (`Vec<i32>`, 4 bytes per node): about $12n$ bytes. At
most $E + P$ of them are alive at once, next to the $O(n + m)$ graph and node
index, so one search needs $O(P\,n + m)$ memory, and the front conversion of
`mocd_q` another $O(E\,n)$. `mocd_d` also holds all $r$ control graphs and the
label vectors of $r + 1$ fronts, $O(r\,(n + m + E\,n))$ more.

**What the paper states.** Shi et al. (2012, §4.3) count the evaluation of one
individual as $O(m + n)$ ($O(m)$ for the objectives, $O(n)$ for the decode),
the PESA-II process as $O(g s^2)$, citing Corne et al. (2001), and the whole
method as $O(g s^2 (m + n))$, "which grows linearly with the scale of the
network", where $g$ is the number of generations and $s$ is both `ipsize` and
`epsize`. They give Max Q as $O(1)$ and Max–Min Distance as
$O(g s^2 (m + n))$, and report that MOCD-D, with one control front, takes
nearly twice as long as MOCD-Q (§4.2.2). This implementation has the paper's per-individual cost, except that
its decode is a union-find with path halving and no union by rank instead of
the paper's linear-time backtracking (§3.1.3). Its overall bound is tighter
than the paper's, because the archive work and the graph work add rather than
multiply: the archive operations compare objective pairs and grid cells and
never read the graph. Both bounds are linear in $n + m$ at fixed $P$ and $G$.
The paper states no memory bound.

**Measured.** The bound was checked on jackal (48 cores, shared; Python
3.10.12; pymocd 6.1.0 from PyPI). Each call ran in a fresh process on the LFR
graph
`pymocd.lfr_graph(n, 2.5, 1.5, 0.3, average_degree=10, max_degree=50, min_community=20, max_community=100, seed=1)`,
and each time below is the median of three calls. Memory is the rise of the
process's peak resident set (`ru_maxrss`) over its resident set just before
the call, in MB of $10^6$ bytes.

| `mocd_q` sweep | Wall time | Log–log slope |
|---|---|---|
| $n$ = 1,000 → 32,000 ($m$ = 4,887 → 160,335), defaults | 0.380 s → 13.5 s | 1.02 against $n + m$ |
| `pop_size` = 25 → 800, $n$ = 2,000, `num_gens=100` | 0.189 s → 6.38 s | 1.03 |
| `num_gens` = 25 → 400, $n$ = 2,000, `pop_size=100` | 0.189 s → 3.28 s | 1.03 |

Across the $n$ sweep, the time per bred and evaluated child stayed between 6.2
and 6.9 ns per unit of $n + m$. The time per unit of `pop_size` stayed between
7.4 and 8.4 ms up to `pop_size=800`, so the $P^2$ archive term did not show at
these sizes. A `num_gens=0` call, which converts the graph, seeds the archive
and converts the front, took under 2% of the time of the default run at every
$n$. At $n$ = 2,000, `mocd_d` took 1.88, 4.00 and 7.16 times as long as
`mocd_q` with `rand_networks` = 1, 3 and 6, and at $n$ = 8,000 it took 4.26
times as long with the default 3. Process CPU time matched wall time to within
2.4% in every call, as expected of a single-threaded engine. Peak memory grew
linearly with $n$, from 7.8 MB at $n$ = 4,000 to 60.9 MB at $n$ = 32,000 (1.9
to 2.0 kB per node at `pop_size=100`; below $n$ = 4,000 the rise was as small
as its run-to-run spread), and at $n$ = 2,000 each unit of
`pop_size` between 200 and 800 added 11.4 bytes per node, close to the 12
bytes per node of one genome and its labels.

## When to use it

- **As a published baseline.** Shi-MOCD is one of the re-implemented
  baselines (ports of other authors' published methods), and the two
  functions implement Shi's two selection rules on the same search. For new
  analyses, start from [Choosing a detector](../guide/choosing.md).
- **To follow the paper, pass its budget.** The rates default to Shi's;
  the population and the number of generations are per graph in the paper
  (see the note at the top of this page).
- **MOCD-Q or MOCD-D.** MOCD-Q needs one search and selects by modularity.
  MOCD-D selects by distance from a random-graph reference rather than by
  modularity, but costs `rand_networks + 1` searches, and the module README
  reports it as the noisier of the two rules (see
  [Selecting a partition](#selecting-a-partition)).
- **Not when you need the front.** Neither function returns it. HP-MOCD
  optimizes the same objective pair and exposes its front.
- **Not when you need reproducible output.** There is no seed parameter, so
  repeated calls can return different partitions. [RIMPSO](rimpso.md), one of
  the library's own methods, returns identical output on repeated calls; see
  [Performance and reproducibility](../guide/performance.md).
- **Not on graphs without edges**, where both functions fail (see
  [Degenerate inputs](#degenerate-inputs)).

## Divergences from the paper

The port was written from the paper text, not from the authors' code (see the
*Original code* row at the top of this page), and the module README lists the
following as deliberate readings of, or departures from, that text:

- **Variation operators from the authors' code.** The paper describes
  mutation only as "randomly select some genes and assign them with other
  randomly selected adjacent nodes" and calls $p_m$ "the ratio of mutation".
  Up to 6.1.0 the port read $p_m$ as a per-gene probability, which at the
  paper's $p_m = 0.4$ redraws 40% of every child and keeps the search from
  converging. The breeding step now follows the authors' C++ implementation
  (see *Original code*): the number of genes a mutation redraws, the copy and
  mutate step, the order of the mutation and crossover steps, uniform
  selection over occupied cells and a 10 × 10 grid. On the paper's
  four-group benchmark (§4.2.1: 128 nodes, $z_{\mathrm{in}} +
  z_{\mathrm{out}} = 16$; `mocd_q` with $p_c = 0.6$, $p_m = 0.4$,
  `pop_size=100`, `num_gens=200`; 50 graphs per point), the mean fraction of
  vertices identified correctly rose from 0.80 to 0.98 at
  $z_{\mathrm{out}} = 4$ and from 0.70 to 0.92 at $z_{\mathrm{out}} = 5$,
  where the paper reports 1.0. On the hierarchical network of §4.1.1, the
  maximum-$Q$ member rose from $Q = 0.634$ to $0.696$ (paper: $0.706$ with
  four communities; mean over 10 graphs).
- **Uniform crossover.** Shi names the operator "uniform two-point crossover"
  but describes it as an independent per-gene choice between the two parents,
  and the authors' code implements uniform crossover. A gene inherited at its
  own position is valid in either parent, so genomes stay valid.
- **PESA-II's archive, not the authors' archive update.** The authors' code
  rejects a candidate only if a member is better in both objectives and, when
  the archive is full, lets a candidate replace a member of a more crowded
  cell. The port keeps the archive of Corne et al. (2001): a candidate is
  rejected if a member dominates it, and an overfull archive is truncated
  from its most crowded cell. On the benchmark above the two archives gave
  the same accuracy.
- **Union-find decoding** instead of the paper's backtracking. The connected
  components, and therefore the partition, are the same.
- **Erdős–Rényi $G(n, m)$ null model.** MOCD-D's control graphs are uniformly
  random simple graphs with the same node set and edge count, Shi's "random
  networks with the same scale" (§3.2). The module README records that a
  degree-preserving null model was considered and rejected because it is not
  what the paper describes. All $n$ nodes, including isolated ones, are
  inserted into each control graph, so its node count matches the input.
- **Full budget for the control fronts.** Each control search uses the same
  `pop_size` and `num_gens` as the real one. The module README explains that
  an under-evolved random front never reaches the fragmented region, which
  makes the real front's fragmented extreme look spuriously deviant and leads
  MOCD-D to the wrong answer.
- **A separate, single-threaded engine.** PESA-II is reimplemented inside the
  module instead of reusing the repository's shared NSGA-II/NSGA-III code.
- **Derived archive size and benchmark rates.** The archive size is not
  exposed (see [Parameters](#parameters)), and the Python functions ship the
  HP-MOCD benchmark rates instead of Shi's.

## Example

The example runs both selection rules on Zachary's karate club at the
defaults and MOCD-Q again at the paper's rates. For each result it prints the
number of communities $k$, the two objectives computed from the partition with
NetworkX, and $1 - \mathrm{intra} - \mathrm{inter}$ next to NetworkX's own
modularity, which should agree (Eq. 3.7). Edge weights are ignored, matching
pymocd.

```python exec="on" source="above" session="shi-mocd" result="text"
import networkx as nx
import pymocd

G = nx.karate_club_graph()
m = G.number_of_edges()

def describe(name, partition):
    groups = {}
    for node, community in partition.items():
        groups.setdefault(community, set()).add(node)
    intra = 1 - sum(G.subgraph(c).number_of_edges() for c in groups.values()) / m
    inter = sum((sum(d for _, d in G.degree(c)) / (2 * m)) ** 2 for c in groups.values())
    q = nx.community.modularity(G, groups.values(), weight=None)
    print(f"{name:<24} k = {len(groups)}   intra = {intra:.4f}   inter = {inter:.4f}"
          f"   1 - intra - inter = {1 - intra - inter:.4f}   Q = {q:.4f}")

describe("mocd_q (defaults)", pymocd.mocd_q(G))
describe("mocd_d (defaults)", pymocd.mocd_d(G))
# Shi et al. (2012) rates; the paper's Table 1 also sets pop_size and num_gens per graph.
describe("mocd_q (pc=0.6, pm=0.4)", pymocd.mocd_q(G, cross_rate=0.6, mut_rate=0.4))
```

MOCD-D draws new control graphs on every call. The next block repeats it five
times and prints the number of communities of each result; the counts need
not agree.

```python exec="on" source="above" session="shi-mocd" result="text"
print([len(set(pymocd.mocd_d(G).values())) for _ in range(5)])
```

## Reference

C. Shi, Z. Yan, Y. Cai and B. Wu, "Multi-objective community detection in
complex networks", *Applied Soft Computing*, vol. 12, no. 2, pp. 850–859,
2012.
[doi:10.1016/j.asoc.2011.10.005](https://doi.org/10.1016/j.asoc.2011.10.005)

```bibtex
@article{Shi_2012,
  title={Multi-objective community detection in complex networks},
  volume={12},
  ISSN={1568-4946},
  url={http://dx.doi.org/10.1016/j.asoc.2011.10.005},
  DOI={10.1016/j.asoc.2011.10.005},
  number={2},
  journal={Applied Soft Computing},
  publisher={Elsevier BV},
  author={Shi, Chuan and Yan, Zhenyu and Cai, Yanan and Wu, Bin},
  year={2012},
  month=Feb,
  pages={850–859}
}
```

Other works cited on this page:

- D. W. Corne, N. R. Jerram, J. D. Knowles and M. J. Oates, "PESA-II:
  Region-based selection in evolutionary multiobjective optimization", in
  *Proc. Genetic and Evolutionary Computation Conference (GECCO 2001)*,
  pp. 283–290, 2001.
- Y. Park and M. Song, "A genetic algorithm for clustering problems", in
  *Proc. 3rd Annual Conference on Genetic Programming*, pp. 568–575, 1998.
- M. E. J. Newman and M. Girvan, "Finding and evaluating community structure
  in networks", *Physical Review E*, vol. 69, no. 2, 026113, 2004.
  [doi:10.1103/PhysRevE.69.026113](https://doi.org/10.1103/PhysRevE.69.026113)
- G. O. Santos, L. S. Vieira, G. Rossetti, C. H. G. Ferreira and
  G. J. P. Moreira, "A high-performance evolutionary multiobjective community
  detection algorithm", *Social Network Analysis and Mining*, vol. 15, no. 1,
  110, 2025.
  [doi:10.1007/s13278-025-01519-7](https://doi.org/10.1007/s13278-025-01519-7).
  Preprint: [arXiv:2506.01752](https://arxiv.org/abs/2506.01752).

## See also

- API: [`mocd_q`][pymocd.mocd_q] and [`mocd_d`][pymocd.mocd_d]
- [HP-MOCD](hpmocd.md), which optimizes the same objective pair with NSGA-II
  and exposes its front
- [Algorithms overview](index.md)
- [Core concepts](../guide/concepts.md), [Choosing a detector](../guide/choosing.md)
  and [Working with Pareto fronts](../guide/fronts.md)
- [Input graphs](../guide/input-graphs.md) and
  [Performance and reproducibility](../guide/performance.md)
- Gallery: [Detector tour](../examples/basics/plot_detector_tour.md) and
  [Comparing detectors](../examples/benchmarks/plot_compare_detectors.md)
- The module README, with the implementation notes, known issues and file
  layout:
  [`src/core/algorithms/mocd/README.md`](https://github.com/csilab-ufop/pymocd/blob/master/src/core/algorithms/mocd/README.md)
