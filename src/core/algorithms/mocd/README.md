# Shi-MOCD — Multi-objective community detection by PESA-II

**Paper**
Shi, C., Yan, Z., Cai, Y. & Wu, B., "Multi-objective community detection in
complex networks", *Applied Soft Computing* 12(2):850–859, 2012.
DOI 10.1016/j.asoc.2011.10.005.

Engine: Corne, D. W., Jerram, N. R., Knowles, J. D. & Oates, M. J., "PESA-II:
Region-based selection in evolutionary multiobjective optimization",
*Proc. GECCO 2001*, pp. 283–290. Shi §3 drives PESA-II, not NSGA-II — the
region-based archive is the part that differs, and it is reimplemented here in
`pesa2/`.

Representation: Park, Y. & Song, M., "A genetic algorithm for clustering
problems", *Proc. 3rd Annual Conf. on Genetic Programming*, pp. 568–575, 1998 —
the locus-based adjacency encoding Shi §3.1.3 adopts.

**Original implementation**
The first author's publications page
(http://www.shichuan.org/html/publications_ch.html) links a `[code]` archive,
`MOGAForClustering28.zip`, next to this paper, along with an article-in-press
copy of the paper (http://www.shichuan.org/doc/11.pdf). The archive is a
Visual C++ project whose `GA` class keeps the IP/EP populations, offers Q's two
components as objectives, and implements both the Max-Q and the Max–Min
Distance selection. This module was written from the paper text; the breeding
step (`GA::runGA`, `GA::mutate`, `GA::crossover2`, `GA::setMaxCluster`) and the
grid resolution (`NICHSIZE`) follow that code where the paper leaves them open.
Every equation below is taken from the paper text; every choice the paper
leaves open is listed under **Divergences**.

**Objectives**
Notation over `G = (V, E)` with `m = |E|` (`graph.edges` holds each undirected
edge exactly once, so `edges.len() == m`), a partition `C = {c_1 … c_k}`,
`l_c` the edges with both endpoints in `c` counted once, and `d_c = Σ_{i∈c}
deg(i)` (which counts every internal edge twice).

Shi Eqs. 3.5/3.6, both **minimized**:

    intra = 1 − Σ_c l_c / m           fraction of edges NOT internal to a community
    inter = Σ_c (d_c / 2m)^2          the degree-concentration (null) term

so that Newman modularity is the exact complement (Eq. 3.7):

    Q = 1 − intra − inter

The names are Shi's and they are easy to misread: `intra` is *one minus* the
internal-edge fraction, so minimizing it maximizes internal edges, and `inter`
is the null term that a single-community partition drives to 1. Fragmenting
lowers `inter` and raises `intra`; coarsening does the reverse. That is the
whole trade-off, and it is why the Pareto front spans the resolution range.

`calculate_objectives` walks node *positions*, accumulating `l_c` and `d_c` into
two `Vec` accumulators indexed by community id, and divides once at the end, so
the float value is reproducible for a given labelling. Internal edges are
counted through `NodeIndex::neighbor_candidates`, which stores both directions
of every edge; the `pos < neighbor` guard is what reduces that to one count per
edge, and it is also what drops the self-allele that `neighbor_candidates`
carries for an isolated node. Removing the guard double-counts `l_c` *and*
gives isolated nodes a phantom internal edge.

**Representation**
Locus-based adjacency (Park & Song 1998; Shi §3.1.3). A `Genome = Vec<usize>`
of length `n`: gene `p` holds a *position* `q`, read as "link position `p` to
position `q`". Positions index `NodeIndex::index_to_node`, which is
`graph.nodes_vec()` — sorted by `Graph::finalize`, so the position order is
fixed for a given graph. `NodeId`s appear only at the module boundary, where
`api.rs::labels_to_partition` zips the dense label vector back against
`nodes_vec()` by position.

`NodeIndex::neighbor_candidates[p]` is the set of legal alleles for position
`p`: its neighbours, or `[p]` itself when it has none. Initialization draws
from that set and both operators redraw from it, so **every genome is legal by
construction and no repair pass exists**. The self-allele is the only
self-referential gene and only an isolated node ever holds it; such a node
becomes a singleton community and `normalize_community_ids` reports it as `-1`.

`decode` is a union-find over the `n` links followed by a relabelling pass that
assigns community ids `0..k` in ascending-position first-seen order, so two
genomes encoding the same partition decode to the identical label array. Shi
describes a backtracking decode; the components are the same set either way.

**Algorithm**
PESA-II keeps two populations: an internal population `IP` of size `ipsize`
that is regenerated every cycle, and an external archive `EP` of capacity
`epsize` that holds the non-dominated set found so far. `EP` is both the result
and the mating pool.

1. `NodeIndex::build` freezes the node order and the neighbour-position lists;
   per-position degrees are read once into a dense `Vec`.
2. `seed_archive` — `ipsize` random legal genomes, decoded and evaluated into
   `[intra, inter]`; each is offered to `EP` through `insert_nondominated`, then
   `EP` is truncated to `epsize`.
3. Per generation `0 .. num_gens−1`, stopping early if `EP` empties:
   a. `assign_cells` bins `EP` into a `GRID_DIVISIONS^2` hyper-grid, every axis
      normalized to the range `EP` *currently* spans, and rebuilds the per-cell
      occupancy `occ` — the squeeze factors. The grid is adaptive: it is rebuilt
      from scratch each time `EP` changes.
   b. `breed`, as the authors' `GA::runGA` does it: `ipsize` clones of archive
      members, each drawn by `select_member` and mutated once; then
      `ceil(ipsize / 2)` variation steps, each of which `mutate`s one random
      clone with probability `mut_rate` or else, with probability
      `cross_rate`, runs `uniform_crossover` in place on two distinct random
      clones. A mutation redraws `mutation_size(n)` random genes from their
      neighbours (`GA::setMaxCluster`: n/2, n/3, n/5, n/50, n/200 or n/1000 for
      n below 20, 50, 100, 1000, 10000 and above).
   c. Evaluate the children and offer each to `EP` via `insert_nondominated`.
   d. If `EP` overflows, `truncate` back to `epsize`.
4. `evolutionary_phase` returns `EP` as the Pareto front.

`select_member` is the region-based selection of the authors' code: a
uniformly random **occupied cell**, then a uniformly random member of it. A cell
with one member is drawn as often as a cell with twenty, which is what spreads
the parents along the front. (Corne et al. 2001 hold a binary tournament
between two occupied cells instead.)

`truncate` is the mirror image: repeatedly pick a most-crowded cell (ties
uniformly), drop a uniformly random member of it, and **rebuild the grid**
before the next removal, until `|EP| == epsize`. `insert_nondominated` rejects a
candidate only if some `EP` member *dominates* it, so objective-wise duplicates
accumulate freely; truncation is what bounds them.

Model selection — the two entry points:

- **MOCD-Q** (Eq. 3.8) — the front member of maximum `Q`. Since
  `Q = 1 − intra − inter`, that is `argmin(intra + inter)`. This selector lives
  in `api/detectors.rs::mocd_q_fn`, **not** in this module; the module only
  produces the front via `Mocd::generate_pareto_front`.
- **MOCD-D** (Eqs. 3.9–3.11) — `Mocd::run`. Generate `rand_networks`
  Erdős–Rényi control graphs of the same scale, evolve a full front on each,
  and return the real front member whose `(intra, inter)` point is furthest from
  all of them: for each real solution take the minimum Euclidean distance to any
  control-front member (over every control front), then keep the solution whose
  minimum is largest. The idea is that a random graph's front is what "no
  community structure" looks like in objective space, so the most deviant real
  point is the most structurally meaningful one.

**Parameters**
`Mocd::new` (constructor arguments; `Mocd` itself is not registered with the
Python module — see **Known issues**):

| Parameter | Default | Effect |
|---|---|---|
| `debug_level` | 0 | `0` silent; `≥1` prints the graph once and the `EP` size every generation. |
| `rand_networks` | 3 | MOCD-D control fronts. Read only by `run()`; ignored by `generate_pareto_front()`. Each one costs a **full** evolutionary run. |
| `pop_size` | 100 | PESA-II `ipsize`. Floored at 1 internally. |
| `num_gens` | 100 | Generations. `0` returns the seeded archive. |
| `cross_rate` | 0.6 | Probability that a variation step which does not mutate crosses two clones. |
| `mut_rate` | 0.4 | Probability that a variation step mutates one clone (every clone is also mutated once when copied). |

`cross_rate = 0.6` / `mut_rate = 0.4` are Shi §4.2.2, introducing Table 1:
*"pc and pm are 0.6 and 0.4 for these four EA based algorithms"* — fixed across
every experiment in the paper (§4.3). `rand_networks = 3` is Shi §3.2, which
generates three control fronts for the artificial networks (the timed
real-network runs use one, §4.2.2).

The Python wrappers `pymocd.mocd_q` / `pymocd.mocd_d` use the same rates.
(Up to 6.1.0 they passed 0.9 and 0.1, the MOCD-baseline setting of the HP-MOCD
benchmark of Santos et al. 2025; pass `cross_rate=0.9, mut_rate=0.1` for that
setting.) To reproduce a column of Shi's Table 1, pass its
`ipsize`/generation counts as `pop_size`/`num_gens`, and `rand_networks=1` for
MOCD-D on these networks:

| Network (Table 1) | n | m | `ip` = `pop_size` | `ep` | `gen` = `num_gens` |
|---|---:|---:|---:|---:|---:|
| Karate (P1) | 34 | 78 | 50 | 50 | 50 |
| Lesmis (P2) | 77 | 254 | 50 | 50 | 50 |
| Polbooks (P3) | 105 | 441 | 50 | 50 | 100 |
| Adjnoun (P4) | 112 | 425 | 50 | 50 | 100 |
| Football (P5) | 115 | 613 | 50 | 50 | 100 |
| Celegansneural (P6) | 297 | 2345 | 100 | 100 | 100 |
| Celegansmetabolic (P7) | 453 | 2025 | 100 | 100 | 100 |
| Netscience (P8) | 1589 | 2742 | 200 | 100 | 200 |
| Power (P9) | 4941 | 6594 | 300 | 100 | 300 |
| Hepth (P10) | 8361 | 15751 | 400 | 100 | 400 |

`n`/`m` are the counts the table prints (it spells P6 "Celegansnearal"). `ep`
equals `min(ip, 100)` in every column, so the derived `ep_size` below matches
the table without an extra parameter.

Internal constants a reader will hit:

| Constant | Value | Where | Meaning |
|---|---|---|---|
| `EPSIZE_CAP` | 100 | `defaults.rs` | Upper bound on the derived `ep_size`. |
| `MOCD_D_RAND_NETWORKS` | 3 | `defaults.rs` | Control-front count for the `mocd_d` wrapper. |
| `GRID_DIVISIONS` | 10 | `pesa2/grid.rs` | Hyper-grid bins per objective axis, so 100 cells for two objectives (`NICHSIZE` of the authors' code). |

`epsize` is **not** a public parameter. It is derived once in `Mocd::new` as
`pop_size.clamp(1, EPSIZE_CAP)`, i.e. `min(pop_size, 100)` with a floor of 1.
Shi Table 1 assigns `ipsize`/`epsize`/generations per graph and tops out at
`epsize = 100`, though §4.1.2 runs `epsize = 200`; the cap encodes the table,
not §4.1.2, and there is no way to reach 200 through the public API.

**Determinism**
**Nondeterministic: two consecutive calls can differ.** `evolutionary_phase`
and `generate_random_networks` both draw from `rand::rng()`, the thread-local
generator, which is seeded from the OS and is neither exposed nor re-seedable
through the public API; there is no seed parameter. Everything downstream is
order-fixed — a sorted node order, `Vec` accumulators, first-seen relabelling,
a single-threaded loop, no rayon anywhere — so the RNG stream is the only source
of run-to-run variation.

Measured with the defaults on karate over 30 runs (AMI against networkx's
`club` split): `mocd_q` AMI 0.5667 ± 0.0000, k = 4.00 ± 0.00; `mocd_d` AMI
0.6161 ± 0.0592, k = 2.47 ± 0.51. MOCD-D is the noisier of the two, because its
answer depends on `rand_networks` freshly drawn random graphs as well as on its
own front. On florentine over 20 runs `mocd_d` returns k ∈ {2, 3, 4} (13/4/3)
while `mocd_q` returns k = 3 every time — the spread is the control fronts, not
the search.

The *selected* partition is often stable on small graphs even though the search
is not: over 10 calls, `mocd_q` returned one partition on florentine and on a
two-clique graph with isolated nodes, two partitions on karate (both with 4
communities and the same AMI) and ten different ones on lesmis. That is the
objective landscape of a small graph, not determinism.

**Divergences**
This module was written from the paper text; where the text is silent on the
breeding step it follows the authors' code (see **Original implementation**).
These are deliberate readings of, or departures from, the paper text:

- **The breeding step follows the authors' code.** The paper says only that
  mutation reassigns "some genes" and that `pm` is "the ratio of mutation". Up
  to 6.1.0 this module read `pm` as a per-gene probability, which at Shi's
  `pm = 0.4` redraws 40% of every child and stalls the search (on the paper's
  four-group benchmark, MOCD-Q identified 0.70 of the vertices correctly at
  z_out = 5 against the paper's 1.0). `GA::runGA` instead mutates every clone
  once with `GA::setMaxCluster(n)` genes, then applies one operator per pair
  of children: mutation with probability `pm`, otherwise crossover with
  probability `pc`. With that step the same benchmark gives 0.92.
- **"Uniform two-point crossover" is implemented as per-gene uniform
  crossover**, as the authors' `GA::crossover2` does it: a random mask swaps
  each position between the two clones in place. An allele kept at its own
  position is legal in either parent, so no repair is needed.
- **The archive is PESA-II's, not the authors' `GA::update`.** That function
  uses a dominance test that is strict in both objectives and, on a full
  archive, swaps the candidate for a member of a more crowded cell. This module
  keeps Corne et al.'s rule (reject a dominated candidate, truncate an overfull
  archive from its most crowded cell); on the four-group benchmark the two
  gave the same accuracy.
- **Union-find decode instead of the paper's backtracking.** Same connected
  components, so the same partition; only the identifiers differ, and those are
  compacted to first-seen order anyway.
- **The MOCD-D null model is Erdős–Rényi `G(n, m)`**, uniformly random simple
  edges over the same node set with the same edge count — Shi's "random networks
  with the same scale" (§3.2). A degree-preserving null (double-edge-swap) was
  considered and rejected: it is a stronger control, but it is not what the
  paper describes, and MOCD-D's decision rule is calibrated against the weaker
  one. All `n` nodes are pre-inserted into the control graph before any edge is
  added, so an isolated node in the original stays a node of the control graph
  and `n` matches; dropping that pre-insertion silently shrinks the control
  graph and shifts the whole control front.
- **Control fronts get the FULL evolutionary budget.** Each of the
  `rand_networks` control runs uses the same `num_gens`/`pop_size`/`ep_size` as
  the real run, so `mocd_d` costs `rand_networks + 1` full searches. This looks
  like an obvious place to economize and is not: an under-evolved random front
  never reaches the fragmented region, which makes the *real* front's fragmented
  extreme look spuriously deviant and hands MOCD-D the wrong answer.
- **No shared engine, no rayon.** PESA-II is reimplemented inside this module
  and runs single-threaded, rather than reusing the project's shared NSGA-II /
  NSGA-III machinery, so this baseline's runtime tracks the published
  single-threaded method rather than the repository's optimizations. Hot-path
  state is dense `Vec`-indexed; the shared `Partition` map is built only at the
  API boundary.
- **`epsize` is derived, not exposed** (see **Parameters**).
- **MOCD-Q lives outside this module.** The max-`Q` argmin is in
  `api/detectors.rs`; the module surface stops at the front. MOCD-D, which needs
  the control-front machinery, is `Mocd::run` here.

**Known issues**
Recorded, deliberately **not** fixed — behaviour is frozen for this module.

- **`mocd_q` raises `RuntimeError("MOCD-Q produced an empty front")` on an empty
  graph and on a single-node graph.** `evolutionary_phase` returns an empty
  archive whenever the graph has no nodes *or* no edges, and the wrapper turns
  the empty front into an error. Any edgeless graph, not just those two, hits it.
- **`mocd_d` *panics* on the same inputs.** `min_max_selection` ends in
  `.expect("Real Pareto front is empty.")`, which surfaces in Python as
  `pyo3_runtime.PanicException`. That is not the `RuntimeError` its sibling
  raises, and because `PanicException` derives from `BaseException` a caller's
  `except Exception` will not catch it.
- **`rand_networks = 0` degenerates silently.** With no control fronts, the fold
  over an empty `random_fronts` leaves every solution scoring `f64::MAX`, the
  strict `>` never fires again, and `run()` returns the *first* member of the
  front — an arbitrary partition, not a selected one. Not reachable from
  `pymocd.mocd_d` (default 3), but reachable through `Mocd::new(.., 0, ..)`.
- **`Mocd` is a `#[pyclass]` that is never registered** with the Python module
  in `lib.rs`, so it cannot be constructed from Python. `pymocd.mocd_q` and
  `pymocd.mocd_d` are the only entry points, and there is no `mocd_fronts`
  exposing `generate_pareto_front` the way the other detectors expose theirs.
- **`Mocd::envolve` is a misspelling of "evolve".** Left alone because the
  method is part of the module's recorded Rust surface.
- **`calculate_objectives`'s `total_edges == 0` early return is unreachable**
  from the live path — `evolutionary_phase` rejects edgeless graphs before any
  evaluation happens. It is defensive only.

**Files**

| Path | Holds |
|---|---|
| `mod.rs` | Module root: submodule wiring and the `Mocd` / defaults re-export. |
| `api.rs` | The `Mocd` pyclass: construction, `generate_pareto_front` (the MOCD-Q candidate set), `run` (MOCD-D), and the dense-labels → `Partition` boundary conversion. |
| `defaults.rs` | The shipped defaults (Shi's rates) and the `EPSIZE_CAP`. |
| `locus.rs` | `NodeIndex` (node order + legal alleles), `Genome`, random initialization, union-find `decode`. |
| `objectives.rs` | `calculate_objectives` → `Metrics { intra, inter }`, Shi Eqs. 3.5/3.6. |
| `operators.rs` | In-place uniform crossover, neighbour-restricted mutation and `mutation_size` (the authors' `GA::setMaxCluster`). |
| `null_model.rs` | `generate_random_networks`: the Erdős–Rényi `G(n, m)` control graphs MOCD-D scores against. |
| `selection.rs` | `min_max_selection`, the MOCD-D decision rule (Eqs. 3.9–3.11). |
| `pesa2/mod.rs` | PESA-II sub-module root. |
| `pesa2/solution.rs` | `Solution` (labels + objectives) and its dominance test, `Member` (genome + solution + grid cell), and `evaluate`. |
| `pesa2/grid.rs` | The adaptive hyper-grid: `GRID_DIVISIONS`, `assign_cells`, `select_member`. |
| `pesa2/archive.rs` | External-archive maintenance: `insert_nondominated` and squeeze-factor `truncate`. |
| `pesa2/engine.rs` | `evolutionary_phase`, the generational loop, plus archive seeding and `breed`. |
