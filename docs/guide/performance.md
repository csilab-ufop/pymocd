---
description: How pymocd uses threads, how to size the Rayon thread pool with max_cores, which detectors are deterministic, how RIMPSO's seed works, the measured time and memory of each detector, and how to time a detector.
---

# Performance and reproducibility

Four detectors run in parallel on a Rayon thread pool that is global to the
process and sized once; the other five run on a single thread. No pymocd
function releases Python's global interpreter lock (GIL), so concurrent runs
need separate processes. Three detectors are deterministic: RIMPSO, whose
output depends only on the graph as given (the order in which it lists its
nodes and edges included) and its parameters, `seed` included, and not on the
thread count; and CDRME and GDPSO, which return identical output on every call
for the same graph object and take no seed. Listing the same nodes or edges in
a different order can change the result of all three (see
[Input order](#input-order)). The other six draw from a random generator the
caller cannot seed, and their output can change from call to call. This page
gives the exact semantics of [`max_cores`][pymocd.max_cores], a per-detector
determinism table with the source of each statement, the measured time and
peak memory of every detector at three graph sizes, and a recipe for timing a
detector.

## Threads

### Which detectors are parallel

| Runs on the Rayon pool | Runs on one thread |
|---|---|
| RIMPSO, HP-MOCD, CDRME, GDPSO | MMCoMO, NSGA-III-CCM, NSGA-III-KRM, Shi-MOCD, MOGA-Net |

Every entry point of the four parallel detectors uses the pool, including
[`rimpso_fronts`][pymocd.rimpso_fronts], [`rimpso_select`][pymocd.rimpso_select],
[`hpmocd_fronts`][pymocd.hpmocd_fronts] and the [`HpMocd`][pymocd.HpMocd] class.
The module READMEs of the five single-threaded baselines state that their
search is sequential, with no Rayon (MMCoMO's: no parallelism at all). The
NSGA-III-CCM, NSGA-III-KRM, Shi-MOCD and MOGA-Net READMEs add that this is
deliberate, so that the baseline's cost tracks the published method rather
than the library's infrastructure. The metric functions
([`gt_metrics`][pymocd.gt_metrics] and the others) also run on one thread.

### Sizing the pool with `max_cores`

[`max_cores(num_threads)`][pymocd.max_cores] builds Rayon's global thread pool
with `num_threads` threads. Its exact behavior, from
`src/core/utils/threading.rs` and tested with the current package:

- **The first call in a process wins.** The pool is built inside a
  once-only guard. On that first call the function prints two lines to the
  process's standard output, tagged `WARNING`: `Global thread pool initialized
  with N threads` and `Using max_cores again has no effect, due to static
  ThreadPoolBuilder initialization`. Both lines are written by Rust, so they
  bypass `sys.stdout` (and `contextlib.redirect_stdout`), and they are printed
  by the first call, not by later ones.
- **Later calls do nothing**, silently, as long as the first call succeeded.
- **`max_cores(0)` gives Rayon's default size.** Rayon treats a thread count
  of 0 as "choose automatically", so the pool gets one thread per logical CPU,
  or `RAYON_NUM_THREADS` threads when that variable is set, although the log
  line reports `0 threads`.
- **Call it before any parallel detector.** If you do not call `max_cores`,
  the first parallel detector to run builds Rayon's default pool, which has
  one thread per logical CPU, or the number given by the `RAYON_NUM_THREADS`
  environment variable when that is set
  ([Rayon documentation](https://docs.rs/rayon/latest/rayon/struct.ThreadPoolBuilder.html)).
  A `max_cores` call after that point raises `pyo3_runtime.PanicException`
  (Rayon reports `GlobalPoolAlreadyInitialized`), and so does every later
  call. `PanicException` derives from `BaseException`, so `except Exception`
  does not catch it. The detectors keep working on the default pool.
  Running only single-threaded detectors does not build the pool.

A first call made after a parallel detector has already run raises instead.
The block below shows that case in a fresh process, so that it does not touch
the pool of the process that builds this page:

```python exec="on" source="above" session="performance" result="text"
import subprocess
import sys

LATE_CALL = """
import networkx as nx
import pymocd

pymocd.rimpso(nx.karate_club_graph())  # builds Rayon's default pool
for attempt in (1, 2):
    try:
        pymocd.max_cores(2)
        print(f"max_cores call {attempt}: no error")
    except BaseException as exc:  # PanicException is not an Exception subclass
        print(f"max_cores call {attempt}: {type(exc).__name__}: {exc}")
"""

run = subprocess.run([sys.executable, "-c", LATE_CALL], capture_output=True, text=True)
print(run.stdout.strip())
```

Because the pool cannot be resized, a program that needs two thread counts
needs two processes; the repository's own test of thread-count independence
(`tests/test_determinism.py`) starts one child process per thread count.

### The GIL

No pymocd function releases the GIL: the Rust sources contain no call to
PyO3's `allow_threads` or `detach`. A detector therefore holds the GIL for its
whole call: other Python threads in the process cannot run Python code while
the Rust code works, and several detectors started from a thread pool such as
`concurrent.futures.ThreadPoolExecutor` do not run their Rust computation at
the same time.

To run many detections at once (several graphs, or repeated runs of a
nondeterministic detector), use processes. Each process has its own Rayon
pool, built at its first parallel call with one thread per logical CPU (or
`RAYON_NUM_THREADS`) unless you size it, so several worker processes running
parallel detectors can start many more threads than there are CPUs. Pass
`initializer=pymocd.max_cores` and
`initargs=(threads_per_worker,)` to `concurrent.futures.ProcessPoolExecutor` to
size each worker's pool before it runs anything.

!!! warning "Start method"
    Use the `spawn` or `forkserver` start method, for example
    `mp_context=multiprocessing.get_context("spawn")`, when the parent process
    has already run a parallel detector. In a test with the current package on
    Linux, workers created with `fork` (the default start method on Linux
    before Python 3.14; see the
    [multiprocessing documentation](https://docs.python.org/3/library/multiprocessing.html))
    after the parent had run [`rimpso`][pymocd.rimpso] hung without an
    initializer. With the `max_cores` initializer, the call raised
    `PanicException` in each worker and the executor raised
    `BrokenProcessPool`. The same workloads completed under `spawn` and
    `forkserver`, with and without the initializer.

## Determinism

The table records, for each detector, whether repeated calls with the same
graph object and parameters return the same partition, whether the thread
count can change the result, and whether a seed is exposed. Each row is taken
from the *Determinism* section of the module README linked in the last column;
the "fixed input order" qualifier comes from the tests under
[Input order](#input-order).

| Detector | Functions | Repeated calls | Thread count | Seed | Source |
|---|---|---|---|---|---|
| RIMPSO | `rimpso`, `rimpso_fronts`, `rimpso_select` | identical for a [fixed input order](#input-order) | no effect | `seed` | [README §8](https://github.com/csilab-ufop/pymocd/blob/master/src/core/algorithms/rimpso/README.md), docstring, `tests/test_determinism.py` |
| HP-MOCD | `hpmocd`, `hpmocd_fronts`, `HpMocd` | can differ | Rayon scheduling is one listed source | none | [README](https://github.com/csilab-ufop/pymocd/blob/master/src/core/algorithms/hpmocd/README.md) |
| CDRME | `cdrme` | bit-identical for a [fixed input order](#input-order) | no effect | none | [README](https://github.com/csilab-ufop/pymocd/blob/master/src/core/algorithms/cdrme/README.md) |
| GDPSO | `gdpso` | bit-identical for a [fixed input order](#input-order) | no effect | none | [README](https://github.com/csilab-ufop/pymocd/blob/master/src/core/algorithms/gdpso/README.md) |
| MMCoMO | `mmcomo`, `mmcomo_fronts` | can differ | single thread | none | [README](https://github.com/csilab-ufop/pymocd/blob/master/src/core/algorithms/mmcomo/README.md) |
| NSGA-III-CCM | `ccm`, `ccm_fronts` | can differ | single thread | none | [README](https://github.com/csilab-ufop/pymocd/blob/master/src/core/algorithms/ccm/README.md) |
| NSGA-III-KRM | `krm`, `krm_fronts` | can differ | single thread | none | [README](https://github.com/csilab-ufop/pymocd/blob/master/src/core/algorithms/krm/README.md) |
| Shi-MOCD | `mocd_q`, `mocd_d` | can differ | single thread | none | [README](https://github.com/csilab-ufop/pymocd/blob/master/src/core/algorithms/mocd/README.md) |
| MOGA-Net | `moga_net`, `moga_net_fronts` | can differ | single thread | none | [README](https://github.com/csilab-ufop/pymocd/blob/master/src/core/algorithms/moganet/README.md) |

### RIMPSO: seeds and threads

RIMPSO is the only detector with a `seed` parameter, accepted by
[`rimpso`][pymocd.rimpso] and [`rimpso_fronts`][pymocd.rimpso_fronts]. Its
docstring states that the same graph and the same parameters, `seed`
included, give the same partition on any number of threads, where "the same
graph" means the same nodes and edges listed in the same order (see
[Input order](#input-order)). The module README
names what enforces this: one independent random stream per (iteration,
particle), every tie broken by the lower index or label, archive candidates
offered in slot order rather than completion order, and particles relabeled to
a shared naming after each flight.

The seed is mixed into every stream (`src/core/algorithms/rimpso/utils/sampling.rs`):

- `seed=0`, the default, contributes nothing to the streams, so it reproduces
  the trajectory RIMPSO followed before the seed became a parameter;
- any other value flies an independent trajectory;
- the seed is an unsigned 64-bit integer, from 0 to $2^{64}-1$; a negative
  value raises `OverflowError`.

Different seeds give independent searches, which may or may not end at the
same partition. On the Les Misérables co-appearance network:

```python exec="on" source="above" session="performance" result="text"
import networkx as nx
import pymocd

lesmis = nx.convert_node_labels_to_integers(nx.les_miserables_graph())

print("repeated call identical:", pymocd.rimpso(lesmis) == pymocd.rimpso(lesmis))
print("default equals seed=0: ", pymocd.rimpso(lesmis) == pymocd.rimpso(lesmis, seed=0))

base = pymocd.rimpso(lesmis, seed=0)
for seed in (1, 2, 3):
    p = pymocd.rimpso(lesmis, seed=seed)
    print(f"seed={seed}: {len(set(p.values()))} communities, ARI vs seed 0 = {pymocd.ari(p, base):.3f}")
```

An ARI of 1 means two seeds ended at the same partition; anything lower means
they did not. To report RIMPSO's spread over runs, vary `seed`, not the call
count: repeated calls with one seed return the same partition.

Thread-count independence can only be shown across processes, because each
process sizes its pool once. The block below runs RIMPSO in two fresh
processes, one with a single thread and one with four, and compares them with
the run in the current process:

```python exec="on" source="above" session="performance" result="text"
import json

FIXED_THREADS = """
import json, sys
import networkx as nx
import pymocd

pymocd.max_cores(int(sys.argv[1]))  # first call in this fresh process
G = nx.convert_node_labels_to_integers(nx.les_miserables_graph())
print(json.dumps(sorted(pymocd.rimpso(G, seed=7).items())))
"""

def rimpso_with_threads(threads):
    out = subprocess.run([sys.executable, "-c", FIXED_THREADS, str(threads)],
                         capture_output=True, text=True, check=True).stdout
    return dict(json.loads(out.strip().splitlines()[-1]))  # max_cores logs first

one, four = rimpso_with_threads(1), rimpso_with_threads(4)
print("1 thread == 4 threads:     ", one == four)
print("fresh process == this one: ", one == pymocd.rimpso(lesmis, seed=7))
```

### CDRME and GDPSO: deterministic without a seed

The CDRME and GDPSO READMEs report both detectors as bit-deterministic and
independent of the thread count, each measured by a test that compares two
runs byte for byte and a test that changes the thread count. Every random draw
comes from a stream identified by a salt (an internal constant, or the
generation number in GDPSO) and a per-task slot index, so no draw depends on
Rayon's scheduling. Neither function takes a seed, so the
public API offers no way to draw a second, independent run: one run per graph
and parameter setting is the whole distribution. The CDRME README reports its
measurements as single runs for this reason.

### Input order

RIMPSO, CDRME and GDPSO are deterministic for a fixed graph object: repeated
calls on the same object return the same partition. A detector call reads the
nodes in the order `G.nodes()` lists them and the edges in the order
`G.edges()` lists them (`src/core/graph/adj.rs`), and these three detectors
number the nodes and fill each node's neighbor list in that order
(`CsrGraph::from_edges` in `src/core/graph/csr.rs`). Two equal
graphs, with the same nodes and the same edges, can therefore give different
partitions when they list their nodes or edges in a different order, for
example when the same network is loaded from another file or built by another
library. In tests with the current package, reordering only the nodes of the
Les Misérables network changed the partition returned by each of the three
detectors; reordering only its edges changed RIMPSO's and GDPSO's, and
reordering only the edges of a planted-partition graph changed CDRME's.

To reproduce a result, keep both orders fixed: add the nodes in sorted order
and then the edges in sorted order. The block below lists the edges of the
Les Misérables graph from above in a shuffled order, then rebuilds both
graphs in sorted order; a `False` on the second line means the order alone
changed RIMPSO's partition:

```python exec="on" source="above" session="performance" result="text"
import random

edges = list(lesmis.edges(data=True))
random.Random(0).shuffle(edges)
shuffled = nx.Graph()
shuffled.add_nodes_from(lesmis.nodes(data=True))  # same node order
shuffled.add_edges_from(edges)  # same edges, another order

def sorted_copy(graph):
    H = nx.Graph()
    H.add_nodes_from(sorted(graph))
    H.add_edges_from(sorted(tuple(sorted(e)) for e in graph.edges()))
    return H

print("equal graphs:              ", nx.utils.graphs_equal(lesmis, shuffled))
print("same partition:            ", pymocd.rimpso(shuffled) == pymocd.rimpso(lesmis))
print("same after sorted rebuild: ",
      pymocd.rimpso(sorted_copy(shuffled)) == pymocd.rimpso(sorted_copy(lesmis)))
```

### Nondeterministic detectors

HP-MOCD, MMCoMO, NSGA-III-CCM, NSGA-III-KRM, Shi-MOCD and MOGA-Net draw from
Rust's thread-local generator (`rand::rng()`), which their READMEs describe as
unseeded or seeded from the operating system, and expose no seed. Seeding
Python's `random` module or NumPy does not affect them. The HP-MOCD README
lists five independent sources of variation, none of which can be switched
off: unseeded generators, Rayon work distribution, front-peeling order, hash
iteration order and float summation order. For the five single-threaded
baselines, the READMEs name the random stream as the source of run-to-run
variation.

Identical output on a small graph does not show that a detector is
deterministic. The READMEs of MMCoMO and NSGA-III-CCM report 30 karate runs
that returned the same score every time (MMCoMO's: the same partition), and
attribute this to the graph rather than to determinism. The MOGA-Net README
reports the same, but in a test with the current package 30 karate runs of
`moga_net` returned two distinct partitions (see
[MOGA-Net › Randomness](../algorithms/moga-net.md#randomness)). The HP-MOCD
README reports that over 1000 runs on karate the detector returned between 3
and 5 communities.

For these detectors, treat one call as one sample. Run each configuration
several times, keep the partitions you report, and give the mean and spread
together with the number of runs; the module READMEs report their reference
measurements this way, typically over 30 runs. For example, on a planted
partition graph with known blocks:

```python exec="on" source="above" session="performance" result="text"
import statistics

G = nx.planted_partition_graph(4, 25, 0.25, 0.1, seed=1)
gt = {v: i for i, block in enumerate(G.graph["partition"]) for v in block}

runs = [pymocd.hpmocd(G) for _ in range(5)]
amis = [pymocd.ami(p, gt) for p in runs]
print("distinct partitions in 5 runs:", len({tuple(sorted(p.items())) for p in runs}))
print(f"AMI over 5 runs: mean {statistics.mean(amis):.3f}, sd {statistics.stdev(amis):.3f}")
```

## What dominates runtime and memory

What the module READMEs state about each detector's cost is summarized here;
each algorithm page's *Complexity* section has the details.

- **MMCoMO** keeps a dense $n \times n$ diffusion kernel and builds a dense
  vote matrix at each exchange, both $O(n^2)$ in memory, and computes the
  kernel with an $O(n^3)$ eigendecomposition. The README identifies this as
  the cap on the graph sizes the detector can handle; the repository's
  synthetic-network benchmark registers MMCoMO with `max_nodes=2000`. See
  [MMCoMO › Complexity](../algorithms/mmcomo.md#complexity), and
  [Measured time and memory](#measured-time-and-memory) below for its cost at
  1,000 and 10,000 nodes.
- **Shi-MOCD**: `mocd_d` runs `rand_networks + 1` full searches, one on the
  graph and one on each random control graph, each with the same `pop_size`
  and `num_gens`. See [Shi-MOCD › Complexity](../algorithms/shi-mocd.md#complexity).
- **CDRME**: the README states that the cost is linear in `pop_size`, and
  reports 6.0 s at `pop_size=100` and 14.3 s at `pop_size=300` on a 50k-node,
  241k-edge stochastic block model with 48 cores. See
  [CDRME › Complexity](../algorithms/cdrme.md#complexity).
- **GDPSO** evaluates modularity in a sparse $O(n + m)$ form rather than the
  reference implementation's dense double loop (README), once per particle per
  generation (`swarm/engine.rs`). See [GDPSO › Complexity](../algorithms/gdpso.md#complexity).
- **RIMPSO** always runs all `num_gens` generations of its `pop_size`
  particles. The README describes incremental objective counts, with a full
  recount only at seeding and after a merge sweep, `n`-sized scratch buffers
  built once per run with one set per worker thread, and a selection step
  that costs one pass over the labels and one over the edges per archive
  member. See
  [RIMPSO › Complexity](../algorithms/rimpso.md#complexity).
- **HP-MOCD** runs objective evaluation, offspring production, the
  non-dominated sort and mutation on the Rayon pool, some of them only above
  size thresholds listed in its README. Objectives supplied as Python
  callables to [`HpMocd`][pymocd.HpMocd] run sequentially. See
  [HP-MOCD › Complexity](../algorithms/hpmocd.md#complexity).
- **NSGA-III-CCM**, **NSGA-III-KRM** and **MOGA-Net**: the module READMEs do
  not state what dominates their cost. See
  [NSGA-III-CCM](../algorithms/ccm.md#complexity),
  [NSGA-III-KRM](../algorithms/krm.md#complexity) and
  [MOGA-Net](../algorithms/moga-net.md#complexity).

## Measured time and memory

The two tables below give the wall time and the peak resident memory (RSS)
of one call of each detector at its default parameters, on LFR graphs of
1,000, 10,000 and 100,000 nodes. They measure cost only, not the quality of
the partitions.

**How it was measured.** The graphs came from
`pymocd.lfr_graph(n, 2.5, 1.5, 0.3, average_degree=10, max_degree=50, min_community=20, max_community=100, seed=1)`;
their mean degree is 9.8 to 10.0, no node has a degree above 50, and they
have 20, 233 and 2,268 ground-truth communities. Each call ran in a fresh
Python process, which built a `networkx.Graph` with the nodes `0` to `n - 1`
in order and the edges in the sorted order `lfr_graph` returns them, and
then called the detector once. The time is `time.perf_counter()` around the
call alone, so it includes reading the networkx graph, building the
detector's internal graph, building Rayon's thread pool (for the parallel
detectors) and building the returned dict. The peak RSS is the process's
`VmHWM` from `/proc/self/status` after the call, with the counter reset just
before the call (by writing `5` to `/proc/self/clear_refs`). It therefore
includes the interpreter, networkx and the graph object, whose resident size
just before the call is the first row of the memory table. A call that took
under 60 s was run three times, each in its own process, and the tables give
the median time and the median peak; a time above 60 s is from a single run.
Each call had a limit of 20 minutes.

The machine was jackal: an AMD Ryzen Threadripper 7960X (24 cores, 48
logical CPUs) under Linux 6.8, with Python 3.10.12, pymocd 6.1.0 from PyPI
and networkx 3.4.2. No call to `max_cores` was made and `RAYON_NUM_THREADS`
was not set, so the four parallel detectors ran on Rayon's default pool of 48
threads (each such process had 49 threads after the call: the main thread
and 48 workers) and the other six on one thread. The machine was shared:
other jobs kept the 1-minute load average between about 11 and 21
whenever no parallel call of this measurement was running, and the
measurement ran at nice level 5. The calls ran one at a time, except the
MMCoMO call at $n$ = 10,000, which ran on one CPU alongside the
single-threaded calls at $n$ = 10,000 and all the calls at $n$ = 100,000.

Wall time of one call, in seconds:

| Call | Threads | $n$ = 1,000 ($m$ = 4,887) | $n$ = 10,000 ($m$ = 49,851) | $n$ = 100,000 ($m$ = 502,196) |
|---|---|---|---|---|
| `rimpso` | 48 | 0.065 | 0.38 | 3.98 |
| `hpmocd` | 48 | 0.27 | 1.44 | 18.9 |
| `cdrme` | 48 | 0.027 | 0.71 | 165 |
| `gdpso` | 48 | 0.10 | 1.16 | 8.40 |
| `mmcomo` | 1 | 55.6 | > 1,200 (limit) | not run |
| `ccm` | 1 | 1.64 | 17.6 | 270 |
| `krm` | 1 | 0.66 | 8.64 | 119 |
| `mocd_q` | 1 | 0.40 | 4.26 | 61.8 |
| `mocd_d` | 1 | 1.55 | 16.5 | 249 |
| `moga_net` | 1 | 0.42 | 4.52 | 70.7 |

Peak RSS of the process during the call, in MiB ($2^{20}$ bytes):

| Call | $n$ = 1,000 | $n$ = 10,000 | $n$ = 100,000 |
|---|---|---|---|
| *(resident before the call)* | 32.9 | 44.8 | 168.6 |
| `rimpso` | 38.8 | 81.0 | 609 |
| `hpmocd` | 51.2 | 158 | 1,105 |
| `cdrme` | 40.9 | 134 | 1,865 |
| `gdpso` | 36.7 | 64.4 | 346 |
| `mmcomo` | 62.6 | ≥ 2,338 when stopped | not run |
| `ccm` | 39.0 | 95.3 | 683 |
| `krm` | 36.3 | 71.3 | 441 |
| `mocd_q` | 36.0 | 61.9 | 344 |
| `mocd_d` | 37.8 | 86.8 | 484 |
| `moga_net` | 41.0 | 117 | 890 |

From $n$ = 10,000 to $n$ = 100,000, ten times the nodes and the edges, the
time of the five single-threaded detectors other than MMCoMO grew by 13.8 to
15.6 times, and the memory a call added above the resident size before it
grew by 7.5 to 10.3 times. `mocd_d` took 3.9 to 4.0 times as long as
`mocd_q` at every size, in line with its four searches at the default
`rand_networks=3`. The time of RIMPSO, HP-MOCD and GDPSO grew by 10.6, 13.1
and 7.3 times. CDRME's grew by 231 times, from 0.71 to 165 s, and the memory
it added by 19 times (see [CDRME › Complexity](../algorithms/cdrme.md#complexity)).
MMCoMO took 55.6 s at $n$ = 1,000, where every other detector took under 2 s,
and did not finish within 20 minutes at $n$ = 10,000. When it was stopped
its process held 2,338 MiB, about three times the 763 MiB of one dense
10,000 × 10,000 matrix of 8-byte floats; it was not run at $n$ = 100,000.

At $n$ = 10,000 and 100,000, the CPU time (user plus system) of the parallel
detectors was 18 to 32 times their wall time; for the single-threaded
detectors the two agreed to within 1%. Between the three runs of a call, the
range of the times (slowest minus fastest) was at most 13% of the median for
the single-threaded detectors and up to 66% for the parallel ones, whose 48
threads compete with the other jobs for the 48 logical CPUs. The widest
range in seconds was HP-MOCD's at $n$ = 100,000, from 15.3 to 24.0 s. The
peak RSS of the three runs differed by less than 10%. These numbers describe
this machine at this load; on another machine, or with another thread count
set through [`max_cores`][pymocd.max_cores], measure again as shown below.

## Timing a detector

Each detector call reads the nodes and edges out of the input graph
(`src/core/graph/adj.rs`) and builds the detector's internal graph before
searching, and builds a Python dict from the result, so a wall-clock
measurement includes both steps.
The first parallel call in a process also builds the Rayon pool. Time several
calls after a warm-up call and report a robust summary such as the median:

```python exec="on" source="above" session="performance" result="text"
import time

def median_seconds(fn, *args, repeats=5, **kwargs):
    fn(*args, **kwargs)  # warm-up
    samples = []
    for _ in range(repeats):
        start = time.perf_counter()
        fn(*args, **kwargs)
        samples.append(time.perf_counter() - start)
    return statistics.median(samples)

G = nx.planted_partition_graph(4, 50, 0.25, 0.02, seed=42)
print(f"graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
for num_gens in (25, 50, 100):
    t = median_seconds(pymocd.rimpso, G, num_gens=num_gens)
    print(f"rimpso, num_gens={num_gens:>3}: {1000 * t:7.1f} ms (median of 5)")
```

The times printed above were measured when this page was built, on the build
machine and with its thread count; they are not benchmark results.

## See also

- API reference: [`max_cores`][pymocd.max_cores], [`rimpso`][pymocd.rimpso],
  [`rimpso_fronts`][pymocd.rimpso_fronts].
- [Choosing a detector](choosing.md): costs and limits of each detector.
- [Evaluating partitions](evaluation.md): the scores used above.
- Algorithm pages: [RIMPSO](../algorithms/rimpso.md),
  [HP-MOCD](../algorithms/hpmocd.md), [CDRME](../algorithms/cdrme.md),
  [GDPSO](../algorithms/gdpso.md), [MMCoMO](../algorithms/mmcomo.md).
