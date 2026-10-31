---
description: A first pymocd session - detect communities in Zachary's karate club with RIMPSO, score them against the known split, look at the Pareto front behind the result and run a second detector.
---

# Getting started

This page is a first session with pymocd. It detects communities in Zachary's
karate club network with [`rimpso`][pymocd.rimpso], inspects and scores the
partition, looks at the Pareto front the partition was selected from, and runs
other detectors on the same graph. Every code block on this page is executed
when the documentation is built, so the printed output comes from the pymocd
build the site was generated with.

## Installation

```bash
pip install pymocd networkx
```

pymocd declares no runtime dependencies, so networkx, which builds the graphs
used below, is installed separately. [Installation](guide/installation.md)
lists the supported Python versions and platforms and explains how to build
from source.

## First detection

[`rimpso`][pymocd.rimpso] is the library's own particle-swarm detector, and
the project README lists it as the recommended default. Pass it a networkx
graph whose node ids are integers:

```python exec="on" source="above" session="getting-started" result="text"
import networkx as nx
import pymocd

G = nx.karate_club_graph()
partition = pymocd.rimpso(G)

print(type(partition).__name__, "with", len(partition), "entries")
print("node 0 is in community", partition[0])
```

The result is a plain `dict` that maps every node of `G` to an integer
community id. A node without edges gets community `-1`; the karate club graph
has none. [Input graphs](guide/input-graphs.md) states which graph types, node
ids and edge attributes the detectors accept.

`rimpso` is deterministic: its docstring states that the same graph and the
same parameters, `seed` included, give the same partition on any number of
threads. "The same graph" means the same nodes and edges listed in the same
order; listing them in another order can change the partition (see
[Input order](guide/performance.md#input-order)).

## Inspecting the partition

Count the communities and list their members:

```python exec="on" source="above" session="getting-started" result="text"
from collections import Counter

sizes = Counter(partition.values())
print("number of communities:", len(sizes))
for community, size in sorted(sizes.items()):
    members = sorted(v for v, c in partition.items() if c == community)
    print(f"community {community}: {size:2d} nodes {members}")
```

networkx represents a partition as a list of node sets. The conversion is one
line, and it lets you use networkx's own tools on the result, for example its
modularity function:

```python exec="on" source="above" session="getting-started" result="text"
communities = [
    {v for v, c in partition.items() if c == k} for k in sorted(sizes) if k != -1
]
print(f"modularity: {nx.community.modularity(G, communities, weight=None):.3f}")
```

`weight=None` makes networkx ignore the `weight` edge attribute that
`karate_club_graph` carries. pymocd's detectors ignore edge weights too (see
[Input graphs](guide/input-graphs.md#edge-weights)), so the two now score the
same unweighted graph.

## Evaluating results

The karate club graph stores a `club` node attribute, `"Mr. Hi"` or
`"Officer"`; the networkx documentation describes it as the club each member
belonged to after the split. Use it as ground truth:

```python exec="on" source="above" session="getting-started" result="text"
truth = {v: 0 if G.nodes[v]["club"] == "Mr. Hi" else 1 for v in G}

nmi, ami, ari, f1 = pymocd.gt_metrics(partition, truth)
print(f"NMI={nmi:.3f}  AMI={ami:.3f}  ARI={ari:.3f}  F1={f1:.3f}")
```

[`gt_metrics`][pymocd.gt_metrics] returns the normalized mutual information,
the adjusted mutual information, the adjusted Rand index and the pairwise F1
score in one call, computed over the nodes the two dicts share. Each score is
also available on its own: [`nmi`][pymocd.nmi], [`ami`][pymocd.ami],
[`ari`][pymocd.ari] and [`f1`][pymocd.f1]. All four equal 1 only when the two
partitions are identical up to renaming the communities, so compare the number
of communities printed above with the two clubs of the ground truth when you
read them. [Evaluating partitions](guide/evaluation.md) covers the metrics in
detail.

!!! note "Two conventions for the karate club split"
    NetworkX's `club` attribute puts node 8 with Mr. Hi. As a comment in
    `src/core/graph/karate.rs` notes, much of the literature puts node 8 with
    the Officers instead, so a karate club score in a paper may not match one
    computed against this attribute.

## Inspecting Pareto fronts

`rimpso` does not search for a single partition. Its swarm builds a Pareto
front of partitions that trade off two objectives, and a label-free rule then
selects one member. [`rimpso_fronts`][pymocd.rimpso_fronts] runs the same
search with the same defaults and returns the whole front together with the
index the rule selects:

```python exec="on" source="above" session="getting-started" result="text"
partitions, points, selected = pymocd.rimpso_fronts(G)

print("partitions on the front:", len(partitions))
print("selected index:", selected)
print("same partition as rimpso(G):", partitions[selected] == partition)
```

Each entry of `points` is the `(cut, pair)` objective point of the matching
partition: `cut` is the fraction of edges that leave their community and
`pair` is the fraction of node pairs that share one. The partition that puts
every node in one community sits at `cut = 0`, `pair = 1`; splitting a
community never lowers `cut` and never raises `pair`. A sample of the front,
ordered by `cut`:

```python exec="on" source="above" session="getting-started" result="text"
def row(i):
    cut, pair = points[i]
    k = len(set(partitions[i].values()))
    return f"{cut:6.3f} {pair:6.3f} {k:6d} {pymocd.ari(partitions[i], truth):8.3f}"

order = sorted(range(len(partitions)), key=lambda i: points[i][0])
print("   cut   pair      k  ARI(truth)")
for i in order[:: max(1, len(order) // 8)]:
    print(row(i))
print("selected:")
print(row(selected))

best = max(range(len(partitions)), key=lambda i: pymocd.ari(partitions[i], truth))
print("highest ARI on the front (uses the ground truth):")
print(row(best))
```

The last line picks a member with the ground truth, which a real analysis does
not have; it shows how far the label-free choice is from the member of the
front that agrees most with the known split. [Core concepts](guide/concepts.md) explains fronts and
selection, and [Working with Pareto fronts](guide/fronts.md) shows the other
front accessors and how to apply the selection rule to partitions of your own.

## Trying another detector

Every detector has the same call shape: a graph in, a `dict` out. HP-MOCD, the
library's other own method, optimizes a decomposition of modularity with
NSGA-II:

```python exec="on" source="above" session="getting-started" result="text"
hp = pymocd.hpmocd(G)

print("communities:", len(set(hp.values())))
print(f"ARI against the ground truth: {pymocd.ari(hp, truth):.3f}")
print(f"ARI against the RIMPSO partition: {pymocd.ari(hp, partition):.3f}")
```

Unlike `rimpso`, [`hpmocd`][pymocd.hpmocd] is not deterministic, so its output
can differ from one run (and one build of this page) to the next. It takes
only the graph; the [`HpMocd`][pymocd.HpMocd] class exposes its budget and
accepts custom objective functions.

The other detectors are re-implementations of published methods, kept as
baselines ([Algorithms](algorithms/index.md) describes each one). They run the
same way:

```python exec="on" source="above" session="getting-started" result="text"
baselines = ["cdrme", "mmcomo", "ccm", "krm", "gdpso", "mocd_q", "mocd_d", "moga_net"]

print(f"{'detector':10s} {'k':>3s}  ARI(truth)")
for name in baselines:
    found = getattr(pymocd, name)(G)
    k = len(set(found.values()))
    print(f"{name:10s} {k:3d}  {pymocd.ari(found, truth):.3f}")
```

<a id="tuning"></a>
Every detector except `hpmocd` takes its own parameters as keyword arguments,
for example `pymocd.rimpso(G, num_gens=50)`; the algorithm pages list every
parameter with its default.

<a id="threads"></a>
RIMPSO, HP-MOCD, CDRME and GDPSO run on a shared Rayon thread pool. To cap its
size, call [`pymocd.max_cores(n)`][pymocd.max_cores] once, before the first
detection: once a parallel detector has run, the pool already exists and the
call raises an error (see
[Performance and reproducibility](guide/performance.md#sizing-the-pool-with-max_cores)).

## Next steps

- [Installation](guide/installation.md): supported platforms, building from
  source, troubleshooting.
- [Input graphs](guide/input-graphs.md): node ids, directed graphs, parallel
  edges, weights, self-loops, isolated nodes and igraph.
- [Core concepts](guide/concepts.md): multi-objective detection, Pareto fronts
  and selecting one partition from a front.
- [Choosing a detector](guide/choosing.md): which detector fits which
  question.
- [Working with Pareto fronts](guide/fronts.md) and
  [Evaluating partitions](guide/evaluation.md).
- [Performance and reproducibility](guide/performance.md): thread count,
  determinism and seeds.
- [Algorithms](algorithms/index.md): one page per detector, with its paper,
  objectives and parameters.
- [Examples](examples/basics/plot_karate_club.md): the karate club example
  with figures, and the rest of the gallery.
- [API reference](api/index.md): every public function and class.
