---
description: What pymocd's detectors accept as input - networkx and igraph graphs, integer node ids, direction, parallel edges, weights, self-loops and isolated nodes.
---

# Input graphs

Every detector takes a graph as its first argument and returns a
`dict[int, int]` from node id to community id. This page states exactly what
the graph readers in `src/core/graph/` accept and how they treat direction,
parallel edges, weights, self-loops and isolated nodes. Every statement below
comes from that code or from running the version of pymocd this site documents;
the executed blocks show the behavior on small graphs.

## Summary

| Input | What pymocd does |
|---|---|
| `networkx.Graph` | Read as given. |
| `networkx.DiGraph` | Direction is discarded: each arc becomes an undirected edge, so a reciprocal pair `u→v`, `v→u` counts as two parallel edges. |
| `networkx.MultiGraph`, `networkx.MultiDiGraph` | Each parallel edge is read as a separate edge. |
| Edge attributes, including `weight` | Ignored: every graph is treated as unweighted ([details](#edge-weights)). |
| Self-loops | Dropped. |
| Nodes without edges | Kept, with community `-1`. |
| Node ids | Integers in the signed 32-bit range. Strings, floats and tuples raise an error. |
| `igraph.Graph` | Read with the same rules, keyed by vertex index ([details](#igraph)). Accepted from 6.1.0; 6.0.1 fails on igraph input ([Changelog](../changelog.md)). |

`mmcomo`, `mmcomo_fronts` and `cdrme` are the exceptions on repeated edges:
they deduplicate edges while building their own adjacency, so a repeated edge,
or a reciprocal pair of arcs, counts once ([details](#parallel-edges)).

## How a graph is read

Two functions in `src/core/graph/adj.rs` read every input graph: `get_nodes`
calls `graph.nodes()` and `get_edges` calls `graph.edges()`, keeping the first
two items of each edge. Only when `graph.nodes()` (or `graph.edges()`) cannot
be called do they fall back to igraph's `vs` (or `get_edgelist()`).

The detectors reach these readers in one of two ways, which matters for the
exception you get on invalid input:

Direct readers
:   [`rimpso`][pymocd.rimpso], [`rimpso_fronts`][pymocd.rimpso_fronts],
    [`rimpso_select`][pymocd.rimpso_select], [`mmcomo`][pymocd.mmcomo],
    [`mmcomo_fronts`][pymocd.mmcomo_fronts], [`gdpso`][pymocd.gdpso] and
    [`cdrme`][pymocd.cdrme] propagate the reader's error as an ordinary Python
    exception (`RuntimeError` or `OverflowError`). `rimpso_select` converts
    its `candidates` before it reads the graph, so a candidate dict with
    non-integer keys or labels raises `TypeError`, and one with keys or labels
    outside the 32-bit range raises `OverflowError`.

Through `Graph::from_python`
:   [`hpmocd`][pymocd.hpmocd], [`hpmocd_fronts`][pymocd.hpmocd_fronts],
    [`HpMocd`][pymocd.HpMocd], [`ccm`][pymocd.ccm],
    [`ccm_fronts`][pymocd.ccm_fronts], [`krm`][pymocd.krm],
    [`krm_fronts`][pymocd.krm_fronts], [`moga_net`][pymocd.moga_net],
    [`moga_net_fronts`][pymocd.moga_net_fronts], [`mocd_q`][pymocd.mocd_q] and
    [`mocd_d`][pymocd.mocd_d] build a Rust `Graph` with `Graph::from_python`,
    which unwraps the reader's result. On invalid input the Rust code panics:
    the panic message is printed to standard error and Python receives a
    `pyo3_runtime.PanicException` that wraps the original error.
    `PanicException` derives from `BaseException`, not from `Exception`, so
    `except Exception` does not catch it.

The examples below use the direct readers, whose errors can be caught normally.

## Node ids

Node ids are stored as `NodeId`, a signed 32-bit integer (`i32` in
`src/core/graph/mod.rs`). A node id must therefore be an integer between
$-2^{31}$ and $2^{31}-1$. Within that range the ids do not need to be
contiguous, start at zero or be non-negative, and the result uses the same ids
as keys. Integer-like objects such as `numpy.int64` are accepted and come back
as Python `int` keys:

```python exec="on" source="above" session="input-graphs" result="text"
import networkx as nx
import numpy as np
import pymocd

K = nx.karate_club_graph()

shifted = nx.relabel_nodes(K, lambda v: 1000 * v - 17)   # -17, 983, 1983, ...
p = pymocd.rimpso(shifted)
print("smallest and largest key:", min(p), max(p))

numpy_ids = nx.relabel_nodes(K, lambda v: np.int64(v))
p = pymocd.rimpso(numpy_ids)
print("key type with numpy.int64 node ids:", type(next(iter(p))).__name__)
```

Any other node type (strings, floats, tuples such as the nodes of
`nx.grid_2d_graph`) is rejected. A direct reader raises `RuntimeError`; the
other detectors raise a `PanicException` that wraps it. The message mentions
`<str>` whatever the actual type is, and asks for "positive integers" although
negative ids work, as shown above:

```python exec="on" source="above" session="input-graphs" result="text"
named = nx.relabel_nodes(K, lambda v: f"member-{v}")

try:
    pymocd.rimpso(named)
except RuntimeError as err:
    print("RuntimeError:", err)
```

To run a detector on such a graph, relabel the nodes to integers with
`nx.convert_node_labels_to_integers` and keep the original labels in a node
attribute (`label_attribute`), then map the result back:

```python exec="on" source="above" session="input-graphs" result="text"
H = nx.convert_node_labels_to_integers(named, label_attribute="name")
p = pymocd.rimpso(H)

by_name = {H.nodes[v]["name"]: c for v, c in p.items()}
print(list(by_name.items())[:4])
```

An id outside the 32-bit range raises `OverflowError` (wrapped in a
`PanicException` by the detectors that use `Graph::from_python`) when the node
has at least one edge:

```python exec="on" source="above" session="input-graphs" result="text"
too_big = nx.relabel_nodes(nx.path_graph(3), lambda v: 2**31 + v)

try:
    pymocd.rimpso(too_big)
except OverflowError as err:
    print("OverflowError:", err)
```

!!! warning "Out-of-range ids on isolated nodes are not rejected"
    `get_nodes` reads a node id as a 64-bit integer and then truncates it to
    32 bits, while `get_edges` reads edge endpoints as 32-bit integers
    directly and fails on overflow. A node without edges whose id is outside
    the 32-bit range but inside the signed 64-bit range ($-2^{63}$ to
    $2^{63}-1$) is therefore not rejected: its id is silently truncated and
    can collide with another node's id. An id outside the signed 64-bit range
    fails the 64-bit read and raises the `RuntimeError` shown above for
    strings. Below, the isolated node `2**32 + 1`
    truncates to `1` and does not appear in the result:

    ```python exec="on" source="above" session="input-graphs" result="text"
    P = nx.path_graph(3)
    P.add_node(2**32 + 1)
    print(pymocd.rimpso(P))
    ```

    Relabel graphs whose ids may exceed the 32-bit range with
    `nx.convert_node_labels_to_integers` before running a detector.

## Edges and objectives

The sections below show how edges are counted with
[`rimpso_select`][pymocd.rimpso_select], which scores given partitions on a
graph and returns their `(cut, pair)` points. `cut` is the fraction of edges
that leave their community, so it reveals how many edges the reader saw. The
reference case is the path `0 - 1 - 2` split into `{0, 1}` and `{2}`: one of
its two edges is cut.

```python exec="on" source="above" session="input-graphs" result="text"
split = {0: 0, 1: 0, 2: 1}

def cut(graph):
    """Fraction of edges that leave their community, as RIMPSO counts it."""
    _, points = pymocd.rimpso_select(graph, [split])
    return points[0][0]

print(f"path 0-1-2: cut = {cut(nx.path_graph(3)):.3f}")
```

### Directed graphs

`get_edges` reads every arc `(u, v)` of a `DiGraph` as an undirected edge. The
direction of a single arc is lost, and a reciprocal pair `u→v`, `v→u` becomes
two parallel edges between `u` and `v`. Below, the pair `0→1`, `1→0` counts
twice inside community 0, so only one edge in three is cut; after
`to_undirected()`, which merges the pair into one edge, the cut is one in two
again:

```python exec="on" source="above" session="input-graphs" result="text"
D = nx.DiGraph([(0, 1), (1, 0), (1, 2)])

print(f"DiGraph 0->1, 1->0, 1->2: cut = {cut(D):.3f}")
print(f"D.to_undirected():        cut = {cut(D.to_undirected()):.3f}")
```

To give every connected pair of nodes one edge, convert a directed graph with
`G.to_undirected()` before detection. Except for `mmcomo`, `mmcomo_fronts` and
`cdrme`, a detector given `G.to_directed()` of an undirected graph `G` reads
every edge twice, so its result can differ from the one it returns for `G`.

### Parallel edges

For a `MultiGraph` (or `MultiDiGraph`), `G.edges()` yields one tuple per
parallel edge, and the readers keep each one. Below, the edge `1 - 2` appears
twice, so two edges in three are cut; `nx.Graph(M)` collapses the parallel
edges:

```python exec="on" source="above" session="input-graphs" result="text"
M = nx.MultiGraph([(0, 1), (1, 2), (1, 2)])

print(f"MultiGraph with 1-2 twice: cut = {cut(M):.3f}")
print(f"nx.Graph(M):               cut = {cut(nx.Graph(M)):.3f}")
```

Which detectors keep repeated edges, per the code:

- `Graph::from_python` (`src/core/graph/adj.rs`) keeps parallel edges as they
  arrive, and `CsrGraph::from_edges` (`src/core/graph/csr.rs`) keeps a
  repeated edge; the GDPSO module README calls this "the multigraph reading".
  This covers every detector except the two below.
- `mmcomo` and `mmcomo_fronts` build their adjacency from sets
  (`src/core/algorithms/mmcomo/graph.rs`), and `cdrme` deduplicates each
  neighbor list (`src/core/algorithms/cdrme/topology.rs`). For these, a
  repeated edge or a reciprocal pair of arcs counts once.

Whether a repeated edge changes a particular detector's result depends on its
objectives. To avoid the question, pass a simple undirected `networkx.Graph`.

### Edge weights

`get_edges` keeps only the two endpoints of each edge, so edge attributes,
including `weight`, never reach the detectors: every graph is treated as
unweighted.

```python exec="on" source="above" session="input-graphs" result="text"
W = nx.path_graph(3)
W.edges[1, 2]["weight"] = 10.0
print(f"path 0-1-2 with weight 10 on 1-2: cut = {cut(W):.3f}")

print("karate club edge (0, 1):", K.edges[0, 1])
unit = K.copy()
for u, v in unit.edges:
    unit.edges[u, v]["weight"] = 1
print("same partition with every weight set to 1:", pymocd.rimpso(K) == pymocd.rimpso(unit))
```

`nx.karate_club_graph()` carries a `weight` attribute on its edges; the
detectors ignore it.

The one place where the original graph object reaches Python code is a custom
objective given to [`HpMocd`][pymocd.HpMocd]: each objective is called as
`objective(graph, partition)` with the graph you passed in
(`src/core/algorithms/hpmocd/api.rs`), so an objective you write can read edge
attributes itself. HP-MOCD's search operators still work on the unweighted
graph built by `Graph::from_python`.

### Self-loops

Self-loops are dropped when the graph is read (`Graph::from_python` skips an
edge whose endpoints are equal, and `CsrGraph::from_edges` and MMCoMO's graph
builder do the same):

```python exec="on" source="above" session="input-graphs" result="text"
L = nx.path_graph(3)
L.add_edge(0, 0)
print(f"path 0-1-2 with a self-loop on 0: cut = {cut(L):.3f}")
```

## Isolated nodes

Every node returned by `G.nodes()` appears in the result. A node of degree
zero is not a member of any community and gets community `-1`. Because
self-loops are dropped, this includes a node whose only edge is a self-loop:

```python exec="on" source="above" session="input-graphs" result="text"
G = nx.karate_club_graph()
G.add_node(100)        # no edges
G.add_edge(101, 101)   # only a self-loop

p = pymocd.rimpso(G)
print("node 100 ->", p[100], "| node 101 ->", p[101])
print("entries in the result:", len(p), "of", G.number_of_nodes(), "nodes")
```

The front accessors (`*_fronts`) apply the same rule to every partition they
return. The metric functions treat `-1` as an ordinary community label, so all
isolated nodes count as one community when a partition is scored; see
[Evaluating partitions](evaluation.md).

## igraph

An `igraph.Graph` is read through the fallback path described in
[How a graph is read](#how-a-graph-is-read): node ids are the vertex indices
`0 … n-1` (vertex names are not used), and edges come from `get_edgelist()`.
The rules above apply unchanged: a directed graph's arcs each count as an
edge, multi-edges are kept, the `weight` attribute and other edge attributes
are ignored, self-loops are dropped, and vertices without edges get community
`-1`.

Each call on an `igraph.Graph` prints one line, wrapped in ANSI colour codes,
to the process's standard output when the edge reader falls back to
`get_edgelist()`:
`[ WARNING: src/core/graph/adj.rs:81 ] networkx.Graph() not found, trying igraph.Graph()`.
It is written by Rust's `println!` (the `debug!` macro in
`src/core/macros.rs`), so it bypasses `sys.stdout`: `contextlib.redirect_stdout`
does not capture it. It is informational and the result is unaffected.

Because the keys are vertex indices, the result converts straight to an
igraph membership list:

```python exec="on" source="above" session="input-graphs" result="text"
import igraph as ig

g = ig.Graph.Famous("Zachary")
partition = pymocd.rimpso(g)
membership = [partition[v] for v in range(g.vcount())]
print(len(membership), "vertices,", len(set(membership)), "communities")

twin = nx.Graph()
twin.add_nodes_from(range(g.vcount()))
twin.add_edges_from(g.get_edgelist())
print("same partition as the NetworkX copy:", partition == pymocd.rimpso(twin))
```

!!! warning "igraph and pymocd 6.0.1"
    `igraph.Graph` input works from pymocd 6.1.0. In 6.0.1, `get_nodes`
    called `vs.__iter__()` as a method, which igraph's `VertexSeq` does not
    have, so every detector failed on an `igraph.Graph` with
    `AttributeError: 'VertexSeq' object has no attribute '__iter__'` (or a
    `PanicException` wrapping it). With 6.0.1, convert the graph first with
    `g.to_networkx()`, whose node ids are the vertex indices.

## See also

- [Getting started](../getting-started.md) for a first detection on the
  karate club graph.
- [Core concepts](concepts.md) for what the detectors optimize.
- The [API reference](../api/index.md) for every detector's signature.
