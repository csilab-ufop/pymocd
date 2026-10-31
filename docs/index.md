---
title: pymocd
description: Multi-objective evolutionary community detection for NetworkX and igraph graphs, with a parallel Rust core.
hide:
  - navigation
  - toc
---

<h1 class="pm-sr-only">pymocd</h1>

<div class="pm-hero" markdown>

![pymocd: multi-objective community detection](assets/logo.svg#only-light)
![pymocd: multi-objective community detection](assets/logo-dark.svg#only-dark)

<p class="pm-tagline">Multi-objective evolutionary community detection for NetworkX and igraph graphs, with a parallel Rust core.</p>

[Get started](getting-started.md){ .md-button .md-button--primary }
[Examples](examples/index.md){ .md-button }
[API reference](api/index.md){ .md-button }

</div>

<div class="pm-start" markdown>

<div markdown>

## Install

```bash
pip install pymocd networkx
```

pymocd itself has no runtime dependencies; NetworkX (or python-igraph)
supplies the graphs. Prebuilt wheels cover CPython 3.10–3.13 on Linux, macOS
and Windows; see [Installation](guide/installation.md) for the platforms of
each release.


</div>

<div markdown>

## Detect communities

```python exec="on" source="above" session="home" result="text"
import networkx as nx
import pymocd

communities = pymocd.rimpso(nx.karate_club_graph())   # {node: community}
print(len(set(communities.values())), "communities")
```

</div>

</div>

<div class="grid cards pm-features" markdown>

-   :material-graph-outline:{ .lg .middle } **Nine algorithms, one interface**

    ---

    Two methods developed in this library, RIMPSO and HP-MOCD, and seven
    re-implemented baselines from the literature. Every detector takes a graph
    and returns a plain `dict[node, community]`.

    [:octicons-arrow-right-24: Algorithms](algorithms/index.md)

-   :material-language-rust:{ .lg .middle } **Rust core, Rayon parallelism**

    ---

    The evolutionary search is written in Rust and exposed to Python through
    PyO3. RIMPSO, HP-MOCD, CDRME and GDPSO run in parallel on a shared Rayon
    thread pool, sized with [`max_cores`][pymocd.max_cores].

    [:octicons-arrow-right-24: Performance](guide/performance.md)

-   :material-file-document-check-outline:{ .lg .middle } **Documented baselines**

    ---

    Each re-implementation records the paper it follows, its parameters, and
    every deliberate divergence from the original method.

    [:octicons-arrow-right-24: Re-implemented baselines](algorithms/index.md#re-implemented-baselines)

-   :material-ruler-square-compass:{ .lg .middle } **Native evaluation metrics**

    ---

    NMI, AMI, ARI and pairwise F1 against ground truth, computed in Rust over
    the nodes both partitions share.

    [:octicons-arrow-right-24: Evaluating partitions](guide/evaluation.md)

</div>

<div class="pm-cite" markdown>

**Using pymocd in your research?** Please cite the HP-MOCD paper: Santos et al.,
*A high-performance evolutionary multiobjective community detection algorithm*,
Social Network Analysis and Mining 15, 110 (2025).

[How to cite](citing.md){ .md-button .md-button--primary }
[Publications](publications.md){ .md-button }

</div>
