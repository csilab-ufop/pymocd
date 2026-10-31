---
description: How to get and use the Pareto fronts behind pymocd's detectors — which functions expose a front, their exact return shapes, selecting a member by modularity, ground truth, community count or RIMPSO's rule, and reproducing papers' best-in-front protocols.
---

# Working with Pareto fronts

A multi-objective search ends with a set of trade-off partitions, a Pareto
front, and each multi-objective detector in pymocd returns one member of it, chosen by the
selection rule published with its algorithm (see
[Core concepts](concepts.md#selecting-one-partition-from-a-front)). The front
accessors expose the whole candidate set instead, so you can make that choice
yourself: against ground truth, for a known number of communities, with your
own quality function, or with RIMPSO's label-free rule applied to another
detector's front. They are also what the best-in-front evaluation protocols of
several re-implemented papers need. This page lists which detectors have a
front accessor, gives the exact return shapes, and shows each way of selecting
a member.

<a id="the-_fronts-functions"></a>

## Which detectors expose a front

Six of the ten detector entry points pair with a front accessor:

| Detector | Returns one partition | Front accessor | Accessor returns |
|---|---|---|---|
| [RIMPSO](../algorithms/rimpso.md) | [`rimpso`][pymocd.rimpso] | [`rimpso_fronts`][pymocd.rimpso_fronts] | `(partitions, points, selected)` |
| [HP-MOCD](../algorithms/hpmocd.md) | [`hpmocd`][pymocd.hpmocd], [`HpMocd.run`][pymocd.HpMocd.run] | [`hpmocd_fronts`][pymocd.hpmocd_fronts], [`HpMocd.generate_pareto_front`][pymocd.HpMocd.generate_pareto_front] | `list[dict]`; `[(partition, objectives), ...]` |
| [MMCoMO](../algorithms/mmcomo.md) | [`mmcomo`][pymocd.mmcomo] | [`mmcomo_fronts`][pymocd.mmcomo_fronts] | `list[dict]` |
| [NSGA-III-CCM](../algorithms/ccm.md) | [`ccm`][pymocd.ccm] | [`ccm_fronts`][pymocd.ccm_fronts] | `list[dict]` |
| [NSGA-III-KRM](../algorithms/krm.md) | [`krm`][pymocd.krm] | [`krm_fronts`][pymocd.krm_fronts] | `list[dict]` |
| [MOGA-Net](../algorithms/moga-net.md) | [`moga_net`][pymocd.moga_net] | [`moga_net_fronts`][pymocd.moga_net_fronts] | `list[dict]` |
| [Shi-MOCD](../algorithms/shi-mocd.md) | [`mocd_q`][pymocd.mocd_q], [`mocd_d`][pymocd.mocd_d] | none | — |
| [GDPSO](../algorithms/gdpso.md) | [`gdpso`][pymocd.gdpso] | none (single objective) | — |
| [CDRME](../algorithms/cdrme.md) | [`cdrme`][pymocd.cdrme] | none (single objective) | — |

The gaps have two different causes:

- **Single objective.** [`gdpso`][pymocd.gdpso] maximizes Newman–Girvan
  modularity and [`cdrme`][pymocd.cdrme] maximizes the linkage sum of its
  paper's Eq. (12). Neither has a Pareto front, so there is no `gdpso_fronts`
  or `cdrme_fronts`.
- **Not exposed.** [`mocd_q`][pymocd.mocd_q] and [`mocd_d`][pymocd.mocd_d] are
  multi-objective, but the Shi-MOCD module README records that the Rust class
  that produces their front is not registered with the Python module, and there
  is no `mocd_fronts`. HP-MOCD optimizes the same $(\mathrm{intra}, \mathrm{inter})$
  pair and exposes its front, so it is the route to that objective space.

Each `*_fronts` function takes the same keyword arguments as its detector, with
the same defaults, and runs the same search; only the return value differs. Two
details matter when you combine a detector with its accessor:

- **Every call is a new search.** Of the six, only RIMPSO is deterministic; for it,
  `rimpso(G)` is exactly `partitions[selected]` from `rimpso_fronts(G)`. The
  other detectors have no seed parameter, so calling `moga_net(G)` and then
  `moga_net_fronts(G)` runs two independent searches, and the first result need
  not be a member of the second. To get a front and a detector-style choice from
  the same search, call the accessor once and apply the rule to its output
  ([below](#by-maximum-modularity)).
- **HP-MOCD's budget.** [`hpmocd`][pymocd.hpmocd] and
  [`hpmocd_fronts`][pymocd.hpmocd_fronts] take only the graph and run with the
  module's fixed defaults (`cross_rate=0.7`, `mut_rate=0.5`), which differ from
  the crossover and mutation probabilities of the paper's experiments (0.8 and
  0.2; see [HP-MOCD](../algorithms/hpmocd.md)). To change the budget, or to read
  each member's objective values, use the [`HpMocd`][pymocd.HpMocd] class
  ([below](#objective-values-with-hpmocd)).

## Return shapes

### `rimpso_fronts`

[`rimpso_fronts`][pymocd.rimpso_fronts] returns a tuple
`(partitions, points, selected)`:

- `partitions` is a `list[dict[int, int]]`, RIMPSO's archive;
- `points` is a list of `(cut, pair)` tuples, one per partition, in the same
  order (both objectives are defined in
  [Core concepts](concepts.md#the-cpm-split-rimpso));
- `selected` is the index of the member that [`rimpso`][pymocd.rimpso] returns.

```python exec="on" source="above" session="fronts" result="text"
import networkx as nx
import pymocd

G = nx.karate_club_graph()
partitions, points, selected = pymocd.rimpso_fronts(G)

print(f"{len(partitions)} partitions, {len(points)} points, selected = {selected}")
print("rimpso(G) is the selected member:", pymocd.rimpso(G) == partitions[selected])

# Walk the resolution profile from coarse (low cut) to fine (high cut).
order = sorted(range(len(points)), key=lambda i: points[i][0])
for i in order[::10]:
    cut, pair = points[i]
    print(f"member {i:>2}: cut = {cut:.3f}  pair = {pair:.3f}  communities = {len(set(partitions[i].values()))}")
```

### The `list[dict]` accessors

[`hpmocd_fronts`][pymocd.hpmocd_fronts], [`mmcomo_fronts`][pymocd.mmcomo_fronts],
[`ccm_fronts`][pymocd.ccm_fronts], [`krm_fronts`][pymocd.krm_fronts] and
[`moga_net_fronts`][pymocd.moga_net_fronts] return a plain list of partitions,
each a `dict[int, int]`, with no objective values. For the NSGA-II and
NSGA-III detectors the list is the rank-1 front of the final population; for
MMCoMO it is the merged rank-1 front of its two populations. Front sizes depend
on the detector and its budget, and for these non-deterministic detectors they
also vary from call to call.

```python exec="on" source="above" session="fronts" result="text"
for name in ("hpmocd_fronts", "mmcomo_fronts", "ccm_fronts", "krm_fronts", "moga_net_fronts"):
    front = getattr(pymocd, name)(G)
    print(f"{name:<16} {type(front).__name__} of {len(front):>3} {type(front[0]).__name__}s")
```

### `HpMocd.generate_pareto_front`

[`HpMocd.generate_pareto_front`][pymocd.HpMocd.generate_pareto_front] returns
`[(partition, objectives), ...]`, where `objectives` is a list of floats in the
order of the configured objectives: `[intra, inter]` for the built-in pair. See
[Objective values with HpMocd](#objective-values-with-hpmocd).

### Isolated nodes and repeated members

Every front member labels isolated nodes `-1`, exactly as the detectors do, so
count communities without that label. A front can also list the same partition
more than once; the last line of the output below compares a MOGA-Net front's
member count with its number of distinct partitions. The helpers defined here
are used in the rest of this page.

```python exec="on" source="above" session="fronts" result="text"
def communities(partition):
    """dict[node, community] -> list of node sets, the form NetworkX expects."""
    groups = {}
    for node, c in partition.items():
        groups.setdefault(c, set()).add(node)
    return list(groups.values())

def n_communities(partition):
    """Number of communities, not counting the -1 label of isolated nodes."""
    return len({c for c in partition.values() if c != -1})

def modularity(G, partition):
    """Unweighted modularity, the Q that pymocd's max-Q rules use."""
    return nx.community.modularity(G, communities(partition), weight=None)

def unique(front):
    """Drop repeated partitions: same node sets, whatever the labels."""
    seen, kept = set(), []
    for p in front:
        key = frozenset(map(frozenset, communities(p)))
        if key not in seen:
            seen.add(key)
            kept.append(p)
    return kept

H = G.copy()
H.add_node(34)                                   # one isolated node
front = pymocd.moga_net_fronts(H)
print("label of the isolated node in every member:", {p[34] for p in front})
print("members:", len(front), "| distinct partitions:", len(unique(front)))
```

<a id="picking-a-member-yourself"></a>

## Selecting a member yourself

The first three examples below use one MOGA-Net front of the karate club; the
last uses an HP-MOCD front. The same code works on any `list[dict]` front, and
on the `partitions` of `rimpso_fronts`.

<a id="max-modularity"></a>

### By maximum modularity

This is the rule [`moga_net`][pymocd.moga_net] and most other detectors apply.
Pass `weight=None` to NetworkX: pymocd's detectors ignore edge weights, and
NetworkX's karate graph carries a `weight` attribute that would otherwise change
$Q$.

```python exec="on" source="above" session="fronts" result="text"
front = pymocd.moga_net_fronts(G)
best = max(front, key=lambda p: modularity(G, p))
print(f"{len(front)} members; max-Q member: {n_communities(best)} communities, Q = {modularity(G, best):.4f}")
```

<a id="against-ground-truth"></a>

### By ground truth

Score every member against a reference partition with
[`gt_metrics`][pymocd.gt_metrics], which returns `(NMI, AMI, ARI, F1)`, or with
[`nmi`][pymocd.nmi], [`ami`][pymocd.ami], [`ari`][pymocd.ari] and
[`f1`][pymocd.f1] individually (see [Evaluating partitions](evaluation.md)).
Here the reference is the two factions of the karate club, from the `club`
node attribute.

```python exec="on" source="above" session="fronts" result="text"
gt = {v: int(G.nodes[v]["club"] != "Mr. Hi") for v in G}

best = max(front, key=lambda p: pymocd.ari(p, gt))
nmi, ami, ari, f1 = pymocd.gt_metrics(best, gt)
print(f"best-ARI member: {n_communities(best)} communities")
print(f"NMI = {nmi:.3f}  AMI = {ami:.3f}  ARI = {ari:.3f}  F1 = {f1:.3f}")
```

This member is the front oracle: it reads the labels, so it measures what the
search found, not what a detector would return.

!!! note "The ground truth may not be on the front"
    Even the best front member can fall short of a perfect match. The
    ground-truth partition may be *dominated* under the search objectives and
    never survive to the final front; the front bounds what any selection rule
    can recover. [Core concepts](concepts.md#the-ground-truth-may-not-be-on-the-front)
    shows how to test this for RIMPSO's objectives.

<a id="target-community-count"></a>

### By a target community count

When the number of communities is known, keep the distinct members closest to
it and break ties with another criterion, here modularity:

```python exec="on" source="above" session="fronts" result="text"
target = 2
distance = min(abs(n_communities(p) - target) for p in front)
closest = [p for p in unique(front) if abs(n_communities(p) - target) == distance]
best = max(closest, key=lambda p: modularity(G, p))
print(f"{len(closest)} distinct member(s) with {n_communities(best)} communities; best Q = {modularity(G, best):.4f}")
```

### With RIMPSO's selection rule

[`rimpso_select`][pymocd.rimpso_select] runs RIMPSO's label-free rule (the
Schwarz-penalized degree-corrected block model of
[Core concepts](concepts.md#rimpso-a-degree-corrected-assortative-block-model))
over candidate partitions produced anywhere. The RIMPSO module README describes
it as the control that separates the search's contribution from the
selector's. It takes the graph and a list of `dict[int, int]` candidates and
returns a two-element list `[selected_index, points]`, where `points` holds
each candidate's `[cut, pair]` as a two-element list; unpack it as
`pick, points = pymocd.rimpso_select(G, candidates)`.

```python exec="on" source="above" session="fronts" result="text"
hp_front = unique(pymocd.hpmocd_fronts(G))
pick, hp_points = pymocd.rimpso_select(G, hp_front)

print(f"{len(hp_front)} distinct HP-MOCD partitions")
print(f"RIMPSO's rule picks member {pick}: {n_communities(hp_front[pick])} communities, "
      f"(cut, pair) = ({hp_points[pick][0]:.4f}, {hp_points[pick][1]:.4f})")
best_q = max(hp_front, key=lambda p: modularity(G, p))
print(f"HP-MOCD's own rule (max Q) picks {n_communities(best_q)} communities")
```

Four details of `rimpso_select` follow from its source:

- It treats all nodes that share a label as one community, including the `-1`
  of isolated nodes, so on a graph with two or more isolated nodes the `pair` it reports
  for a pymocd partition can be larger than `rimpso_fronts` reports for the same
  partition. The selection score skips isolated nodes, so the chosen index is
  not affected.
- If every candidate is rejected as not assortative, it returns index `0`
  without signalling it.
- An empty candidate list returns `[0, []]`.
- Each candidate should label every node of the graph. A node missing from a
  candidate is not treated as a singleton: it keeps a placeholder label, its
  position in the graph's node list, which can coincide with a community id
  and silently merge it into that community.

## Objective values with HpMocd

The [`HpMocd`][pymocd.HpMocd] class runs HP-MOCD with a configurable budget
(`pop_size`, `num_gens`, `cross_rate`, `mut_rate`) and its
[`generate_pareto_front`][pymocd.HpMocd.generate_pareto_front] method returns
each member with its objective values. With the built-in objectives these are
$(\mathrm{intra}, \mathrm{inter})$, both minimized, and
$Q = 1 - \mathrm{intra} - \mathrm{inter}$. Like [`HpMocd.run`][pymocd.HpMocd.run],
each call runs a new search.

```python exec="on" source="above" session="fronts" result="text"
model = pymocd.HpMocd(G, pop_size=100, num_gens=200)
members = model.generate_pareto_front()          # [(partition, [intra, inter]), ...]

partition, (intra, inter) = max(members, key=lambda pm: 1 - sum(pm[1]))
print(f"{len(members)} members")
print(f"max-Q member: intra = {intra:.4f}, inter = {inter:.4f}, {n_communities(partition)} communities")
print(f"1 - intra - inter = {1 - intra - inter:.6f}")
print(f"NetworkX Q        = {modularity(G, partition):.6f}")
```

`HpMocd` also accepts `objectives`, a list of Python callables
`(graph, partition) -> float`, all minimized, that replace the built-in pair.
The HP-MOCD module README warns that `run` then still returns the member that
maximizes $M - \sum_i f_i$ for $M$ objectives, which orders members by
modularity only under the built-in pair (there it equals $1 + Q$); with custom
objectives, select from `generate_pareto_front` yourself. The module README
also notes that a single partition dict is cleared and refilled for each
individual evaluated, so an objective must not keep a reference to it between
calls. Here the objectives are RIMPSO's $(\mathrm{cut}, \mathrm{pair})$,
searched by NSGA-II and reduced by maximum modularity:

```python exec="on" source="above" session="fronts" result="text"
def cut(graph, partition):
    internal = sum(partition[u] == partition[v] for u, v in graph.edges())
    return 1 - internal / graph.number_of_edges()

def pair(graph, partition):
    n = graph.number_of_nodes()
    sizes = {}
    for c in partition.values():
        sizes[c] = sizes.get(c, 0) + 1
    return sum(s * (s - 1) / 2 for s in sizes.values()) / (n * (n - 1) / 2)

custom = pymocd.HpMocd(G, pop_size=50, num_gens=50, objectives=[cut, pair])
members = custom.generate_pareto_front()         # [(partition, [cut, pair]), ...]
partition, values = max(members, key=lambda pm: modularity(G, pm[0]))
print(f"{len(members)} members; max-Q member has {n_communities(partition)} communities, "
      f"(cut, pair) = ({values[0]:.4f}, {values[1]:.4f})")
```

[`HpMocd.set_on_generation`][pymocd.HpMocd.set_on_generation] registers a
callback `(generation, num_gens, front_size)` that runs once per generation,
which lets you follow the size of the front during the search. The gallery
example [custom objectives](../examples/objectives/plot_custom_objectives.md)
and the [HP-MOCD page](../algorithms/hpmocd.md) cover the class further.

## Reproducing best-in-front protocols

The re-implemented baselines that have a front accessor return the
max-modularity member of their front, but several of their papers report the
best member of the front against ground truth. Reproducing those numbers needs
the whole front, which is why the accessors exist:

| Paper | Accessor | Protocol, per the module README or docstring |
|---|---|---|
| [Pizzuti (2012)](https://doi.org/10.1109/TEVC.2011.2161090), MOGA-Net | [`moga_net_fronts`][pymocd.moga_net_fronts] | Table 1: the best-NMI member of the front. Sec. VI-C fixes `r = 2`, the default. |
| [Shaik, Ravi and Deb (2021)](https://doi.org/10.1007/s42979-020-00382-x), NSGA-III-CCM and NSGA-III-KRM | [`ccm_fronts`][pymocd.ccm_fronts], [`krm_fronts`][pymocd.krm_fronts] | Tables 1–2: the best-NMI and the best-modularity members of the front. |
| [Zhang et al. (2023)](https://doi.org/10.1109/MCI.2023.3277773), MMCoMO | [`mmcomo_fronts`][pymocd.mmcomo_fronts] | Table IV: the best-NMI member; Table III uses maximum modularity, which is what `mmcomo` returns. |
| [Santos et al. (2025)](https://doi.org/10.1007/s13278-025-01519-7), HP-MOCD | [`hpmocd_fronts`][pymocd.hpmocd_fronts] | The docstring recommends best-in-front (selector-free) comparison when other detectors are scored by their front oracle. |

```python exec="on" source="above" session="fronts" result="text"
front = pymocd.moga_net_fronts(G, r=2.0)    # Pizzuti (2012), Sec. VI-C: r = 2
print("MOGA-Net, best NMI over the front:", round(max(pymocd.nmi(p, gt) for p in front), 4))

front = pymocd.ccm_fronts(G)
print("CCM, best NMI over the front:     ", round(max(pymocd.nmi(p, gt) for p in front), 4))
print("CCM, best Q over the front:       ", round(max(modularity(G, p) for p in front), 4))
```

Each call is one independent run of a non-deterministic search, so repeat it
as many times as the paper's protocol requires and aggregate the per-run values.
When one detector is scored best-in-front, score the others the same way: the
`hpmocd_fronts` docstring notes that comparing HP-MOCD's single selected
partition against another detector's front oracle silently handicaps HP-MOCD.
For RIMPSO, use the `partitions` of `rimpso_fronts`.

## See also

- [Core concepts](concepts.md): Pareto dominance, the objective pairs and the
  selection rules.
- [Evaluating partitions](evaluation.md): the metrics used to score members.
- Gallery: [the resolution profile](../examples/fronts/plot_resolution_profile.md)
  and [choosing from a front](../examples/fronts/plot_choosing_from_front.md),
  which plot the trade-off a front spans.
- API: [Pareto fronts](../api/fronts.md), the full signatures of the accessors.
