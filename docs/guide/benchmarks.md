---
description: Generate LFR benchmark graphs with pymocd.lfr_graph, which takes NetworkX's parameters and defaults but wires every node to its drawn degree and the requested mixing, and see how it differs from networkx.LFR_benchmark_graph.
---

# Benchmark graphs

A detector is judged on graphs whose communities are known in advance. The LFR
benchmark of Lancichinetti, Fortunato and Radicchi
([Phys. Rev. E 78, 046110, 2008](https://doi.org/10.1103/PhysRevE.78.046110))
plants communities in a graph whose degrees and community sizes both follow
power laws. One parameter, the mixing $\mu$, sets the fraction of each node's
edges that leave its community. [`lfr_graph`][pymocd.lfr_graph] generates LFR
graphs with crisp communities: undirected, unweighted, and without overlap. It
takes the parameters, defaults and validation of NetworkX's
`LFR_benchmark_graph`. Unlike NetworkX, it wires the edges so that every node
ends with exactly the degree it drew and with the requested share of its edges
inside its community, and it draws degrees whose expected mean is
`average_degree` exactly. When a draw of degrees and community sizes cannot be
wired that way, it draws again instead of failing. NetworkX's generator
overshoots both the degrees and the mixing
([networkx#6809](https://github.com/networkx/networkx/issues/6809)), and
[the comparison below](#how-it-differs-from-networkx) runs both generators side
by side.

## Generating a graph

`lfr_graph` returns `(edges, ground_truth)`. `edges` lists every undirected edge
once as `(u, v)` with `u < v`, in sorted order, over the nodes `0` to `n - 1`.
`ground_truth` maps every node to its community, numbered from `0`. The block
below generates a graph and measures what it got:

```python exec="on" source="above" session="benchmarks" result="text"
import networkx as nx
import pymocd

params = dict(n=1000, tau1=2.5, tau2=1.5, mu=0.3, average_degree=20,
              max_degree=50, min_community=20, max_community=100)
edges, truth = pymocd.lfr_graph(**params, seed=1)


def measure(n, edges, community):
    """Degrees, and each node's count of edges leaving its community."""
    degree, outside = [0] * n, [0] * n
    for u, v in edges:
        degree[u] += 1
        degree[v] += 1
        if community[u] != community[v]:
            outside[u] += 1
            outside[v] += 1
    return degree, outside


degree, outside = measure(params["n"], edges, truth)
sizes = {}
for c in truth.values():
    sizes[c] = sizes.get(c, 0) + 1
node_mixing = sum(o / d for o, d in zip(outside, degree)) / params["n"]
edge_mixing = sum(truth[u] != truth[v] for u, v in edges) / len(edges)
print(f"{len(sizes)} communities of {min(sizes.values())} to {max(sizes.values())} nodes, "
      f"{len(edges)} edges")
print(f"mixing: requested {params['mu']}, realized {node_mixing:.4f} averaged over nodes, "
      f"{edge_mixing:.4f} as the fraction of edges between communities")
print(f"degree: requested mean {params['average_degree']}, realized mean "
      f"{sum(degree) / params['n']:.3f}, range {min(degree)} to {max(degree)}")
```

The edge list builds a NetworkX graph directly, and the ground truth is in the
shape the [evaluation metrics](evaluation.md) take:

```python exec="on" source="above" session="benchmarks" result="text"
G = nx.Graph()
G.add_nodes_from(range(params["n"]))
G.add_edges_from(edges)
found = pymocd.rimpso(G)
print(f"AMI of rimpso against the planted communities: {pymocd.ami(found, truth):.4f}")
```

## Parameters

The names and defaults are NetworkX's, and so are the meanings, except for
`tol` and `max_iters`, which in NetworkX also drive a numerical search for the
minimum degree that `lfr_graph` does not need.

| Parameter | Meaning |
|---|---|
| `n` | Number of nodes, labelled `0` to `n - 1`. |
| `tau1` | Exponent of the degree distribution, $p(k) \propto k^{-\tau_1}$; must be finite and greater than 1. |
| `tau2` | Exponent of the community-size distribution, $p(s) \propto s^{-\tau_2}$; must be finite and greater than 1. |
| `mu` | Mixing parameter in $[0, 1]$: the fraction of each node's edges that leave its community. |
| `average_degree` | Requested mean degree. Exactly one of `average_degree` and `min_degree` is given. |
| `min_degree` | Smallest degree a node can draw. |
| `max_degree` | Largest degree a node can draw, in $(0, n]$; `None` means `n`. Degrees are capped at `n - 1`, the most a node of a simple graph can have. |
| `min_community` | Smallest community size; `None` means the smallest drawn degree. |
| `max_community` | Largest community size; `None` means the largest drawn degree. |
| `tol` | Slack on `average_degree`: a request up to `tol` outside the reachable range (see below) is accepted and uses the nearer end of that range. A negative `tol` lowers the top of the range by `abs(tol)` and refuses requests below `1 + abs(tol)`, but does not raise the bottom end $E(1)$. Default `1e-07`. In NetworkX it is the tolerance of the bisection that finds the minimum degree. |
| `max_iters` | Draws of degrees, community sizes and assignment before giving up (see [When a draw cannot be completed](#when-a-draw-cannot-be-completed)), and, within a draw, attempts at a degree sequence with an even sum and at a size sequence that sums to `n`, as in NetworkX. Fewer draws are made once the repairs of the rejected draws have used up the work of 32 draws; draws that fail before the repairs or in the wiring use none of that work, and only `max_iters` limits them. In NetworkX it also caps the search for the minimum degree and, multiplied by `10 * n`, the steps of the community assignment; `lfr_graph` has neither. Default `500`. |
| `seed` | Integer seed in $[0, 2^{64})$, or `None` to seed from the operating system. Other integers raise `OverflowError`; NetworkX also accepts negative integers and `random.Random` instances. |

`n` and `max_iters` take Python integers; a float such as `300.0` raises
`TypeError`, as in NetworkX. `min_degree`, `max_degree`, `min_community` and
`max_community` also take a float with an integer value, such as `20.0`, and
use that integer. NetworkX takes floats for them too, except `max_degree`
together with `average_degree`, where it raises `TypeError` and `lfr_graph`
does not. A fractional value such as `5.5` raises `ValueError` here. NetworkX
accepts it, and a fractional lower bound lets its sampler draw values below the
bound: `zipf_rv` computes `x = int(xmin * u ** -(1.0 / a1))` with
`a1 = alpha - 1`, so with `xmin = 5.5` it can return 5. A
value that is not finite raises `ValueError` for every float parameter,
including some that NetworkX accepts (see [Errors](#errors)).

### The degree law behind `average_degree`

Degrees are independent draws from $p(k) \propto k^{-\tau_1}$ on the integers
$[k_{\min}, k_{\max}]$, where $k_{\max}$ is the smaller of `max_degree` and
`n - 1`. With `min_degree`, $k_{\min}$ is `min_degree`. With `average_degree`,
no integer $k_{\min}$ would do: the expected degree $E(m)$ of the law on
$[m, k_{\max}]$ grows with $m$ in steps, and at low averages the steps are
large. For $\tau_1 = 2.5$ and $k_{\max} = 100$, $E(1) = 1.80$, $E(2) = 4.15$ and
$E(3) = 6.46$. `lfr_graph` therefore takes the $j$ with
$E(j) \le$ `average_degree` $< E(j + 1)$ and draws each degree from

$$
(1 - a)\,\text{law}[j, k_{\max}] + a\,\text{law}[j + 1, k_{\max}],
\qquad a = \frac{E(j) - \text{average\_degree}}{E(j) - E(j + 1)},
$$

whose mean is `average_degree` exactly. Above $j$ this is the power law with
exponent $\tau_1$; only the probability of degree $j$ is reduced. A request
outside $[E(1), k_{\max}]$, widened by `tol`, raises `RuntimeError`. A
negative `tol` narrows the range to
$[\max(E(1), 1 + |\text{tol}|), k_{\max} - |\text{tol}|]$: it lowers the top
by $|\text{tol}|$ but leaves $E(1)$ where it is. The block
below computes the law for three requests and compares its mean, and the mean
degree of five generated graphs, with the request and with the nearest integer
minimum degree:

```python exec="on" source="above" session="benchmarks" result="text"
def truncated_mean(tau, lo, hi):
    weights = [k ** -tau for k in range(lo, hi + 1)]
    return sum(k * w for k, w in zip(range(lo, hi + 1), weights)) / sum(weights)


for tau1, max_degree, average in ((2.5, 100, 3), (2.0, 100, 5), (2.5, 50, 20)):
    E = {m: truncated_mean(tau1, m, max_degree) for m in range(1, max_degree + 1)}
    j = max(m for m, e in E.items() if e <= average)
    a = (E[j] - average) / (E[j] - E[j + 1])
    nearest = min(E, key=lambda m: abs(E[m] - average))
    total = 0
    for seed in range(5):
        edges, _ = pymocd.lfr_graph(n=2000, tau1=tau1, tau2=1.5, mu=0.5, average_degree=average,
                                    max_degree=max_degree, min_community=60, max_community=200,
                                    seed=seed)
        total += 2 * len(edges)
    print(f"tau1={tau1}, max_degree={max_degree}, average_degree={average}: j={j}, a={a:.4f}, "
          f"law mean {(1 - a) * E[j] + a * E[j + 1]:.6f}, realized {total / (5 * 2000):.3f}; "
          f"nearest integer minimum {nearest} gives {E[nearest]:.3f} "
          f"({100 * (E[nearest] / average - 1):+.1f}%)")
```

To fix the degree law to a pure power law, pass `min_degree` instead.

## What the generator guarantees

For every graph `lfr_graph` returns:

- **A simple graph.** There are no self-loops and no repeated edges.
- **The drawn degrees, exactly.** Each node's degree is the value it drew from
  the degree law, whose mean is `average_degree` when that is given. The
  degree sequence is redrawn until its sum is even, as in NetworkX.
- **The requested mixing, node by node.** A node of degree $k$ has
  $\lfloor (1-\mu) k \rfloor$ or $\lceil (1-\mu) k \rceil$ edges inside its
  community, and exactly $(1-\mu)k$ when that is an integer. A node is
  rounded up with probability equal to the fractional part, so the expected
  number of internal edge ends is $(1-\mu)\sum k$, unless the communities are
  too small for their members' ceilings (see
  [When nothing can compensate](#when-nothing-can-compensate)). The repairs
  below can move a node's rounding to another node; where that leaves some
  nodes further from $\mu$ than chance, [Limits](#limits) measures it. A
  community whose internal degrees have an odd sum cannot be wired. When the total over all
  communities is odd, one node is rounded the other way, up or down by a fair
  coin; the odd communities are then paired, and each pair exchanges a member
  of odd internal degree for one of even internal degree, which changes no
  node's internal degree. A pair that no exchange fixes rounds one member up
  and one down. None of these steps moves the expected internal total.
- **One exception to the rounding rule.** When a community's internal degrees
  add up to an odd number and none of its members can be rounded the other
  way (every $(1-\mu)k$ in it is an integer, for instance), one of its nodes
  moves one edge across and ends one edge more or one fewer inside than its
  rounding allows. A later exchange of members can move that node to another
  community. With `n=99`, every degree 10 and $\mu = 0.5$, for instance,
  three communities of 33 nodes each need $33 \times 5$ internal edge ends, an
  odd number. At $\mu = 0$ the community's nodes are assigned again instead.
- **Community sizes in range.** Sizes are drawn from the truncated power law with
  exponent `tau2` and sum to `n` exactly. Every node's internal degree is
  smaller than its community's size.
- **The two ends of $\mu$.** With $\mu = 0$ no edge joins two communities, and
  with $\mu = 1$ no edge lies inside one.

A node goes to a community larger than $\lceil (1-\mu) k \rceil$, so that
either rounding fits. NetworkX instead admits a node to a community larger than
`round((1 - mu) * k)`. When no assignment exists on the ceilings, `lfr_graph`
lets the fewest nodes needed into communities larger than only their floor,
choosing among the nodes whose $(1-\mu)k$ NetworkX itself rounds down.
Whenever NetworkX can place a draw of degrees and sizes, so can `lfr_graph`.
The two do not admit the same nodes, though: when every node fits on its
ceiling, `lfr_graph` never seats a node in a community no larger than
$\lceil (1-\mu)k \rceil$, which NetworkX does. A node admitted on its floor
that rounded up and lands where its ceiling does not fit is rounded down, and
another node that rounded down is rounded up in its place: a member of the
same community when one can be, which keeps the community's expected internal
total; otherwise a member of another community of the same size; otherwise
any node of the communities being placed, which are all of them except when a
repair places a few communities again. When none can be, the internal edge end
is lost, and the fraction of edges between communities ends above $\mu$ (see
[When nothing can compensate](#when-nothing-can-compensate)). Among the
placements in which every node fits, the one used is drawn uniformly at
random.

## Limits

The rounding rule above holds for every node but not, in every case, with the
probabilities it states: some draws force a node to one side, and the repairs
move rounding between nodes. This section measures where that shows. Every
number was measured on jackal (48 cores) with this version of pymocd; $z$ is a
deviation divided by its standard error under the stated rounding.

### Small communities at the default `min_community`

A community of $s$ nodes can hold a node of degree $k$ on its ceiling only if
$\lceil (1-\mu)k \rceil < s$. With `min_community` left at its default, the
smallest drawn degree $k_{\min}$, communities of $k_{\min}$ nodes are drawn,
and when $0 < \mu < 1/k_{\min}$ the ceiling of $(1-\mu)k_{\min}$ is $k_{\min}$
itself: such a community cannot hold a node of degree $k_{\min}$ on its
ceiling. Where NetworkX's `round((1 - mu) * k)` is the floor, as in the
requests below, `lfr_graph` admits such nodes on their floor (see above), and
every one of them gets $\lfloor (1-\mu)k \rfloor$ internal edges. No member of
their community can be rounded up in their place, so nodes of other
communities are rounded up instead. In the requests below enough of them can
be: that keeps the expected number of internal edge ends, so the fraction of
edges between communities stays at $\mu$, but the forced nodes mix more than
$\mu$ and the others less, and the mixing averaged over nodes ends above
$\mu$. Communities one node larger can force nodes one degree up in the same
way. When too few nodes can be rounded up, the fraction itself ends above
$\mu$ (see [When nothing can compensate](#when-nothing-can-compensate)).

Ten graphs per request (seeds 0 to 9), $n = 10^5$, $\mu = 0.1$, every other
parameter at its default. A node is counted as forced when its community
cannot hold its ceiling.

| Request | Forced nodes | Their mixing | Fraction of edges between communities | Mixing averaged over nodes | Mixing of the other nodes |
|---|---|---|---|---|---|
| `tau1=3.5, tau2=3, min_degree=5` | 18.2%: degree 5 in communities of 5 | 0.2000 | 0.10003 | 0.10403 | 0.0826; by community size from 0.0733 (6 nodes) to 0.0887 (above 50) |
| `tau1=3.5, tau2=3, average_degree=10` ($k_{\min} = 6$) | 23.8%: degree 6 in communities of 6, degree 7 in communities of 7 | 0.1582 | 0.10000 | 0.10348 | 0.0864; 0 for the nodes of degree 6 in communities of 7 (all on their ceiling), from 0.0841 (8 to 10 nodes) to 0.0927 (above 50) for the rest |
| `tau1=3.5, tau2=1.5, min_degree=5` | 1.1%: degree 5 in communities of 5 | 0.2000 | 0.09997 | 0.10015 | 0.0990; by community size from 0.0986 to 0.0991 |

Among the nodes that are not forced, in the first request, those of degree 5
sit on their ceiling 63.4% of the time against a fractional part of 0.5, those
of degree 6 53.4% against 0.4, degree 7 42.3% against 0.3, degree 8 29.3%
against 0.2 and degree 9 15.7% against 0.1; in the second, every node of
degree 6 sits on its ceiling (fraction 0.4) and those of degree 7 43.8% of the
time (0.3). In the third, degrees 5 to 9 sit on their ceiling within 0.0074 of
their fractional part. The first
request at $\mu = 0.3$, where $\lceil 0.7 \times 5 \rceil = 4 < 5$, forces no
node: the mixing of every community-size class lies between 0.2997 and 0.3002,
and each degree sits on its ceiling as often as its fractional part (0.500
against 0.5, 0.199 against 0.2, 0.900 against 0.9).

### When nothing can compensate

A node held to its floor is compensated by rounding up another node that
rounded down and has room for its ceiling. When the communities are so small
that few nodes anywhere have that room, the internal edge ends of the forced
nodes are lost, and the fraction of edges between communities ends above
$\mu$. `lfr_graph` keeps such a graph, since every draw of the same request
has communities just as small, and emits one `RuntimeWarning` per call when
the loss raises the expected fraction by more than 0.001. The warning gives
$\mu$, the realized fraction, the shift, and the remedy: a larger
`min_community`, with `max_community` at least as large. For example, every degree 2 at $\mu = 0.3$ with the default
community bounds draws communities of 2 nodes, which hold one internal edge
per node, and the 40% of nodes that round 1.4 up are held to 1:

```python exec="on" source="above" session="benchmarks" result="text"
import warnings

with warnings.catch_warnings(record=True) as caught:
    warnings.simplefilter("always")
    pairs, pair_truth = pymocd.lfr_graph(n=10_000, tau1=2.5, tau2=1.5, mu=0.3, min_degree=2,
                                         max_degree=2, seed=1)
between = sum(pair_truth[u] != pair_truth[v] for u, v in pairs) / len(pairs)
print(f"fraction of edges between communities: {between:.4f} (mu=0.3)")
for w in caught:
    print(f"{w.category.__name__}: {w.message}")
```

The table gives the requests of this kind measured, $n = 10^4$ unless stated,
every other parameter at its default unless stated. The last column is the
realized fraction less $\mu$ and the shift the warning states, which leaves
rounding noise only.

| Request | Graphs (seeds) | Fraction of edges between communities | Warned | Shift stated | Realized − ($\mu$ + shift) |
|---|---|---|---|---|---|
| every degree 2, `tau1=2.5, tau2=1.5`, $\mu = 0.3$ | 5 (0 to 4) | 0.5000 in each | 5 | 0.1958 to 0.2020 | −0.0020 to +0.0042 |
| `tau1=3.5, tau2=3`, degrees 2 to 3, $\mu = 0.3$ | 5 (0 to 4) | 0.3970 to 0.4076 | 5 | 0.0934 to 0.1072 | +0.0004 to +0.0036 |
| `tau1=3, tau2=3`, degrees 3 to 5, $\mu = 0.2$ | 5 (0 to 4) | 0.2365 to 0.2523 | 5 | 0.0373 to 0.0514 | −0.0015 to +0.0016 |
| `tau1=3.5, tau2=3`, degrees 5 to 8, $\mu = 0.1$, $n = 10^5$ | 10 (0 to 9) | 0.1250 to 0.1289 | 10 | 0.0249 to 0.0287 | −0.0004 to +0.0003 |
| every degree 10, communities of 10, `tau1=2, tau2=1.5`, $\mu = 0.06$, $n = 2000$ | 10 (0 to 9) | 0.1000 in each | 10 | 0.0387 to 0.0420 | −0.0020 to +0.0013 |

The fourth request differs from the first request of
[Small communities at the default `min_community`](#small-communities-at-the-default-min_community)
only in `max_degree=8`. With `max_degree` at 10 or 20 instead, enough nodes
have room, and over the same ten seeds the fraction was 0.10008 and 0.09992
with no warning; the three requests of that section gave 0.10003, 0.10000 and
0.09997, also with no warning. Over the 160 cells of a sweep at $n = 2000$ (`min_degree` of 3, 5,
10 and 20, `max_degree` of one more, twice, five times or the default,
`tau1` of 2 and 3, `tau2` of 1.5 and 2, every $\mu$ of 0.02, 0.05, 0.1 and 0.2
below $1/k_{\min}$, five seeds each), 48 cells returned graphs, 211 in all.
In 9 of those cells the internal total missed $(1-\mu)\sum k$ by more than 4.4
standard errors. In 8 of them the miss was large ($z = -22.9$ to $-63.0$) and
every graph warned; with the edge ends the warnings state added back, their
$|z|$ was at most 1.6. The ninth (`min_degree=3`, `max_degree=6`, `tau1=2`,
`tau2=2`, $\mu = 0.2$) sat just past the line at $z = -4.4$: 3 of its 5 graphs
warned, and the other 2 had no unpaid edge end, so the miss there is within
the rounding noise of a cell that small. The 504 graphs of the
acceptance grid (see
[When a draw cannot be completed](#when-a-draw-cannot-be-completed)) raised no
warning.

### Mixing by community size

The repairs that make communities wirable move nodes between communities
according to their internal degrees, and a node's internal degree depends on
how it was rounded. A community that is not graphical, typically one that
holds hubs, takes in members of high internal degree, which rounded up more
often than not. `lfr_graph` therefore looks for exchange partners among the
communities of the same size first, which moves rounding between communities
but not between sizes. The rounding that exchanges between different sizes do
move is booked per community and given back at the end, by exchanging members
of equal internal degree and different degree, which changes no community's
internal degrees and no node's.

What is left was measured on a heavy-tailed set, `tau1=2`, `tau2=1.5`,
`average_degree=5`, `max_degree=50`, communities of 10 to 50 nodes, at
$\mu$ of 0.05, 0.3 and 0.7 and $n$ of 2,000 (20,000 graphs per $\mu$), 5,000
(8,000), 10,000 (4,000) and $10^5$ (400), about 40 million nodes for each $n$
and $\mu$, each block from its own seeds. Nodes are grouped by the size of
their community; each cell gives the mean of $k_{in} - (1-\mu)k$ per node and
its $z$ against the rounding variance $\sum f(1-f)$, $f$ the fractional part
of $(1-\mu)k$. One standard error is about 0.00013 in every cell of a single
$n$ and 0.00006 in the pooled rows.

| $\mu$ | $n$ | 10 to 14 | 15 to 19 | 20 to 29 | 30 to 39 | 40 to 50 |
|---|---|---|---|---|---|---|
| 0.05 | 2,000 | +0.00008 ($z$ +0.64) | +0.00014 (+0.93) | +0.00060 (+5.34) | −0.00001 (−0.11) | −0.00083 (−6.54) |
| 0.05 | 5,000 | +0.00004 (+0.30) | +0.00034 (+2.33) | +0.00020 (+1.79) | +0.00022 (+1.79) | −0.00041 (−3.22) |
| 0.05 | 10,000 | +0.00012 (+0.94) | +0.00004 (+0.30) | +0.00012 (+1.07) | +0.00020 (+1.64) | −0.00025 (−1.99) |
| 0.05 | $10^5$ | +0.00015 (+1.17) | +0.00003 (+0.20) | −0.00007 (−0.64) | +0.00005 (+0.39) | −0.00019 (−1.52) |
| 0.05 | pooled | +0.00010 (+1.53) | +0.00014 (+1.88) | +0.00021 (+3.76) | +0.00011 (+1.86) | −0.00042 (−6.68) |
| 0.3 | pooled | +0.00004 (+0.47) | −0.00005 (−0.55) | +0.00020 (+2.92) | −0.00002 (−0.30) | −0.00001 (−0.15) |
| 0.7 | pooled | +0.00006 (+0.77) | −0.00008 (−0.88) | −0.00005 (−0.70) | −0.00001 (−0.18) | −0.00003 (−0.35) |

At $\mu = 0.05$ the communities of 40 to 50 nodes, which hold the hubs, end
with fewer internal edges than $(1-\mu)k$ and those of 20 to 29 nodes with
more. The largest residual is 0.00083 internal edges per node, which at the
mean degree of 5 moves those nodes' mixing by about 0.0002. It shrinks with
$n$: in the communities of 40 to 50 nodes it is beyond five standard errors at
$n = 2000$, about three at $n = 5000$ and within two at $n = 10{,}000$ and
$10^5$. At $\mu = 0.3$ the largest $|z|$ for a single $n$ is 2.88 (20 to 29
nodes at $n = 5000$), and at $\mu = 0.7$ it is 2.16.

### Overall mixing at the extremes

The pooled mixing was measured on 35 settings: the five parameter sets below,
each with $n = 2000$ and $\mu$ of 0.001, 0.01, 0.05, 0.5, 0.95, 0.99 and
0.999.

| Set | `tau1` | `tau2` | `average_degree` | `max_degree` | Communities |
|---|---|---|---|---|---|
| A | 2.5 | 1.5 | 20 | 50 | 20 to 100 |
| G | 3.5 | 1.05 | 5 | 50 | 10 to 50 |
| G2 | 2 | 1.5 | 5 | 50 | 10 to 50 |
| S | 2.5 | 2 | 10 | 30 | 15 to 40 |
| F | 2 | 1.5 | 20 | 50 | 51 to 60 |

For each setting, graphs from seeds $10^6 i$, $10^6 i + 1$, and so on ($i$ the
position of $\mu$ in the list above) were pooled in batches of 100 until the
standard error of the mean realized fraction of edges between communities was
at most $2.5 \times 10^{-5}$, with at least 200 and at most 8,000 graphs;
44,400 graphs in all, every one returned. The mean realized fraction was within
2.73 standard errors of $\mu$ in every setting, and the largest difference was
$6.7 \times 10^{-5}$.

Nodes were also grouped by degree and by community size within each setting.
For each group, $z$ divides the summed deviation of $k_{in}$ from $(1-\mu)k$
by the square root of its variance: the rounding variance $\sum f(1-f)$, plus
one unit per community for the coin flip that makes each graph's internal
total even, shared among the groups in proportion to their nodes. Over the
1,557 groups the largest $|z|$ was 3.03; the largest of 1,557 independent
standard normal $|z|$ has mean 3.55 and median 3.51. Against the rounding
variance alone, over the 1,436 groups where it is positive, the largest $|z|$
was 3.37, against a mean of 3.53 and a median of 3.49 for the largest of that
many.

A node admitted on its floor cannot always keep its expected internal degree:
with `tau1=2`, degrees 3 to 12, communities of 10 to 12 nodes and
$\mu = 0.05$, a node of degree 12 needs $0.95 \times 12 = 11.4$ internal
edges and has room for 11 at most. Its community gives the difference to
another member, whose lower degrees leave room for their ceilings, so the
community's total and the overall mixing are kept, but the nodes of degree 12
have 0.4 fewer internal edges than $(1-\mu)k$. Where no node has room, the
total is not kept, as
[When nothing can compensate](#when-nothing-can-compensate) shows.

Apart from that case, which the warning reports, a request for which no draw
meets these constraints raises an error instead of returning a graph that
misses them (see [Errors](#errors)).

### Draws rejected that could have been wired

The repairs are heuristics: a draw they give up on may have had an exact
wiring, and `lfr_graph` has no exact test for that. An independent search for
such wirings (random roundings, a greedy placement with local repairs,
Havel–Hakimi inside communities and a greedy construction between them, each
wiring found checked edge by edge) was run on the first draws of 160 seeds of
six small requests. It cannot prove that a draw has no wiring, but it found
one for 7 draws that `lfr_graph` rejected, all with two communities of
unequal size:

- `n=250, tau1=3, tau2=1.1, mu=0.2, average_degree=10, max_degree=30`,
  communities of 50 to 200: 48 of 160 first draws rejected, 3 of them with an
  exact wiring (seeds 30, 127 and 158). In seed 30 the 163 members of the
  larger community with the fewest external edge ends hold 249 of the 494, more
  than half, so no exchange of members balances it under the rounding drawn;
  a different rounding would.
- `n=400, tau1=2, tau2=1.5, mu=0.05, min_degree=3, max_degree=60`,
  communities of 40 to 400: 26 of 160 first draws rejected; the search could
  reconstruct 11 of them and found an exact wiring for 4 (seeds 34, 69, 91 and
  143). In seed 34 the exchanges ran out before the larger community gave up
  enough external edge ends.
- Four requests with two to four communities of nearly equal size (`n` of 200
  to 1,000): every first draw of 160 seeds was wired.

Such draws are replaced like any other failed draw, so they add to the
selection effect below.

## When a draw cannot be completed

Some draws of degrees and community sizes cannot be wired exactly. A hub may
need a community larger than any drawn, or one community may hold more than
half of all the edge ends that leave a community, so that they cannot all be
paired with ends outside it. `lfr_graph` then draws a new degree sequence, new
community sizes and a new assignment, up to `max_iters` draws, and raises
`RuntimeError` with the message of the stage that failed in the last one.
NetworkX raises on the first such draw, or loops forever.

The repairs that try to save a draw have a budget: 5,000 steps per edge end
of the draw, plus a constant, and 32 such draws' worth for the whole request.
A request whose draws keep failing in the repairs therefore stops after the
equivalent of about 32 of them, even when `max_iters` allows more, and fails
in time roughly proportional to the size of the graph (see
[Errors](#errors)). The budget covers the repairs only: a draw that fails
before them, because no placement exists, uses none of it, a draw that fails
afterwards, in the wiring, has used only the repairs it ran first, and
`max_iters` alone limits how many such draws are made. In the
successful draws measured, the largest spent 3,834 steps per edge end
(`tau1=2`, `average_degree=10`, `min_community=20`, $n = 10^5$, $\mu = 0.3$).

The graphs returned are therefore conditioned on their draw having been
completed, and where many draws fail, their degrees lean away from the
requested law: draws with large hubs fail more often. The table measures the
effect at $\mu = 0.1$ for NetworkX's docstring example and at $\mu = 0.3$ for
the others, with every other parameter at its default. The mean degree pools
all nodes of all graphs; its $z$ is against the degree law's standard
deviation over that many nodes.

| Request | Graphs (seeds) | Draws rejected | Mean degree | Fraction of edges between communities | Mixing averaged over nodes |
|---|---|---|---|---|---|
| NetworkX's docstring example: `n=250, tau1=3, tau2=1.5, average_degree=5, min_community=20` | 5,000 (0 to 4,999) | 4,135 of 9,135 (45.3%) | 4.885, −2.30% ($z$ −22.3) | 0.0999 | 0.0998 |
| `n=10000, tau1=2, tau2=1.5, average_degree=10, min_community=20` | 400 (0 to 399) | 3,529 of 3,929 (89.8%) | 8.909, −10.9% ($z$ −20.9) | 0.3000 | 0.3000 |
| `n=100000`, otherwise the same | 100 (0 to 99) | 810 of 910 (89.0%) | 9.063, −9.4% ($z$ −10.1) | 0.3000 | 0.3000 |
| `n=10000, tau1=2.5, tau2=1.5, average_degree=20, max_degree=50`, communities of 20 to 100 | 400 (0 to 399) | none | 20.007, +0.04% ($z$ +1.6) | 0.3000 | 0.3000 |

The most draws one graph needed was 12 in the first row, 70 in the second and
39 in the third. A lower `max_degree` lowers the largest degrees, and with
them the rejections: the request of the second row with `max_degree=1000`
had 12 of 212 draws rejected over seeds 0 to 199 (5.7%), and with
`max_degree=100` none of 200.

A draw's community sizes are rejected together with its degrees, so they lean
too. With `n=20000, tau1=2.5, tau2=3, min_degree=5` and the other bounds at
their defaults, seeds 0 to 47 gave 46 graphs at $\mu = 0.4$ (two calls raised
an error) whose mean community size was 10.03, against 9.05 under the
size law given each graph's degree range (+10.9%; a chi-square test of the
pooled size counts against the law gives $p = 5 \times 10^{-246}$ when it also
counts the number of communities, and $p = 2 \times 10^{-69}$ for the shape
alone), and a mean largest degree of 2,651. At $\mu = 1$, where no node
needs room in its community, the same seeds gave 48 graphs with a mean size of
9.04 against 9.05 ($p = 0.87$) and a mean largest degree of 4,908.

A `min_degree` request leans the same way, and since it states no mean, no
warning reports it. With `n=2000, tau1=2, tau2=1.05, min_degree=3`,
communities of 50 to 500 and $\mu = 0.3$, seeds 0 to 199 gave 172 graphs (the
others raised `RuntimeError`) whose pooled mean degree was 14.08, 16.9% below
the law's 16.93; at $\mu = 1$ the same seeds gave 200 graphs 3.2% below it.
With `n=500` and communities of 50 to 125 the two were −20.3% and −3.8%, and
with `n=10000`, `tau2=1.5`, `min_community=20` and the other bounds at their
defaults −9.4% and −2.6%. None of these graphs raised a warning.
Over the grid of the acceptance tests (`n` of 1,000 and 10,000, `tau1` of 2
and 3, `tau2` of 1.1, 1.5 and 2, `average_degree` of 10 and 20,
`max_degree=50`, communities of 20 to 100 nodes, $\mu$ from 0 to 1, three
seeds each, 168 requests), the mean degree was within 2.4 standard errors of
`average_degree` in every request.

When `average_degree` is given and the mean degree of the returned graph
differs from it by more than four standard errors of the degree law (its
standard deviation over $\sqrt{n}$) and by more than 2%, `lfr_graph` emits a
`RuntimeWarning` that states both numbers and how many draws were rejected.
The standard error is that of a single graph, and with `tau1=2` and the
default `max_degree` the law's standard deviation is large: 104.6 at
$n = 10^4$ and 293.4 at $n = 10^5$, so four standard errors are 4.19 and 3.71
degrees. None of the graphs of the two middle rows crossed that, though their
means ranged from 7.24 to 10.59 and from 7.56 to 9.91: the warning is meant
for misses well beyond what a single draw of the law produces. A graph from a
single draw still crosses both margins about once in 11,800 graphs: 12 of the
141,600 graphs of the size-class and extremes measurements in
[Limits](#limits) did, each from its first draw, and the warning says so.

## How it differs from NetworkX

`lfr_graph` keeps NetworkX's parameters, defaults and validation. It also keeps
the community-size law and the degree law for `min_degree`. It changes eight
things.

**The wiring.** NetworkX
[wires each node in turn](https://github.com/networkx/networkx/blob/networkx-3.4.2/networkx/generators/community.py#L1058-L1070).
It adds random neighbours while `G.degree(u)` is below the node's internal
target, then while it is below the node's degree. `G.degree(u)` counts every
edge `u` already has, including the external edges that nodes wired earlier
attached to it, so a node wired late can reach its internal target with few
internal edges. Each new neighbour `v` is drawn without checking whether `v`
already has its full degree, so a node wired early is pushed past its degree.
An internal neighbour is drawn from the whole community, `u` included, so
self-loops appear. `lfr_graph` follows the original algorithm instead. It fixes
every node's internal and external degree first, then wires them exactly:
Havel–Hakimi inside each community and a configuration model between
communities, each randomized by degree-preserving edge swaps. The internal
degrees of a community have to form a sequence that a simple graph can
realize, and the external edge ends have to pair up across communities;
NetworkX checks neither. `lfr_graph` repairs a community that fails by
exchanging members with other communities. When one community holds more than
half of the external edge ends, it exchanges that community's members of
highest external degree for outsiders of lower external degree. When the
external edges are dense and spread over few communities, the configuration
model's repair stalls (with two communities of 50 nodes and every node of
external degree 49, nearly every rewiring repeats an edge), and the external
edges are then built greedily, as Havel–Hakimi does: the community with the
most edge ends left gives its node with the most left, which is joined to the
nodes with the most edge ends left in other communities. With two communities
this construction succeeds whenever a graph exists, whichever node is joined
first.

**The rounding of $(1-\mu)k$.** NetworkX uses Python's `round`
([line 792](https://github.com/networkx/networkx/blob/networkx-3.4.2/networkx/generators/community.py#L792)),
which rounds half to even and always the same way. When many nodes share a
degree, the realized mixing then sits off $\mu$. `lfr_graph` rounds at random,
as described above.

**Where a node may sit.** NetworkX seats a node in a community larger than
`round((1 - mu) * k)`
([line 796](https://github.com/networkx/networkx/blob/networkx-3.4.2/networkx/generators/community.py#L796)).
`lfr_graph` asks for a community larger than the ceiling, so that either
rounding fits, and falls back on the floor only for the fewest nodes needed,
among those NetworkX rounds down (see
[What the generator guarantees](#what-the-generator-guarantees)).

**How nodes are placed.** NetworkX follows the paper's homeless-node
procedure: a node joins a random community large enough for it, and a
community that overflows evicts a random member, for at most
`10 * n * max_iters` steps. `lfr_graph` draws the placement directly,
uniformly among those in which every node fits, seating the most demanding
nodes first, so it has no step budget to run out of. With `tau1=2, tau2=1.1,
mu=0.5, min_degree=1, max_degree=9, min_community=1, max_community=50`, the
communities of one node only take nodes admitted on their floor, and exactly
as many nodes are admitted that way as there are such seats. `lfr_graph`
generated graphs of 25,000 to 200,000 nodes for it in time proportional to
`n`, a quarter to a third of a second of CPU time each at 200,000 nodes (0.23
to 0.33 s over two runs of three seeds on jackal, a 48-core server shared with
other jobs).

**Failed draws.** NetworkX stops at the first draw it cannot place, or loops
forever on it. `lfr_graph` draws again (see
[When a draw cannot be completed](#when-a-draw-cannot-be-completed)).

**The degree law for `average_degree`.** NetworkX's
[`_generate_min_degree`](https://github.com/networkx/networkx/blob/networkx-3.4.2/networkx/generators/community.py#L726-L755)
bisects a real-valued minimum against an approximate mean, and then rounds it.
The approximation divides a sum truncated at `max_degree` by an untruncated
Hurwitz zeta function, and the rounded minimum then fixes a pure power law, so
the mean degree misses the request. For its docstring example it picks a
minimum degree of 2, whose law has mean 3.17 against the requested 5; for some
feasible requests it raises `ExceededMaxIterations`. `lfr_graph` computes the
exact means of the truncated laws and draws from the mixture described
[above](#the-degree-law-behind-average_degree), whose mean is the request.

**Infeasible requests.** NetworkX loops forever on some requests, such as
`min_degree` above `max_degree`, or a node that needs more external neighbours
than there are nodes outside its community
([networkx#6763](https://github.com/networkx/networkx/issues/6763)). `lfr_graph`
checks for these and raises an error (see [Errors](#errors)).

**Community sizes.** NetworkX draws sizes until their sum reaches `n` and keeps
the draw only if the sum is exactly `n`, at most `max_iters` times. With
communities of thousands of nodes every try can miss. `lfr_graph` makes the
same tries, then samples the same conditioned law another way. It tabulates,
for each remainder $r$, the total probability of the size sequences summing to
$r$, draws sizes freely until the remainder is within the table, and keeps that
draw with probability proportional to the total probability of its
completions. The table stops once its values flatten, or after $2 \times 10^9$
products. In a test at `n` $= 10^6$, size ranges up to `[100, 100000]` drew
their sizes on every one of 10 seeds, each in under a second of CPU time on
jackal.

The block below generates graphs with both generators, from the same parameters
as the example above, and measures what each returns. A self-loop counts twice
toward its node's degree and as an edge inside the community, as NetworkX counts
it.

```python exec="on" source="above" session="benchmarks" result="text"
def summary(n, edges, community):
    degree, outside = [0] * n, [0] * n
    loops = 0
    for u, v in edges:
        if u == v:
            loops += 1
            degree[u] += 2
            continue
        degree[u] += 1
        degree[v] += 1
        if community[u] != community[v]:
            outside[u] += 1
            outside[v] += 1
    mixing = sum(o / d for o, d in zip(outside, degree)) / n
    return sum(degree) / n, min(degree), max(degree), mixing, loops


print(f"networkx {nx.__version__}; requested mean degree {params['average_degree']}, "
      f"max_degree {params['max_degree']}")
print(f"{'mu':>4}  {'generator':<9} {'mean degree':>11} {'min':>4} {'max':>4} "
      f"{'mixing':>7} {'self-loops':>10}")
for mu in (0.1, 0.3, 0.5, 0.7):
    p = dict(params, mu=mu)
    G = nx.LFR_benchmark_graph(**p, seed=1)
    graphs = {
        "networkx": (list(G.edges()), {v: min(G.nodes[v]["community"]) for v in G}),
        "pymocd": pymocd.lfr_graph(**p, seed=1),
    }
    for name, (e, community) in graphs.items():
        mean, low, high, mixing, loops = summary(p["n"], e, community)
        print(f"{mu:>4}  {name:<9} {mean:>11.3f} {low:>4} {high:>4} {mixing:>7.4f} {loops:>10}")
```

The mixing column is each node's fraction of edges leaving its community,
averaged over nodes. The `min` and `max` columns are the smallest and largest
degree in the graph.
[LFR fidelity](../examples/benchmarks/plot_lfr_fidelity.md) in the gallery
plots this comparison over the whole range of $\mu$.

## Errors

Invalid parameters raise `ValueError`. The first five checks are NetworkX's, in
its order and with its messages. The others guard inputs on which NetworkX
3.4.2 loops forever, fails inside its sampler with a message that does not name
the parameter, or builds a graph from a value no power law is defined by; the
bound on `n` is pymocd's own. The NetworkX column for NaN and infinite values
was measured on jackal with `n=500, tau1=2.5, tau2=1.5, mu=0.3,
max_degree=50`, communities of 60 to 200 nodes, `seed=1`, and either
`average_degree=10` or `min_degree=5`; it is the same for both unless the cell
says otherwise, and "no result" means none after 60 seconds.

| Condition | Message | NetworkX 3.4.2 |
|---|---|---|
| `tau1` not greater than 1 (NaN and $-\infty$ included) | `tau1 must be greater than one` | the same |
| `tau2` not greater than 1 (NaN and $-\infty$ included) | `tau2 must be greater than one` | the same |
| `mu` outside $[0, 1]$ (or NaN, or infinite) | `mu must be in the interval [0, 1]` | the same |
| `max_degree` outside $(0, n]$ (or NaN, or infinite) | `max_degree must be in the interval (0, n]` | the same |
| both or neither of `average_degree` and `min_degree` | `Must assign exactly one of min_degree and average_degree` | the same |
| `tau1` or `tau2` $+\infty$ | `tau1 must be finite` (or `tau2`) | no result |
| `n` below 1 | `n must be a positive integer` | depends on the other parameters (measured for `n` of 0, −1, −2, −3 and −10, `max_degree` at its default): with `min_degree=5` and default community bounds, `ValueError: min() arg is an empty sequence`; with `min_degree=5` and communities of 60 to 200, `ExceededMaxIterations: Could not assign communities` for `n=0` and `Could not create power law sequence` for negative `n`; with `average_degree=5`, `ZeroDivisionError: 0.0 cannot be raised to a negative power` for `n=0` and `TypeError: '>' not supported between instances of 'complex' and 'int'` for negative `n` |
| `n` of $2^{32}$ or more | `n must be smaller than 2**32` | no limit |
| `min_degree` NaN or infinite | `min_degree must be finite` | NaN: `ValueError: cannot convert float NaN to integer`; $+\infty$: `OverflowError: cannot convert float infinity to integer`; $-\infty$: `ValueError: xmin < 1` |
| `min_degree`, `max_degree`, `min_community` or `max_community` fractional | `min_degree must be an integer` (and so on) | a graph |
| `min_degree` below 1 | `min_degree must be at least 1` | `ValueError: xmin < 1` |
| `min_degree` above `max_degree` | `min_degree must not exceed max_degree` | loops forever |
| `average_degree` NaN or infinite | `average_degree must be finite` | NaN: a graph; $\pm\infty$: `ExceededMaxIterations: Could not match average_degree` |
| `tol` NaN or infinite | `tol must be finite` | NaN and $+\infty$: a graph; $-\infty$: `ExceededMaxIterations: Could not match average_degree` with `average_degree`, a graph with `min_degree` |
| `min_community` NaN or infinite | `min_community must be finite` | NaN: `ValueError: cannot convert float NaN to integer`; $+\infty$: `OverflowError: cannot convert float infinity to integer`; $-\infty$: `ValueError: xmin < 1` |
| `max_community` NaN or infinite | `max_community must be finite` | NaN and $+\infty$: `ExceededMaxIterations: Could not create power law sequence`; $-\infty$: no result |
| `min_community` below 1 | `min_community must be at least 1` | `ValueError: xmin < 1` |
| `max_community` below 1 | `max_community must be at least 1` | loops forever |
| `min_community` above `max_community` | `min_community must not exceed max_community` | loops forever |

Valid parameters for which no graph is found raise `RuntimeError`, with the
message of the stage that failed in the last draw. A few failures would repeat
in every draw and are raised without drawing again: a degree law whose every
sequence has an odd sum (a single possible degree, odd, with `n` odd), no
size sequence summing to `n` when both `min_community` and `max_community`
are given or `min_community` exceeds `n`, `min_degree` equal to `n`, and a
lowest degree that fits no community any draw can have. Every other failed
draw is drawn again, including one whose `max_iters` attempts at a degree
sequence with an even sum all failed, up to `max_iters` draws, or fewer once
the repairs of the failed draws have used up the work of 32 draws. When both
generators stop for the same reason, `lfr_graph` uses the message of
NetworkX's `ExceededMaxIterations`; the first three rows are those messages.

| Message | Raised when |
|---|---|
| `Could not match average_degree` | `average_degree` is outside the reachable range. |
| `Could not create power law sequence` | No degree sequence with an even sum was drawn in the last draw's `max_iters` attempts, or no community-size sequence sums to `n`: every size is above `n`, or the sizes cannot add up to `n`. Also `min_degree` equal to `n` (possible when `max_degree` is `n`), which no node of a simple graph can have; where NetworkX raises `Could not assign communities` for it, so does `lfr_graph`. |
| `Could not assign communities; try increasing min_community` | Some node's internal degree fits no community, even on the terms NetworkX places nodes on, or the communities' internal degrees could not be made realizable. |
| `Could not wire the inter-community edges; try decreasing max_community` | The external edges could not be placed: one community holds more than half of all external edge ends even after exchanging its members of highest external degree, the external degrees fail the Erdős–Gallai test, or neither the configuration model nor the greedy construction found a simple graph for them. NetworkX has no counterpart: on such draws it returns a graph whose degrees overshoot, or loops forever when a node finds no new neighbour outside its community. |

On some inputs the outcome or the message differs from NetworkX 3.4.2. Each
row below was run with both generators on jackal (`n=500`, `tau1=2.5`,
`tau2=1.5`, `mu=0.3`, `average_degree=10`, `max_degree=50`, communities of 60
to 200 nodes, `seed=1`, unless the row says otherwise):

| Input | `lfr_graph` | NetworkX 3.4.2 |
|---|---|---|
| `average_degree=0` | `RuntimeError: Could not match average_degree` | a graph: its bisection loop never runs and it uses `round((max_degree + 1) / 2)` as the minimum degree (mean degree 44.75 here); the same for `average_degree` or `tol` NaN |
| `tol` negative (`-1`) | a graph: requests above `max_degree - abs(tol)` or below `1 + abs(tol)` are refused, and $E(1)$ stays the bottom of the range when it is higher (with `max_degree=50`, $E(1) = 1.74$: 2.0 and 49.0 give graphs, 1.99 and 49.01 `RuntimeError: Could not match average_degree`) | `ExceededMaxIterations: Could not match average_degree` |
| `average_degree` equal to `max_degree` | a graph in which every degree is `max_degree` | `ExceededMaxIterations: Could not match average_degree` |
| `max_iters=0` | `RuntimeError: Could not create power law sequence` | `ExceededMaxIterations: Could not match average_degree` |
| `max_iters` too small for the size draws (`max_iters=1`, `min_degree=5`) | a graph: the sizes come from the exact sampler | `ExceededMaxIterations: Could not create power law sequence` |
| `max_degree=50.0` with `average_degree` | a graph | `TypeError: 'float' object cannot be interpreted as an integer` |
| `min_degree` above `max_degree` | `ValueError: min_degree must not exceed max_degree` | no result after 60 seconds |
| `min_community` above `n` (`max_community` left at its default) | `RuntimeError: Could not create power law sequence` | no result after 60 seconds |
| a single community and $\mu > 0$ (`n=60`, communities of 60, `min_degree=5`, `max_degree=10`, $\mu = 0.5$) | `RuntimeError: Could not wire the inter-community edges; try decreasing max_community`, after `max_iters` draws | no result after 60 seconds |
| `seed=-1` | `OverflowError` | a graph |

The two generators draw different random numbers and, for `average_degree`,
different degree laws, so even on a request that both can meet they return
different graphs for the same seed, and on requests that only some draws can
meet, NetworkX succeeds on some seeds only.

A request that no draw can meet makes `lfr_graph` try `max_iters` draws, or
the equivalent of 32 draws of repairs, before it raises, so the time it takes
grows with the size of the graph. The times below are CPU time per call, with
the calls run side by side on jackal, a 48-core server shared with other
jobs; a second run on a busier machine took up to 1.6 times as long. With
`tau1=1.5, tau2=1.5, mu=0.5, min_degree=2` and `max_degree` and the community
bounds at their defaults, seeds 1 to 4 all failed: in about 2 seconds each at
`n=5000` (1.9 to 2.4), about 5 at `n=20000` (3.8 to 5.6) and about a minute
at `n=100000` (55 to 57). A request whose draws use up their repair work
takes longest: with `tau1=1.409, tau2=1.637, mu=0.5, min_degree=2` and
`max_degree` a third of `n`, seeds 1 and 3 failed in about a minute at
`n=3000` (53 to 55 s) and two and a half minutes at `n=6000` (140 s), and in
the other run in 81 to 83 s, 165 to 169 s and 414 to 417 s at `n` of 3,000,
6,000 and 12,000: from 18 to 35 ms per node, more the larger `n`. Draws that
fail before using much repair work are bounded by `max_iters` alone:
`n=3000, tau1=1.393, tau2=2.919, mu=0.021, min_degree=3, max_degree=1000`,
communities of 747 to 1,122 nodes, failed in about a second (1.00 to 1.15 s)
on each of five seeds at the default 500 draws, and in proportion to
`max_iters` for seed 1, from 0.10 s at 50 draws to 3.8 s at 2,000. A smaller
`max_iters` fails sooner.

None of these calls can be interrupted with Ctrl-C: `lfr_graph` builds the
graph without holding the GIL and checks no signals, so Python raises
`KeyboardInterrupt` only once the call returns. A SIGINT sent one second into
the `n=20000` call above (seed 1) raised it 4.7 seconds in, when the call
ended.

```python exec="on" source="above" session="benchmarks" result="text"
for bad in (dict(params, mu=1.5), dict(params, min_degree=5)):
    try:
        pymocd.lfr_graph(**bad, seed=1)
    except ValueError as error:
        print(f"ValueError: {error}")
try:
    pymocd.lfr_graph(**dict(params, average_degree=60), seed=1)
except RuntimeError as error:
    print(f"RuntimeError: {error}")
```

## Reproducibility

The same parameters and `seed` give the same graph on the same pymocd build;
`seed=None` seeds from the operating system. Each draw has its own random
stream, seeded from `seed`, the draw's number and the parameters that fix the
degree and size laws: `n`, `tau1`, the degree law (its lowest degree and the
weight $a$ above), the degree cap, `tau2`, `min_community` and
`max_community`. Two calls that agree on these make the same first draw of
degrees and community sizes, the same second draw, and so on. In a sweep over
$\mu$ with one seed, the graphs therefore share their degrees and community
sizes whenever they come from the same draw, which is the first one unless a
draw failed. They also share the random numbers that decide each node's
rounding, so they are not independent replicates: pool graphs across $\mu$
only with a different seed for each $\mu$.

Calls whose laws differ draw independently. Two values of `average_degree`
give different laws unless both are within `tol` of the same end of the
reachable range, where they share that end's law and return the same graph.
Calls that differ only in `tol` return the same graph when both accept the
request, and so do calls that differ in `max_iters`, as long as the call
succeeds within the smaller number of draws and attempts. A seed does not
reproduce NetworkX's graph for the same seed, and unlike NetworkX the seed must
be an integer in $[0, 2^{64})$.

```python exec="on" source="above" session="benchmarks" result="text"
low_mu, _ = measure(params["n"], *pymocd.lfr_graph(**dict(params, mu=0.1), seed=1))
high_mu, _ = measure(params["n"], *pymocd.lfr_graph(**dict(params, mu=0.6), seed=1))
print("same degree sequence at mu=0.1 and mu=0.6:", low_mu == high_mu)
print("same graph from the same seed:",
      pymocd.lfr_graph(**params, seed=1) == pymocd.lfr_graph(**params, seed=1))
```

## See also

- API reference: [`lfr_graph`][pymocd.lfr_graph].
- [Evaluating partitions](evaluation.md): the scores to compare a detector's
  output with the ground truth.
- Gallery: [LFR fidelity](../examples/benchmarks/plot_lfr_fidelity.md), and
  [Recovery as mixing increases](../examples/benchmarks/plot_mixing_sweep.md)
  for a detector sweep on planted-partition graphs.
