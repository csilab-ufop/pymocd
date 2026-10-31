---
description: GDPSO in pymocd — a clean-room port of Cai et al.'s greedy discrete particle swarm, which maximizes Newman–Girvan modularity with label-propagation seeds and greedy single-node moves.
---

# GDPSO

<span class="pm-pill pm-pill--baseline">Re-implemented baseline</span>
<span class="pm-pill pm-pill--neutral">Single-objective</span>
<span class="pm-pill pm-pill--neutral">Deterministic</span>

GDPSO (Greedy Discrete Particle Swarm Optimization) maximizes Newman–Girvan
modularity with a swarm of label vectors. Each particle starts from a short run
of asynchronous label propagation. Once per generation, a binary particle-swarm
velocity decides which nodes are offered a move, and each offered node moves to
the neighboring community with the largest modularity gain when that gain is
strictly positive (nodes of degree 1 are the one exception; see
[Divergences from the paper](#divergences-from-the-paper)). A coarsening
mutation is applied to a fixed fraction of the swarm. The function
returns the single highest-modularity partition the swarm ever held. GDPSO has
one objective, so it produces no Pareto front.

<div class="pm-facts" markdown>

| | |
|---|---|
| **Entry points** | [`gdpso`][pymocd.gdpso] |
| **Objectives** | Newman–Girvan modularity $Q$, maximized (one objective) |
| **Search** | Swarm of label vectors seeded by label propagation; a binary PSO velocity masks which nodes take greedy single-node modularity moves; coarsening mutation |
| **Selection rule** | The best position the swarm ever held (the global best) |
| **Front accessor** | None (single objective); there is no `gdpso_fronts` |
| **Deterministic** | Yes, for a fixed input, and independent of the thread count (module README, measured by tests); the function takes no seed. The same graph (same nodes and edges, listed in the same order) and the same parameters give the same partition on every call; listing the nodes or edges in another order can change it (see [Input order](../guide/performance.md#input-order)). |
| **Original code** | [doctor-cai/GDPSO](https://github.com/doctor-cai/GDPSO) — C++, no license file; this port was written clean-room from a description of it |
| **Paper** | Cai et al. (2015), *Information Sciences* 316:503–516 — [DOI](https://doi.org/10.1016/j.ins.2014.09.041) |

</div>

## Objectives

GDPSO is the method of
[Cai, Gong, Ma, Ruan, Yuan & Jiao (2015)](https://doi.org/10.1016/j.ins.2014.09.041).
It maximizes one objective, the modularity of
[Newman & Girvan (2004)](https://doi.org/10.1103/PhysRevE.69.026113):

$$
Q = \sum_{c} \left[ \frac{e_c}{m} - \left( \frac{K_c}{2m} \right)^{2} \right]
$$

Here $m$ is the number of edges, $e_c$ is the number of edges inside community
$c$ (each counted once), and $K_c = \sum_{i \in c} k_i$ is the degree sum of
$c$, with $k_i$ the degree of node $i$. Edges are counted, not weighted. An
edgeless graph is defined to have $Q = 0$.

The greedy step ranks the candidate moves of a node $i$ with the comparable
part of the exact single-node gain. Let $\kappa_{i,l}$ be the number of edges
from $i$ into community $l$ and $\Sigma_l$ the degree sum of $l$ (the $K_l$
above, `sigma` in the code). Then

$$
\operatorname{score}(l) = \kappa_{i,l} - \frac{k_i \, \Sigma_l}{2m},
\qquad
\Delta Q(i \to l) = \frac{\operatorname{score}(l) - \operatorname{score}(\mathrm{own})}{m},
$$

where $\Sigma_{\mathrm{own}}$ excludes $k_i$ itself. For a fixed node, the score
orders moves the same way as $\Delta Q$, so taking the highest score and
accepting only a strict increase is best improvement with a strictly positive
gain (module README).

## Search

**Representation.** A particle's *position* is a vector of $n$ community
labels, one per node, and the partition is read directly from it. After every
write the labels are renumbered by first appearance, so they are always
$0, \dots, k-1$. A particle's *velocity* is not a displacement: it is a
per-node boolean mask that says which nodes are offered a greedy move in the
current generation. Each particle also keeps its personal best position
(`pbest`); the swarm shares one global best (`gbest`).

**Initialization.** Each particle is seeded by `lpa_sweeps` sweeps of
asynchronous label propagation
([Raghavan, Albert & Kumara, 2007](https://doi.org/10.1103/PhysRevE.76.036106))
from all-singleton labels, in a fixed node order, with ties among the most
frequent neighbor labels broken uniformly at random. Each particle draws from
its own random stream, so different particles can reach different seeds. The
seed is the particle's first personal best, and the best seed becomes the first
global best.

**Velocity.** The update follows the binary particle swarm of
[Kennedy & Eberhart (1997)](https://doi.org/10.1109/ICSMC.1997.637339). For
each node $i$, with $x_i$ the current label, $p_i$ and $g_i$ the personal-best
and global-best labels, $v_i \in \{0, 1\}$ the previous mask bit, $[\cdot]$ the
indicator function, $\sigma$ the logistic sigmoid and $r_1, r_2, u$ fresh
uniform draws on $[0, 1)$:

$$
s_i = w \, v_i + c_1 r_1 \, [x_i \neq p_i] + c_2 r_2 \, [x_i \neq g_i],
\qquad
v_i \leftarrow [\, u < \sigma(s_i) \,].
$$

The defaults for $w$, $c_1$ and $c_2$ are the constriction constants of
[Clerc & Kennedy (2002)](https://doi.org/10.1109/4235.985692). The module
README notes that they carry no convergence theory here, because the velocity
is re-binarized every generation. With non-negative coefficients, $s_i$ lies
in $[0, w + c_1 + c_2]$, which at the defaults is $[0, 3.722]$, so at the
defaults a node's move probability always lies in $[0.5, 0.976]$ (module
README).

!!! note "The swarm terms carry no label information"
    The docstring of [`gdpso`][pymocd.gdpso] notes that `pbest` and `gbest`
    carry no label information: they enter only as the two indicator bits
    above, which shift a node's move probability. The docstring concludes that
    in practice the method behaves as the best of `pop_size` label-propagation
    seeds, each polished by Louvain local moving
    ([Blondel et al., 2008](https://doi.org/10.1088/1742-5468/2008/10/P10008)),
    and that `lpa_sweeps`, not `num_gens`, is the lever on seed diversity. The
    module README adds that the comparison $x_i \neq p_i$ is made on
    renumbered labels, so it is not a distance between partitions; it only
    gates move attempts.

**Generation loop.** In every generation, each particle independently:

1. updates its velocity mask as above;
2. runs a greedy sweep: every masked node, in increasing index order and in
   place, joins the neighboring community with the highest score, with ties
   going to the smallest label, but only on a strictly positive gain (a
   node of degree 1 joins its neighbor's community without this test);
3. if it is one of the first `floor(pop_size * mut_frac)` particles, is
   mutated: each node, with probability `mut_rate`, copies its
   own label onto all of its neighbors;
4. recomputes $Q$ and replaces its personal best if $Q$ improved. There is no
   rejection step; the particle keeps the position it reached.

After all particles have moved, the global best is updated from the swarm.
It is held fixed while the particles of a generation move, which is what lets
them run in parallel without changing the result.

## Selecting a partition

There is no front to select from. [`gdpso`][pymocd.gdpso] returns the global
best: the highest-modularity position any particle held, at initialization or
after any generation. When several positions tie, the one found first is kept,
and within a generation the lowest-numbered particle wins.

The result is a `dict[int, int]` mapping each node to a community numbered
from 0. Isolated nodes are not touched by any operator and are reported as
community `-1`; on a graph with no edges every node is `-1`. With
`num_gens=0` the function returns the best label-propagation seed.

## Parameters

<!-- pymocd:defaults pymocd.gdpso -->
| Parameter | Default | Description |
|---|---|---|
| `pop_size` | `100` | Number of particles. `0` is raised to 1. |
| `num_gens` | `250` | Number of generations. `0` returns the best seed. |
| `w` | `0.7298` | Inertia weight on the previous move mask. |
| `c1` | `1.4961` | Cognitive weight, applied to the indicator $[x_i \neq p_i]$ (node disagrees with its personal best). |
| `c2` | `1.4961` | Social weight, applied to the indicator $[x_i \neq g_i]$ (node disagrees with the global best). |
| `mut_rate` | `0.1` | Probability that a node broadcasts its label inside a mutated particle. |
| `mut_frac` | `0.1` | Fraction of the swarm mutated each generation. |
| `lpa_sweeps` | `5` | Asynchronous label-propagation sweeps that seed each particle. |

Coefficients that are negative or not finite are treated as `0.0`; the two
rates are clamped into $[0, 1]$, and a non-finite rate becomes `0.0`.

The module README records the following about the original implementation's
budgets:

- it ships 100 particles and 250 generations, the defaults here;
- its in-source comments recommend 20–60 particles and 50–300 generations;
- the paper's published tables used 100 particles and 100 generations.

It also reports that on Zachary's karate club the reference's best
modularity over 30 runs is $Q = 0.41979$, the four-group modularity optimum.
With 100 generations, this port stops short of that value at the default
`pop_size=100` and reaches it at `pop_size=300`; a Rust test
(`karate_reaches_the_reference_modularity`) checks both outcomes.

The Rust `Config` has two further fields, `mut_by_index` and
`const_move_prob`. The module README describes them as ablation switches that
are not part of the method. The Python function does not expose them; it
always mutates by index (`mut_by_index = true`, the reference's rule) and
always uses the velocity update described above.

## Complexity

The module README states that modularity is computed in the sparse
$O(n + m)$ form, for $n$ nodes and $m$ edges, instead of the reference's dense
double loop. The code evaluates it once per particle for the seed and once per
particle per generation. Each particle stores a position, a velocity mask and a
personal best, each of length $n$. Neither the module README nor the function
docstring states an overall bound on the running time or memory of a run.
<!-- TODO(docs): overall time/memory complexity of GDPSO — take it from the complexity analysis in Cai et al. (2015), doi:10.1016/j.ins.2014.09.041, or have the maintainer state a bound for this port in src/core/algorithms/gdpso/README.md -->

## When to use it

- **As a published baseline.** GDPSO is one of the re-implemented baselines,
  which pymocd provides for papers and benchmarks. For new analyses, start
  from [Choosing a detector](../guide/choosing.md).
- **When one modularity-maximizing partition is enough.** GDPSO optimizes a
  single scalar and returns one partition, with no front and no
  `gdpso_fronts`. If you need several resolutions or a Pareto front, use a
  multi-objective detector such as [RIMPSO](rimpso.md) or
  [HP-MOCD](hpmocd.md).
- **Not when small communities matter on large sparse graphs.** The module
  README states that GDPSO inherits modularity's resolution limit
  ([Fortunato & Barthélemy, 2007](https://doi.org/10.1073/pnas.0605965104))
  whole, with no mitigation, and that over-coarsening should be expected on
  large sparse graphs. The function docstring also states that GDPSO inherits
  the resolution limit.
- **Tuning.** The docstring and the module README name `lpa_sweeps`, not
  `num_gens`, as the lever on how many distinct seeds the swarm explores. The
  karate result under [Parameters](#parameters) also depends on `pop_size`.

## Divergences from the paper

The port was written clean-room. The module README lists the following as
departures from the paper and from the described behavior of the reference
implementation, not from source code that was read.

- **No pre-shuffle.** The reference shuffles before each label-propagation
  sweep. The port omits it; on the first sweep it would permute pairwise
  distinct labels, so it has no effect there.
- **Seed as first personal best.** The seed counts as the particle's first
  personal best. The reference discards it and starts `pbest` from the first
  position after a velocity update.
- **Two mutation knobs.** `mut_rate` and `mut_frac` are separate parameters.
  The reference uses a single value, 0.1, for both the per-node broadcast
  probability and the fraction of the swarm that is mutated.
- **Which particles mutate.** The module README calls this a documented
  ambiguity. The port follows the reference: the first
  `floor(pop_size * mut_frac)` particles are mutated in every generation.
  Mutating each particle with probability `mut_frac` is available only in
  Rust, as an ablation.
- **Sparse modularity.** $Q$ is computed in $O(n + m)$ instead of the
  reference's dense double loop, and the community degree sums are recomputed
  for each particle in each generation rather than carried across operators,
  which keeps them exact.
- **Graph input.** Self-loops are dropped. A repeated edge is kept and counted
  twice (the multigraph reading); $Q$ and its gain stay consistent either way.
- **Isolated nodes.** They are left out of every operator and reported as
  community `-1` rather than as singleton communities.

The module README also flags the following behaviors, which can look like
bugs:

- In the greedy sweep, a node of degree 1 joins its neighbor's community with
  no gain test. When the neighbor is in another community, this move always
  increases $Q$: with $\kappa_{i,\mathrm{to}} = 1$,
  $\kappa_{i,\mathrm{own}} = 0$ and
  $\Sigma_{\mathrm{to}} \le 2m - 1 - \Sigma_{\mathrm{own}}$, the gain formula
  under [Objectives](#objectives) gives
  $\Delta Q \ge (1 + 2\Sigma_{\mathrm{own}})/(2m^2) > 0$. The missing test
  therefore changes no outcome; the README states that it reproduces the
  reference and must not be changed. The step that can lower $Q$ is the
  mutation, which broadcasts labels without a gain test, and the particle
  keeps the result (there is no rejection step).
- In the mutation, a node of degree 1 adopts its neighbor's label instead of
  broadcasting its own, the opposite of the operator's normal direction.
- The cognitive and social terms transfer no structure (see the note under
  [Search](#search)); the README states that the code deliberately does not
  copy labels. The inertia term is a one-bit memory, because the stored
  velocity is re-binarized every generation.

## Example

The example runs GDPSO on Zachary's karate club at the defaults, at the
larger swarm used by the Rust test above, and with fewer label-propagation
sweeps. It scores each result with the unweighted modularity, because
[`gdpso`][pymocd.gdpso] reads only the edge endpoints, not the weights that
networkx stores on this graph. The last
line checks that a repeated call returns the same partition.

```python exec="on" source="above" session="gdpso" result="text"
import networkx as nx
import pymocd

G = nx.karate_club_graph()

def summarize(partition):
    groups = {}
    for node, community in partition.items():
        groups.setdefault(community, set()).add(node)
    q = nx.community.modularity(G, groups.values(), weight=None)
    return len(groups), q

settings = [{}, {"pop_size": 300, "num_gens": 100}, {"lpa_sweeps": 2}]
for kwargs in settings:
    k, q = summarize(pymocd.gdpso(G, **kwargs))
    label = ", ".join(f"{key}={value}" for key, value in kwargs.items()) or "defaults"
    print(f"{label:<28} communities = {k}   Q = {q:.4f}")

print("repeated call identical:", pymocd.gdpso(G) == pymocd.gdpso(G))
```

Compare the number of communities and $Q$ across the three settings to see
how the swarm size and the seeding sweeps affect the result.

## Reference

Q. Cai, M. Gong, L. Ma, S. Ruan, F. Yuan and L. Jiao, "Greedy discrete
particle swarm optimization for large-scale social network clustering",
*Information Sciences*, vol. 316, pp. 503–516, 2015.
[doi:10.1016/j.ins.2014.09.041](https://doi.org/10.1016/j.ins.2014.09.041)

```bibtex
@article{Cai_2015,
  title={Greedy discrete particle swarm optimization for large-scale social network clustering},
  volume={316},
  ISSN={0020-0255},
  url={http://dx.doi.org/10.1016/j.ins.2014.09.041},
  DOI={10.1016/j.ins.2014.09.041},
  journal={Information Sciences},
  publisher={Elsevier BV},
  author={Cai, Qing and Gong, Maoguo and Ma, Lijia and Ruan, Shasha and Yuan, Fuyan and Jiao, Licheng},
  year={2015},
  month=sep,
  pages={503–516}
}
```

This is the publisher's record from doi.org, except for the month: doi.org
writes `month=Sept`, which is not a BibTeX month macro, so it is given here as
the standard macro `sep`.

Other works cited on this page:

- M. E. J. Newman and M. Girvan, "Finding and evaluating community structure
  in networks", *Physical Review E*, vol. 69, no. 2, 026113, 2004.
  [doi:10.1103/PhysRevE.69.026113](https://doi.org/10.1103/PhysRevE.69.026113)
- J. Kennedy and R. C. Eberhart, "A discrete binary version of the particle
  swarm algorithm", in *Proc. 1997 IEEE International Conference on Systems,
  Man, and Cybernetics*, vol. 5, pp. 4104–4108, 1997.
  [doi:10.1109/ICSMC.1997.637339](https://doi.org/10.1109/ICSMC.1997.637339)
- M. Clerc and J. Kennedy, "The particle swarm — explosion, stability, and
  convergence in a multidimensional complex space", *IEEE Transactions on
  Evolutionary Computation*, vol. 6, no. 1, pp. 58–73, 2002.
  [doi:10.1109/4235.985692](https://doi.org/10.1109/4235.985692)
- U. N. Raghavan, R. Albert and S. Kumara, "Near linear time algorithm to
  detect community structures in large-scale networks", *Physical Review E*,
  vol. 76, no. 3, 036106, 2007.
  [doi:10.1103/PhysRevE.76.036106](https://doi.org/10.1103/PhysRevE.76.036106)
- V. D. Blondel, J.-L. Guillaume, R. Lambiotte and E. Lefebvre, "Fast
  unfolding of communities in large networks", *Journal of Statistical
  Mechanics: Theory and Experiment*, vol. 2008, no. 10, P10008, 2008.
  [doi:10.1088/1742-5468/2008/10/P10008](https://doi.org/10.1088/1742-5468/2008/10/P10008)
- S. Fortunato and M. Barthélemy, "Resolution limit in community detection",
  *Proceedings of the National Academy of Sciences*, vol. 104, no. 1,
  pp. 36–41, 2007.
  [doi:10.1073/pnas.0605965104](https://doi.org/10.1073/pnas.0605965104)

## See also

- API: [`gdpso`][pymocd.gdpso]
- [Algorithms overview](index.md) and the other single-objective baseline,
  [CDRME](cdrme.md)
- [Core concepts](../guide/concepts.md) and
  [Choosing a detector](../guide/choosing.md)
- [Evaluating partitions](../guide/evaluation.md) and
  [Performance and reproducibility](../guide/performance.md)
- Gallery: [Detector tour](../examples/basics/plot_detector_tour.md),
  [Comparing detectors](../examples/benchmarks/plot_compare_detectors.md) and
  [Recovery as mixing increases](../examples/benchmarks/plot_mixing_sweep.md)
- The module README, with the full derivation:
  [`src/core/algorithms/gdpso/README.md`](https://github.com/csilab-ufop/pymocd/blob/master/src/core/algorithms/gdpso/README.md)
