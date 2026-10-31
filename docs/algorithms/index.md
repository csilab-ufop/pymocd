---
description: Overview of the nine community-detection algorithms in pymocd, the library's own RIMPSO and HP-MOCD and seven re-implemented baselines, comparing their objectives, search engines, selection rules, Pareto-front access, determinism and original code.
---

# Algorithms

pymocd exposes **nine community-detection algorithms** through **ten detector
entry points**: Shi-MOCD ships under two selection rules,
[`mocd_q`][pymocd.mocd_q] and [`mocd_d`][pymocd.mocd_d]. Two of the nine are
this library's own methods, **RIMPSO** and **HP-MOCD**. The other seven are
**re-implemented baselines**: other authors' published methods, ported to this
repository from their papers. Every detector returns one crisp partition as a
`dict` from node to community, and isolated nodes are assigned community `-1`.
This page compares the nine side by side and links to one page per algorithm.
For a decision guide, see [Choosing a detector](../guide/choosing.md).

## Overview

The first table lists what each algorithm optimizes and how it searches.

| Algorithm | Entry points | Paper | Objectives | Search engine |
|---|---|---|---|---|
| [RIMPSO](rimpso.md) | [`rimpso`][pymocd.rimpso] | Santos, in preparation (2026) | `cut` and `pair`, both minimized: an exact split of the Constant Potts Model (CPM) | memetic particle swarm, one particle per rung of a geometric resolution ladder |
| [HP-MOCD](hpmocd.md) | [`hpmocd`][pymocd.hpmocd], [`HpMocd`][pymocd.HpMocd] | [Santos et al. (2025)](https://doi.org/10.1007/s13278-025-01519-7), *Social Network Analysis and Mining* | `intra` and `inter`, both minimized: a split of Newman modularity | parallel NSGA-II |
| [CDRME](cdrme.md) | [`cdrme`][pymocd.cdrme] | [Dabaghi-Zarandi et al. (2025)](https://doi.org/10.1016/j.jnca.2024.104070), *Journal of Network and Computer Applications* | one maximized scalar, the paper's Eq. (12): inner linkage plus outer linkage | softmax-weighted random walks, then stochastic agglomerative merge chains and a similarity-driven mutation |
| [MMCoMO](mmcomo.md) | [`mmcomo`][pymocd.mmcomo] | [Zhang et al. (2023)](https://doi.org/10.1109/MCI.2023.3277773), *IEEE Computational Intelligence Magazine* | kernel *k*-means (KKM) and ratio cut (RC), both minimized | macro/micro co-evolutionary NSGA-II over a dense diffusion kernel |
| [NSGA-III-CCM](ccm.md) | [`ccm`][pymocd.ccm] | [Shaik et al. (2021)](https://doi.org/10.1007/s42979-020-00382-x), *SN Computer Science* | community score, community fitness and modularity, all maximized | NSGA-III |
| [NSGA-III-KRM](krm.md) | [`krm`][pymocd.krm] | [Shaik et al. (2021)](https://doi.org/10.1007/s42979-020-00382-x), *SN Computer Science* | KKM and RC minimized, modularity maximized | NSGA-III |
| [GDPSO](gdpso.md) | [`gdpso`][pymocd.gdpso] | [Cai et al. (2015)](https://doi.org/10.1016/j.ins.2014.09.041), *Information Sciences* | Newman–Girvan modularity, maximized (single objective) | greedy discrete particle swarm |
| [Shi-MOCD](shi-mocd.md) | [`mocd_q`][pymocd.mocd_q], [`mocd_d`][pymocd.mocd_d] | [Shi et al. (2012)](https://doi.org/10.1016/j.asoc.2011.10.005), *Applied Soft Computing* | `intra` and `inter`, both minimized: a split of Newman modularity | PESA-II |
| [MOGA-Net](moga-net.md) | [`moga_net`][pymocd.moga_net] | [Pizzuti (2012)](https://doi.org/10.1109/TEVC.2011.2161090), *IEEE Transactions on Evolutionary Computation* | community score and community fitness, both maximized | genetic algorithm with NSGA-II ranking and elitist generational replacement |

The second table lists how each algorithm turns its search into one partition,
whether the candidate set is exposed, whether repeated calls agree, and whether
the original authors released code. The *Deterministic* column follows the
*Determinism* section of each module README: "yes" means that the same graph,
with its nodes and edges in the same order, and the same parameters give the
same partition on every call, at any thread count. Listing the same nodes or
edges in a different order can change the partition, so keep the input order
fixed when results must be reproduced (see
[Input order](../guide/performance.md#input-order)).

| Algorithm | Selection rule | Front accessor | Deterministic | Original code |
|---|---|---|---|---|
| [RIMPSO](rimpso.md) | the member that a degree-corrected assortative block model fits best, with a Schwarz penalty; no labels, no parameter | [`rimpso_fronts`][pymocd.rimpso_fronts]; [`rimpso_select`][pymocd.rimpso_select] applies the rule to external candidates | yes, with a `seed` argument | this library |
| [HP-MOCD](hpmocd.md) | maximum modularity $Q$ on the rank-1 front | [`hpmocd_fronts`][pymocd.hpmocd_fronts]; [`HpMocd.generate_pareto_front`][pymocd.HpMocd.generate_pareto_front] | no | this library |
| [CDRME](cdrme.md) | maximum modularity $Q$; ties broken by Eq. (12) | none (single objective) | yes; no `seed` argument | a private notebook supplied by the authors, vendored in this repository at [`res/original_algs/cdrme`](https://github.com/csilab-ufop/pymocd/tree/master/res/original_algs/cdrme) |
| [MMCoMO](mmcomo.md) | maximum modularity $Q$ on the merged rank-1 front | [`mmcomo_fronts`][pymocd.mmcomo_fronts] | no | none published |
| [NSGA-III-CCM](ccm.md) | maximum modularity $Q$ on the rank-1 front | [`ccm_fronts`][pymocd.ccm_fronts] | no | none published |
| [NSGA-III-KRM](krm.md) | maximum modularity $Q$ on the rank-1 front | [`krm_fronts`][pymocd.krm_fronts] | no | none published |
| [GDPSO](gdpso.md) | the best position the swarm ever held | none (single objective) | yes; no `seed` argument | [doctor-cai/GDPSO](https://github.com/doctor-cai/GDPSO), C++, no license file |
| [Shi-MOCD](shi-mocd.md) | `mocd_q`: maximum modularity $Q$ (Shi Eq. 3.8); `mocd_d`: the member whose minimum distance to Erdős–Rényi control fronts is largest (Shi Eqs. 3.9–3.11) | none exposed | no | none published |
| [MOGA-Net](moga-net.md) | maximum modularity $Q$ on the rank-1 front (Pizzuti Sec. V-E) | [`moga_net_fronts`][pymocd.moga_net_fronts] | no | [Moganet2016.zip](https://staff.icar.cnr.it/pizzuti/codice/Moganet2016.zip), MATLAB shipped as obfuscated P-code |

Each algorithm also has a module README next to its source code, with the full
derivation, the parameter table and every deliberate divergence from its paper:
[`rimpso`](https://github.com/csilab-ufop/pymocd/blob/master/src/core/algorithms/rimpso/README.md) ·
[`hpmocd`](https://github.com/csilab-ufop/pymocd/blob/master/src/core/algorithms/hpmocd/README.md) ·
[`cdrme`](https://github.com/csilab-ufop/pymocd/blob/master/src/core/algorithms/cdrme/README.md) ·
[`mmcomo`](https://github.com/csilab-ufop/pymocd/blob/master/src/core/algorithms/mmcomo/README.md) ·
[`ccm`](https://github.com/csilab-ufop/pymocd/blob/master/src/core/algorithms/ccm/README.md) ·
[`krm`](https://github.com/csilab-ufop/pymocd/blob/master/src/core/algorithms/krm/README.md) ·
[`gdpso`](https://github.com/csilab-ufop/pymocd/blob/master/src/core/algorithms/gdpso/README.md) ·
[`mocd`](https://github.com/csilab-ufop/pymocd/blob/master/src/core/algorithms/mocd/README.md) (Shi-MOCD) ·
[`moganet`](https://github.com/csilab-ufop/pymocd/blob/master/src/core/algorithms/moganet/README.md) ·
[index](https://github.com/csilab-ufop/pymocd/blob/master/src/core/algorithms/README.md).

## Which one should I use?

[`rimpso`][pymocd.rimpso] is the recommended default. One run covers every
resolution, and the partition is chosen with no ground truth and no resolution
parameter to set. [`hpmocd`][pymocd.hpmocd] gives HP-MOCD with fixed defaults
and maximum-modularity selection, and the
[`HpMocd`][pymocd.HpMocd] class accepts your own objective functions. The seven
re-implemented baselines are meant for reproducing papers and for benchmarks.
[Choosing a detector](../guide/choosing.md) expands on each case.

## Library methods

RIMPSO and HP-MOCD were developed in this library, and this repository holds
their reference implementations. Their full descriptions are on their own
pages: [RIMPSO](rimpso.md) and [HP-MOCD](hpmocd.md).

### RIMPSO

RIMPSO (Resolution-Indexed Memetic Particle Swarm Optimisation) minimizes two
objectives over a partition $\mathcal{C}$ of a graph with $n$ nodes and $m$
edges, where $e_c$ is the number of internal edges and $n_c$ the number of
nodes of community $c$:

$$
\operatorname{cut}(\mathcal{C}) = 1 - \frac{1}{m} \sum_{c} e_c,
\qquad
\operatorname{pair}(\mathcal{C}) = \binom{n}{2}^{-1} \sum_{c} \binom{n_c}{2}.
$$

Maximizing the Constant Potts Model at any resolution $\gamma$ is the same as
minimizing a positive weighted sum of these two terms,
$\operatorname{cut} + (\gamma/\gamma_d)\,\operatorname{pair}$ with
$\gamma_d = m/\binom{n}{2}$, so the Pareto front over
$(\operatorname{cut}, \operatorname{pair})$ is the graph's resolution profile,
and $\gamma$ is not a parameter the caller has to set. A particle swarm builds that front, with each
particle pinned to one rung of a geometric resolution ladder. Every
`ls_period` iterations a particle is repaired back to a CPM local optimum and a
community-merge sweep runs. A label-free rule then selects one member: the
partition that a degree-corrected assortative block model fits best once its
free densities are penalized.

[`rimpso`][pymocd.rimpso] returns the selected partition,
[`rimpso_fronts`][pymocd.rimpso_fronts] returns the whole front with each
member's $(\operatorname{cut}, \operatorname{pair})$ point and the selected
index, and [`rimpso_select`][pymocd.rimpso_select] runs the selection rule
alone on partitions produced elsewhere. RIMPSO is deterministic: the same graph
(same nodes and edges, listed in the same order) and the same parameters,
`seed` included, give the same partition on any number of threads; listing the
nodes or edges in another order can change it. See the
[RIMPSO page](rimpso.md) for the search, the selector and the parameters.

### HP-MOCD

HP-MOCD minimizes the pair $(\mathrm{intra}, \mathrm{inter})$ that its paper
takes from [Shi et al. (2010)](https://doi.org/10.1109/CEC.2010.5585987), an
exact decomposition of Newman modularity,
$Q = 1 - \mathrm{intra} - \mathrm{inter}$, with a parallel NSGA-II, and
returns the member of the rank-1 front with the highest $Q$. It was published
by [Santos et al. (2025)](https://doi.org/10.1007/s13278-025-01519-7).

[`hpmocd`][pymocd.hpmocd] and [`hpmocd_fronts`][pymocd.hpmocd_fronts] take the
graph and nothing else. They run with fixed defaults (`pop_size=100`,
`num_gens=100`, `cross_rate=0.7`, `mut_rate=0.5`). The population and
generation counts match the paper's experiments, but the paper used a
crossover probability of 0.8 and a mutation probability of 0.2
(Santos et al., 2025, Sec. 5.1 and Table 1);
`HpMocd(G, cross_rate=0.8, mut_rate=0.2).run()` runs that setting. The
[HP-MOCD page](hpmocd.md#parameters) compares the defaults with the paper's
setting. The [`HpMocd`][pymocd.HpMocd] class exposes the same search with a tunable budget,
a per-generation callback, the objective values of every front member, and
user-supplied Python objective functions. HP-MOCD is not deterministic and has
no seed parameter. See the [HP-MOCD page](hpmocd.md) for details.

<a id="the-re-implemented-baselines"></a>

## Re-implemented baselines

The seven detectors below are other authors' algorithms, ported from their
papers so that benchmarks in this repository compare against a faithful
implementation rather than a paraphrase. Each module README lists what was
read, what the paper leaves open, and every place where the port deliberately
diverges from it.

- [**CDRME**](cdrme.md) (Dabaghi-Zarandi et al., 2025): random walks build
  primary communities, and agglomerative merge chains optimize one linkage
  scalar.
- [**MMCoMO**](mmcomo.md) (Zhang et al., 2023): a macro population of
  community-center genomes and a micro population of label vectors co-evolve
  under kernel *k*-means and ratio cut.
- [**NSGA-III-CCM**](ccm.md) (Shaik et al., 2021): NSGA-III over community
  score, community fitness and modularity.
- [**NSGA-III-KRM**](krm.md) (Shaik et al., 2021): NSGA-III over kernel
  *k*-means, ratio cut and modularity.
- [**GDPSO**](gdpso.md) (Cai et al., 2015): a greedy discrete particle swarm
  that maximizes modularity.
- [**Shi-MOCD**](shi-mocd.md) (Shi et al., 2012): PESA-II over the
  $(\mathrm{intra}, \mathrm{inter})$ split of modularity, under two published
  selection rules (`mocd_q`, `mocd_d`).
- [**MOGA-Net**](moga-net.md) (Pizzuti, 2012): a genetic algorithm over
  community score and community fitness.

### Original code

Reference code exists for three of the seven, according to their module
READMEs:

- **MOGA-Net**: [Moganet2016.zip](https://staff.icar.cnr.it/pizzuti/codice/Moganet2016.zip)
  is MATLAB, shipped as obfuscated P-code. The search loop cannot be read, but
  its result files can, and they were used to settle the objectives.
- **GDPSO**: [doctor-cai/GDPSO](https://github.com/doctor-cai/GDPSO) is C++
  with no license file, so this port was written clean-room from a description
  of it; no reference source was read while writing it.
- **CDRME**: the reference is a private Python notebook supplied by the
  authors. No public repository exists, so it is vendored in this repository at
  [`res/original_algs/cdrme`](https://github.com/csilab-ufop/pymocd/tree/master/res/original_algs/cdrme),
  converted to plain Python. It was consulted only to disambiguate what the
  paper leaves open, and four of its choices were deliberately not reproduced.

No original code was published for MMCoMO, NSGA-III-CCM, NSGA-III-KRM or
Shi-MOCD. For MMCoMO, NSGA-III-KRM and Shi-MOCD, the module READMEs state that
the authors' sites and GitHub were searched, so only the papers' text and
tables constrain those ports.

### Fronts and single-objective baselines

Two of the seven are single-objective, so they have no Pareto front and no
`*_fronts` accessor: [`gdpso`][pymocd.gdpso] maximizes Newman–Girvan
modularity directly, and [`cdrme`][pymocd.cdrme] maximizes the sum that its
paper's Eq. (12) forms from inner and outer linkage. [`mocd_q`][pymocd.mocd_q]
and [`mocd_d`][pymocd.mocd_d] are multi-objective but expose no front accessor
either; the two names run the same PESA-II search under Shi's two published
model-selection rules.

!!! warning "Shi-MOCD defaults are not Shi's"
    [`mocd_q`][pymocd.mocd_q] and [`mocd_d`][pymocd.mocd_d] default to this
    repository's HP-MOCD-parity benchmark budget (`cross_rate=0.9`,
    `mut_rate=0.1`), not to the paper's `pc=0.6`, `pm=0.4` with per-graph
    population and generation counts from its Table 1. Pass those explicitly
    to reproduce the paper.

## Deprecated aliases

Earlier versions of pymocd exposed RIMPSO under the function names
`mr_mocd*`, and before that `scale*`. Those names remain as aliases so that code pinned to an
older API keeps working:

| Deprecated name | Same object as |
|---|---|
| `pymocd.mr_mocd` | [`rimpso`][pymocd.rimpso] |
| `pymocd.mr_mocd_fronts` | [`rimpso_fronts`][pymocd.rimpso_fronts] |
| `pymocd.mr_mocd_select` | [`rimpso_select`][pymocd.rimpso_select] |
| `pymocd.scale` | [`rimpso`][pymocd.rimpso] |
| `pymocd.scale_fronts` | [`rimpso_fronts`][pymocd.rimpso_fronts] |

Each alias is the same function object under a former name. The aliases emit
no warning and do not appear in the type stubs. New code should use the
`rimpso*` names. The check below runs against the installed package when this
page is built:

```python exec="on" source="above" session="algorithms" result="text"
import pymocd

aliases = {
    "mr_mocd": "rimpso",
    "mr_mocd_fronts": "rimpso_fronts",
    "mr_mocd_select": "rimpso_select",
    "scale": "rimpso",
    "scale_fronts": "rimpso_fronts",
}
for old, new in aliases.items():
    print(f"pymocd.{old:<15} is pymocd.{new:<14} {getattr(pymocd, old) is getattr(pymocd, new)}")
```

## Citation

To cite pymocd, see [Citing pymocd](../citing.md). Each algorithm page has a
*Reference* section that gives the full citation of the paper defining the
method, with a BibTeX entry. RIMPSO is the exception: its paper is in
preparation and has no DOI yet.
