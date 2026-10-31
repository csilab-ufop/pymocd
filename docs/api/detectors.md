# Detectors

Every detector takes a NetworkX or igraph graph (see
[Input graphs](../guide/input-graphs.md)) and returns one partition as a
`dict[int, int]` mapping node to community id. Isolated nodes are assigned
community `-1`. The [Algorithms](../algorithms/index.md) pages describe what
each method optimizes, the paper it comes from and when to use it.

<a id="this-librarys-algorithms"></a>

## This library's methods

::: pymocd.rimpso

::: pymocd.hpmocd

The `HpMocd` class runs the same search as [`hpmocd`][pymocd.hpmocd] with the
budget as constructor arguments. It also accepts Python objective functions in
place of the built-in ones and reports progress through a per-generation
callback; see [HP-MOCD](../algorithms/hpmocd.md) and the
[custom objectives example](../examples/objectives/plot_custom_objectives.md).

::: pymocd.HpMocd

## Re-implemented baselines

These entry points port methods published by other authors. Each takes its
paper's parameters, plus any the port adds where the paper is silent; its
[algorithm page](../algorithms/index.md#re-implemented-baselines) lists every
default and says which ones the paper states.

::: pymocd.cdrme

::: pymocd.mmcomo

::: pymocd.ccm

::: pymocd.krm

::: pymocd.gdpso

::: pymocd.mocd_q

::: pymocd.mocd_d

::: pymocd.moga_net

## Deprecated aliases

`pymocd.mr_mocd`, `pymocd.mr_mocd_fronts` and `pymocd.mr_mocd_select` are the
names RIMPSO carried before it was renamed; `pymocd.scale` and
`pymocd.scale_fronts` are older still. All five are the same function objects
as [`rimpso`][pymocd.rimpso], [`rimpso_fronts`][pymocd.rimpso_fronts] and
[`rimpso_select`][pymocd.rimpso_select], kept so pinned callers keep working.
They emit no warning and are not in the type stub. Use the new names.
