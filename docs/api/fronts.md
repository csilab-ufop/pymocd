# Pareto fronts

These functions return the candidate set a detector selects its partition
from, so you can inspect it or apply your own selection rule. They take the
same keyword arguments as their detector. `rimpso_select` is the exception: it
takes candidate partitions you supply and returns the index RIMPSO's rule
selects together with each candidate's `(cut, pair)` point. See
[Working with Pareto fronts](../guide/fronts.md) for worked examples.

Six detectors have a front accessor. [`gdpso`][pymocd.gdpso] and
[`cdrme`][pymocd.cdrme] optimize a single scalar, so they have no front;
[`mocd_q`][pymocd.mocd_q] and [`mocd_d`][pymocd.mocd_d] are multi-objective but
do not expose theirs.

## RIMPSO

::: pymocd.rimpso_fronts

::: pymocd.rimpso_select

## HP-MOCD

::: pymocd.hpmocd_fronts

## Re-implemented baselines

::: pymocd.mmcomo_fronts

::: pymocd.ccm_fronts

::: pymocd.krm_fronts

::: pymocd.moga_net_fronts
