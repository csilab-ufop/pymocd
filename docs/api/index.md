# API reference

Everything below lives in the top-level `pymocd` namespace. The reference is
generated from `pymocd.pyi`, the type stub that `cargo run --bin stub_gen`
writes from the Rust sources, so every signature and default here is the one the
installed package exposes. The docstrings themselves live next to the code, in
[`src/api/`](https://github.com/csilab-ufop/pymocd/tree/master/src/api) and
[`src/core/algorithms/hpmocd/api.rs`](https://github.com/csilab-ufop/pymocd/blob/master/src/core/algorithms/hpmocd/api.rs).

Two conventions hold throughout:

- **Graphs** are NetworkX graphs with integer node ids, or igraph graphs
  (see [Input graphs](../guide/input-graphs.md)).
- **Partitions** are plain `dict[int, int]` mapping node to community id.
  Isolated nodes are assigned community `-1`.

<!-- pymocd:api-summary -->

## Deprecated aliases

`pymocd.mr_mocd`, `pymocd.mr_mocd_fronts` and `pymocd.mr_mocd_select` are the
names RIMPSO carried before it was renamed; `pymocd.scale` and
`pymocd.scale_fronts` are older still. All five are the same function objects as
[`rimpso`][pymocd.rimpso], [`rimpso_fronts`][pymocd.rimpso_fronts] and
[`rimpso_select`][pymocd.rimpso_select], kept so pinned callers keep working.
They emit no warning and are not in the type stub. Use the new names.
