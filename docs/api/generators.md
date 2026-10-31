# Benchmark graphs

Generators of graphs with a planted community structure, returned with their
ground truth so a detector's output can be scored against it. The generator
needs no other package: it returns a plain edge list and a
`dict[int, int]`, which [`nmi`][pymocd.nmi], [`ami`][pymocd.ami] and the other
metrics accept as ground truth. See [Benchmark graphs](../guide/benchmarks.md)
for what the generator guarantees and how it differs from NetworkX's.

## LFR

::: pymocd.lfr_graph
