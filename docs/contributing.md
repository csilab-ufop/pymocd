---
description: How to report issues, set up a development build of pymocd, run the tests as CI does, add a detector, and write documentation and gallery examples.
---

# Contributing

Contributions to pymocd are welcome: bug reports, fixes, new detectors,
documentation and examples. This page describes how the repository is set up,
so that a change passes the same checks locally that it meets in continuous
integration (CI). It covers reporting issues, building pymocd from source,
running the tests, how the API reference is generated from the Rust sources,
the steps for adding a detector, writing documentation, and the commit
conventions of the project.

## Reporting issues

Issues are tracked on
[GitHub](https://github.com/csilab-ufop/pymocd/issues). Opening a new issue
offers the three templates in `.github/ISSUE_TEMPLATE/`:

| Template | Its description | It asks for |
|---|---|---|
| **Bug report** | "Create a report to help us improve" | A description of the bug, steps to reproduce it, the expected behavior, screenshots if applicable, your operating system, and any other context |
| **Feature request** | "Suggest an idea for this project" | The problem it solves, the solution you would like, alternatives you considered, and any other context |
| **Custom issue template** | None (placeholder text) | Free text: the template body is empty |

For a bug in a detector, a minimal script is the most useful reproduction: the
graph (ideally one networkx can build, such as `nx.karate_club_graph()`), the
exact call with its keyword arguments, and the full traceback or the
unexpected output. Include your pymocd and Python versions;
`importlib.metadata.version("pymocd")` returns the former.

## Development setup

pymocd is a Rust crate compiled into a Python extension module with
[PyO3](https://pyo3.rs) and [maturin](https://www.maturin.rs/). To build it
you need:

- a **Rust toolchain**. The crate uses Rust edition 2024 (`Cargo.toml`). The
  workflows that build the Python package (`python.yml`, `docs.yml` and
  `release.yml`) install the stable toolchain with
  `dtolnay/rust-toolchain@stable`; `rust.yml` has no toolchain step;
- **Python 3.10 or newer** (`requires-python = ">=3.10"` in `pyproject.toml`).
  The Makefile route below needs **Python 3.11 or newer** as `python3`: some
  of the versions pinned in `res/requirements.txt` (networkx 3.6.1,
  numpy 2.4.2, pandas 3.0.1, scikit-learn 1.8.0) declare
  `Requires-Python >=3.11` on PyPI.

Clone the repository and let the Makefile create a development environment:

```bash
git clone https://github.com/csilab-ufop/pymocd
cd pymocd
make build
```

The Makefile targets used during development:

| Target | What it does |
|---|---|
| `make dependencies` | Creates a virtual environment in `.venv` with `python3 -m venv`, upgrades its pip, and installs the pinned development and benchmark requirements of `res/requirements.txt` into it (maturin, networkx, igraph, matplotlib and others). |
| `make stubs` | Runs `cargo run --bin stub_gen`, which writes the type stub `pymocd.pyi` to the repository root. |
| `make build` | `make dependencies`, then generates the stub, installs pymocd into `.venv` with `maturin develop --release`, and deletes `pymocd.pyi` again. |
| `make test` | `make build`, then the Rust test suite with `cargo test`. |
| `make all` | `make build` and `make test`. |
| `make docs` / `make docs-serve` | Build or serve this documentation site (see [Writing documentation](#writing-documentation)). |
| `make clean` | Runs `cargo clean` and removes `.venv`, `target`, `build`, `dist`, `*.egg-info`, `pymocd.pyi` and every `__pycache__` directory. |

Without the Makefile, build with `cargo run --bin stub_gen` followed by
`pip install .`; CI runs the same two steps with uv
(`uv pip install --system .`). See
[Building from source](guide/installation.md#building-from-source).

!!! warning "Leave the `extension-module` feature off for cargo"
    `Cargo.toml` deliberately keeps PyO3's `extension-module` feature out of
    the default features. Its comment explains why: the feature omits the
    link against libpython, so with it on by default every `cargo build`,
    `cargo test` and `cargo run` fails to link the `stub_gen` binary. maturin
    turns it on when it builds the extension (`[tool.maturin] features` in
    `pyproject.toml`).

The `benchmark*` and `hard-*` targets drive the benchmark harness in
`benchmarks/`. They are not part of CI.

## Running the tests

There are two test suites. Rust unit tests live in `#[cfg(test)]` modules next
to the code they test: in the detector modules under `src/core/algorithms/`
and in the graph and metrics code under `src/core/`. The Python suite in
`tests/` exercises the installed package:

| File | What it checks |
|---|---|
| `tests/test_partitions.py` | The detectors (all but `gdpso` and `cdrme`) and the `rimpso_fronts` and `mmcomo_fronts` accessors return valid partitions of the karate club graph (every node assigned, integer community ids), and `hpmocd` and `rimpso` recover both triangles of a graph made of two triangles joined by one edge. |
| `tests/test_determinism.py` | `rimpso` and `rimpso_fronts` return identical results on repeated calls, `rimpso` returns the same partition under 1 and 4 threads, and repeated `hpmocd` runs return valid partitions. |
| `tests/test_sanity.py` | A modularity floor on karate for `rimpso` and `hpmocd`, perfect self-agreement of the metrics, community `-1` for isolated nodes, `rimpso_fronts` objective values that match their partitions, and the complete-graph and empty-graph edge cases. |
| `tests/test_igraph.py` | Every detector accepts an `igraph.Graph` (Zachary's karate club) and assigns every vertex, `rimpso` returns the same partition for an igraph graph and a NetworkX copy with the same node and edge order, and an isolated vertex gets community `-1`. Skipped when python-igraph is not installed. |

Two workflows run them on pull requests:

| Workflow | Runs on | Steps |
|---|---|---|
| `.github/workflows/rust.yml` | pull requests to `master` | `cargo build --verbose`, then `cargo test --verbose` |
| `.github/workflows/python.yml` | every pull request | Python 3.10; `cargo run --bin stub_gen`; installs networkx, python-igraph and pymocd from the checkout; runs the Python suite from `tests/` |

To run the same checks locally, run the Rust suite from the repository root:

```bash
cargo test --verbose
```

and the Python suite from `tests/`, with a Python interpreter that has pymocd
and networkx installed (after `make build`, the one in `.venv`):

```bash
cd tests
python -m unittest discover -q
```

`test_determinism.py` starts child processes to compare thread counts,
because [`max_cores`][pymocd.max_cores] sizes the thread pool only once per
process.

## Docstrings and the type stub

The [API reference](api/index.md) is generated; no signature on it is written
by hand. The pipeline is:

1. **Docstrings are Rust doc comments.** Each Python-facing function, and the
   `HpMocd` class with its methods, carries a `///` comment, which
   mkdocstrings parses as a Google-style docstring
   (`docstring_style: google` in `mkdocs.yml`). Most detector docstrings
   have `Args:` and `Returns:` sections; the metric functions and
   `max_cores` have short summaries only. They live in
   `src/api/detectors.rs`, `src/api/metrics.rs`,
   `src/core/utils/threading.rs` (`max_cores`) and
   `src/core/algorithms/hpmocd/api.rs` (the `HpMocd` class). The first
   sentence of a docstring becomes its summary on the API index page.
2. **Annotations select what is stubbed.** `#[gen_stub_pyfunction]` (and
   `#[gen_stub_pyclass]`, `#[gen_stub_pymethods]` for `HpMocd`) make
   [pyo3-stub-gen](https://github.com/Jij-Inc/pyo3-stub-gen) include an item.
   Parameter defaults come from the `#[pyo3(signature = (...))]` attribute.
   Some functions that return a generic Python object carry a
   `#[gen_stub(override_return_type(...))]` attribute that states the precise
   return type; without one, the stub declares `typing.Any`, as it does for
   [`rimpso_select`][pymocd.rimpso_select].
3. **`cargo run --bin stub_gen` writes `pymocd.pyi`.** The binary
   (`src/bin/stub_gen.rs`) calls the gatherer that
   `define_stub_info_gatherer!` defines in `src/lib.rs`. The stub is not
   tracked in git (`.gitignore`); maturin packages it, together with
   `py.typed`, into every wheel and source distribution.
4. **The docs build reads the stub.** `docs/hooks.py` copies `pymocd.pyi` to
   `.docs-stub/pymocd.py`, rewriting PyO3's `__new__` constructors as
   `__init__`, and mkdocstrings renders that copy. The build stops if
   `pymocd.pyi` is missing (run `make stubs`), if a public name in the stub
   has no `::: pymocd.<name>` directive in `docs/api/detectors.md`,
   `docs/api/fronts.md` or `docs/api/metrics.md`, or if a directive names
   something the stub does not contain.

To change a signature, default or docstring, edit the Rust source and
regenerate the stub. Never edit `pymocd.pyi` by hand: the next
`cargo run --bin stub_gen` overwrites it.

## Adding a detector

Each detector is a self-contained module under `src/core/algorithms/`: it owns
its representation, objectives, operators, search engine and entry point. The
steps below follow the existing modules; the
[`mocd_q`][pymocd.mocd_q] entry point serves as the example.

1. **Write the module.** Create `src/core/algorithms/<name>/` and declare it
   with `pub mod <name>;` in `src/core/algorithms/mod.rs`. The existing
   modules each have a `mod.rs`, an `api.rs` with the entry point, and their
   shipped defaults as named constants (`defaults.rs` or `config/`).
2. **Expose it to Python.** Add a function to `src/api/detectors.rs` with a
   Google-style docstring and the stub annotation. An abridged excerpt of the
   existing `mocd_q` entry point:

    ```rust
    /// Run Shi-MOCD (Shi, Yan, Cai, Wu 2012) — PESA-II over Shi's
    /// ...
    /// Returns:
    ///     ``dict[node, community]``. Isolated nodes get community ``-1``.
    #[gen_stub_pyfunction]
    #[pyfunction]
    #[pyo3(name = "mocd_q", signature = (graph, pop_size = mocd::DEFAULT_POP_SIZE, num_gens = mocd::DEFAULT_NUM_GENS, cross_rate = mocd::BENCH_CROSS_RATE, mut_rate = mocd::BENCH_MUT_RATE))]
    pub fn mocd_q_fn(
        graph: &Bound<'_, PyAny>,
        pop_size: usize,
        num_gens: usize,
        cross_rate: f64,
        mut_rate: f64,
    ) -> PyResult<Partition> {
    ```

    Every detector returns one partition as `dict[node, community]`, with
    isolated nodes in community `-1`. A front accessor, if the detector has
    one, is a second function named `<name>_fronts`.
3. **Register it** in the `pymocd` module function in `src/lib.rs`, next to
   the other detectors:

    ```rust
    m.add_function(wrap_pyfunction!(mocd_q_fn, m)?)?;
    ```

4. **Regenerate the stub** with `cargo run --bin stub_gen` and check that the
   new name appears in `pymocd.pyi`.
5. **Add it to the API reference.** Put a `::: pymocd.<name>` directive in
   `docs/api/detectors.md` and, if the detector has a front accessor, one
   for `<name>_fronts` in `docs/api/fronts.md`. The docs build fails while a
   public name is undocumented.
6. **Write its algorithm page**, `docs/algorithms/<name>.md`, and add it to
   `nav` in `mkdocs.yml` under *Library methods* or *Re-implemented
   baselines*. Follow the structure of the existing pages: the paper with full
   citation and DOI, objectives, search, selection, parameters, complexity,
   when to use it, and divergences from the paper. The parameter table must
   be marked with `<!-- pymocd:defaults pymocd.<name> -->`; the build checks
   it against the stub (see [Writing documentation](#writing-documentation)).
   Add the method to the [Algorithms overview](algorithms/index.md) too.
7. **Write the module README**, `src/core/algorithms/<name>/README.md`. Like
   the existing ones, as `src/core/algorithms/README.md` describes them, it
   carries the full citation, the objective equations, the representation,
   the search loop, the parameter table, whether the detector is
   bit-deterministic, and every known divergence from its paper. Add a row
   for it to the table in `src/core/algorithms/README.md`. The top-level
   `README.md` only points to the documentation, so it needs no change.
8. **Test it.** Add a case to `tests/test_partitions.py` (a valid partition of
   the karate club graph, plus the front accessor if there is one), add the
   entry point to the `DETECTORS` list in `tests/test_igraph.py` so that it is
   also tested on igraph input, and add Rust unit tests for the module's
   pieces.
9. **Cite it.** For a re-implemented baseline, add its paper, DOI and BibTeX
   to [Citing pymocd](citing.md) and [Publications](publications.md).

To include the detector in the benchmark harness, register a wrapper with the
`@algorithm(...)` decorator in `benchmarks/_exp_synt_net/algorithms.py`, next
to the existing wrappers. That file wraps every detector except `gdpso` and
`cdrme`.

## Writing documentation

The site is built with [MkDocs](https://www.mkdocs.org/) and the
[Material theme](https://squidfunk.github.io/mkdocs-material/). Its
configuration is `mkdocs.yml`, its pages are under `docs/`, and its build
dependencies are listed in `docs/requirements.txt`. Add any new docs
dependency there, because both the Makefile and the deploy workflow install
from it.

```bash
make docs        # build the site into site/ with mkdocs build --strict
make docs-serve  # serve it locally with live reload
```

Both targets first build pymocd into `.venv` and regenerate the stub, because
the code on the pages runs during the build. `--strict` turns every warning
into an error, including a broken internal link, a missing anchor and a page
that is not listed in `nav`. The site is deployed to GitHub Pages by
`.github/workflows/docs.yml` on every push to `master`; pull requests do not
build it, so run `make docs` before opening one that touches `docs/`,
`examples/` or a docstring.

### Page conventions

- **Front matter.** The guide, algorithm and research pages start with YAML
  front matter holding a `description:` line, then one `#` title.
- **Executed code.** Python examples on the pages are written as
  [markdown-exec](https://pawamoy.github.io/markdown-exec/) blocks with
  `exec="on"`, so they run at build time and their output is that of the
  current package. Blocks that share a `session` share their variables, in
  page order:

    ````markdown
    ```python exec="on" source="above" session="<page-slug>" result="text"
    import networkx as nx
    import pymocd

    G = nx.karate_club_graph()
    print(len(set(pymocd.rimpso(G).values())))
    ```
    ````

    Keep these blocks fast (small graphs, small budgets), because the whole
    site is rebuilt in CI. Shell commands go in plain `bash` blocks, which are
    not executed.

- **API links.** Refer to the API with cross-references such as
  `` [`rimpso`][pymocd.rimpso] ``; any name in the stub works as
  `pymocd.<name>`.
- **Parameter tables.** A table placed right after an HTML comment of the form
  `<!-- pymocd:defaults pymocd.<name> -->` is compared with the stub on every
  build. Its first column is the parameter name in backticks and its second
  the default in backticks, written as in the stub or as an exact fraction
  such as `1/68`. Every parameter except `graph` must appear, and no other
  row may. Several targets can share one table if their shared parameters
  agree.
- **Missing facts.** Do not guess. When a fact is not known,
  write a sentence that is true without it and put an HTML comment starting
  with `TODO(docs):` right after it, saying what is missing and where it could
  come from. The build lists every such marker at the end.
- **Moved pages.** Keep URLs stable. When a page moves, add the old path to
  `redirect_maps` under the `redirects` plugin in `mkdocs.yml`.

### Gallery examples

The [Examples](examples/index.md) gallery is generated by
[mkdocs-gallery](https://smarie.github.io/mkdocs-gallery/) from the scripts in
`examples/`. Every file whose name starts with `plot_` is executed during the
build, and its printed output and figures are rendered into its page. To add
an example:

- put `plot_<topic>.py` in a section folder of `examples/` (`basics`,
  `benchmarks`, `fronts` or `objectives`); each folder has a `README.md` whose
  heading titles the section. Sections appear in alphabetical order of their
  folder names (the `subsection_order` option of mkdocs-gallery 0.10.4 recurses
  forever, so it is not used), and examples within a section by file name;
- open the script with a module docstring whose first line is a Markdown
  `#` title, followed by a sentence on what the example shows;
- split the script into cells with `# %%` lines; the comment lines right after
  a `# %%` are rendered as Markdown text between the code cells;
- keep it self-contained and quick, and show figures with `plt.show()`.

The generated pages land in `docs/examples/`, which is not tracked in git. The
edit button of a gallery page opens the example's source script.
[Communities in Zachary's karate club](examples/basics/plot_karate_club.md)
is a short example to start from.

## Commits and pull requests

Propose changes as a pull request against `master`, where the CI workflows
above run on them. The repository has no written commit policy; most commit
subjects in its history follow these conventions:

- **A scope, a colon, then a short imperative summary in lower case**, with
  no trailing period: `rimpso: rename mr_mocd and replace the front selector`,
  `tests: move the python suite onto mr_mocd`,
  `ci: use setup-python for the interpreter, uv for installs`. Some subjects
  have no scope (`vendor the cdrme reference implementation`).
- **The scope names what changed**: a detector module (`rimpso`, `cdrme`,
  `gdpso`, `mmcomo`), or an area such as `docs`, `readme`, `tests`, `ci`,
  `bench`, `build(deps)` or `chore`. Older commits also use conventional
  types with a scope in parentheses, such as `refactor(ccm):` and
  `feat(api):`.
- **One logical change per commit.** CDRME, for example, was added as a
  series of commits, one per part (`cdrme: add the topology, objective and
  entry points`, `cdrme: add the configuration and shipped defaults`, and so
  on).
- **A body only when it is needed.** Most commits are a single line; a body
  explains changes the subject cannot, such as a replaced selection rule.

Dependabot checks the Cargo dependencies weekly and opens pull requests for
updates (`.github/dependabot.yml`).

### Releases

`make bump V=<version>` sets the version in `Cargo.toml` and
`pyproject.toml`, commits it as `chore: bump version to <version>`, tags the
commit `v<version>` and pushes `master` with the tag.
The tag triggers `.github/workflows/release.yml`, which builds wheels for
Linux, Windows and macOS on Python 3.10 to 3.13, plus a source distribution,
and publishes them to PyPI.

## Licence

pymocd is licensed under the GNU General Public License, version 3 or later
(GPL-3.0-or-later), as declared by `license` in `Cargo.toml` and by the
licence classifier in `pyproject.toml`; the full text is in
[LICENSE](https://github.com/csilab-ufop/pymocd/blob/master/LICENSE). Rust
source files open with a `//!` summary of the module followed by a licence
header; new files follow the same layout.

## See also

- [Installation](guide/installation.md): installing from PyPI and building
  from source.
- [API reference](api/index.md): the generated reference this pipeline
  produces.
- [Algorithms](algorithms/index.md): the existing detector pages to model a
  new one on.
