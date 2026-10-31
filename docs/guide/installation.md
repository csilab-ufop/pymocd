---
description: Install pymocd from PyPI or build it from source with Rust and maturin, check that the installation works, and fix the problems the repository documents.
---

# Installation

pymocd is distributed on PyPI as prebuilt wheels plus a source distribution.
The package itself has no runtime dependencies, so the graph library you pass
graphs from (networkx) is installed separately. Building from source needs a
Rust toolchain, because the detectors are written in Rust and compiled into a
Python extension module with [PyO3](https://pyo3.rs) and
[maturin](https://www.maturin.rs/).

## Install from PyPI

```bash
pip install pymocd
```

`pyproject.toml` declares no runtime dependencies, so this installs pymocd and
nothing else. The detectors take a NetworkX (or igraph) graph as input, and
the gallery examples draw their figures with matplotlib, so a typical
installation is:

```bash
pip install pymocd networkx matplotlib
```

## Supported Python versions and platforms

`pyproject.toml` sets `requires-python = ">=3.10"`, and its classifiers list
CPython 3.10, 3.11, 3.12 and 3.13 on Linux, macOS and Windows.

The release workflow (`.github/workflows/release.yml`) builds wheels on three
GitHub-hosted runners for four Python versions each, and uploads one source
distribution from the `ubuntu-latest` / Python 3.10 job. It publishes all of
them to PyPI when a tag is pushed or the workflow is started by hand.

| Runner | Python versions |
|---|---|
| `ubuntu-latest` | 3.10, 3.11, 3.12, 3.13 |
| `windows-latest` | 3.10, 3.11, 3.12, 3.13 |
| `macos-latest` | 3.10, 3.11, 3.12, 3.13 |

The workflow does not pass a target architecture to maturin, and the
repository does not list the CPU architectures the published wheels cover. For
6.0.1, PyPI's record of the release (`https://pypi.org/pypi/pymocd/6.0.1/json`)
lists wheels for CPython 3.10–3.13 tagged `manylinux_2_34_x86_64`,
`macosx_11_0_arm64` and `win_amd64`, plus the source distribution; the
[project page on PyPI](https://pypi.org/project/pymocd/) shows the files of
every release.
<!-- TODO(docs): state the CPU architecture (and Linux glibc floor) of each published wheel. release.yml passes no --target to `maturin build`; its "Build universal2 wheels on macOS" comment is not matched by a universal2 target. PyPI's file list for 6.0.1 (checked 2026-10-30) shows macosx_11_0_arm64, manylinux_2_34_x86_64 and win_amd64 wheels for cp310-cp313 plus pymocd-6.0.1.tar.gz; the maintainer should confirm and decide what to promise. -->
When no wheel matches your interpreter or platform, pip downloads the source
distribution and compiles it, which needs a Rust toolchain; see
[Building from source](#building-from-source).

## Verify the installation

Import pymocd and run a detector on a small graph:

```python exec="on" source="above" session="installation" result="text"
from importlib.metadata import version

import networkx as nx
import pymocd

print("pymocd", version("pymocd"))
print("networkx", nx.__version__)

partition = pymocd.rimpso(nx.karate_club_graph())
print("rimpso found", len(set(partition.values())), "communities in the karate club graph")
```

If this prints a version and a community count, the extension module loaded
and ran. The [Getting started](../getting-started.md) page continues from
here.

## Building from source

You need:

- **A Rust toolchain.** The crate uses Rust edition 2024 (`Cargo.toml`), and
  the CI workflows that build the package install the current stable
  toolchain. [rustup](https://rustup.rs/)
  installs and updates it.
- **Python 3.10 or newer** with pip.
- **maturin.** It is the build backend declared in `pyproject.toml`
  (`maturin>=1,<2`), so `pip install .` fetches it into the isolated build
  environment on its own. The Makefile route below installs it into a virtual
  environment instead.

### With pip

The CI workflows `python.yml` and `docs.yml` run the same two build steps,
with uv's pip interface (`uv pip install --system .`) in place of `pip install .`:

```bash
git clone https://github.com/csilab-ufop/pymocd
cd pymocd
cargo run --bin stub_gen   # generates the type stubs, pymocd.pyi
pip install .
```

`pymocd.pyi` is generated from annotations in the Rust sources and is not
tracked in git. `pyproject.toml` packages it, together with `py.typed`, into
the wheel and the source distribution, and the workflows generate it before
they install the package.

### With the Makefile

The Makefile drives a development setup in a local `.venv`:

| Target | What it runs |
|---|---|
| `make build` | Creates `.venv` with `python3 -m venv` and installs the pinned development and benchmark requirements of `res/requirements.txt` into it (this step reruns only when that file changes), runs `cargo run --bin stub_gen`, installs pymocd into `.venv` with `maturin develop --release`, then deletes `pymocd.pyi` from the checkout. |
| `make stubs` | `cargo run --bin stub_gen`, which writes `pymocd.pyi` to the repository root. |
| `make test` | `make build`, then the Rust test suite with `cargo test`. |
| `make docs` / `make docs-serve` | Builds (or serves) this documentation site; see [Contributing](../contributing.md). |

Because `make build` deletes `pymocd.pyi` after installing, run `make stubs`
when you need the stub file in the checkout; the `docs` targets do this
themselves.

The Python test suite lives in `tests/`. The CI workflow `python.yml` installs
networkx, python-igraph and the package, then runs it from that directory:

```bash
cd tests
python -m unittest discover -q
```

## Troubleshooting

### `ModuleNotFoundError: No module named 'networkx'`

pymocd does not install networkx (it declares no runtime dependencies).
Install it alongside pymocd:

```bash
pip install networkx
```

### pip compiles pymocd instead of installing a wheel

If pip downloads `pymocd-<version>.tar.gz` and starts a Rust build, no
published wheel matched your Python version, operating system or CPU
architecture (see the table and the PyPI file list above). Either install
pymocd into a Python version and platform for which a wheel exists, or install
a Rust toolchain so the source build can complete.

### `cargo` fails to link `stub_gen`

`Cargo.toml` deliberately leaves PyO3's `extension-module` feature out of the
default features: that feature omits the link against libpython, and with it
enabled every `cargo build`, `cargo test` and `cargo run` fails to link the
`stub_gen` binary. maturin turns the feature on itself when it builds the
extension (`[tool.maturin] features` in `pyproject.toml`). Do not enable it for
plain cargo commands.

### A detector fails on an `igraph.Graph`

The 6.0.1 release on PyPI cannot read `igraph.Graph` input; this is not an
installation issue, and the fix is not yet released. See
[Input graphs](input-graphs.md#igraph) for the workaround with that release.
