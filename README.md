<div align="center">
  <img src="https://raw.githubusercontent.com/csilab-ufop/pymocd/master/res/logo.svg" alt="pymocd logo" width="75%">
</div>

<div align="center">

[![PyPI Publish](https://github.com/csilab-ufop/pymocd/actions/workflows/release.yml/badge.svg)](https://github.com/csilab-ufop/pymocd/actions/workflows/release.yml)![Rust Compilation](https://img.shields.io/github/actions/workflow/status/csilab-ufop/pymocd/rust.yml)
![PyPI - Version](https://img.shields.io/pypi/v/pymocd)
![PyPI - License](https://img.shields.io/pypi/l/pymocd)

</div>

**pymocd**: multi-objective evolutionary community detection for NetworkX and
igraph graphs, with a parallel Rust core.

```bash
pip install pymocd
```

**Documentation: <https://csilab.ufop.br/pymocd/>**

- [Getting started](https://csilab.ufop.br/pymocd/getting-started/) and
  [installation](https://csilab.ufop.br/pymocd/guide/installation/)
- [Algorithms](https://csilab.ufop.br/pymocd/algorithms/): the library's
  methods and the re-implemented baselines, with their papers
- [Evaluating partitions](https://csilab.ufop.br/pymocd/guide/evaluation/)
  (NMI, AMI, ARI, F1) and
  [LFR benchmark graphs](https://csilab.ufop.br/pymocd/guide/benchmarks/)
- [Examples](https://csilab.ufop.br/pymocd/examples/) and the
  [API reference](https://csilab.ufop.br/pymocd/api/)
- [Changelog](https://csilab.ufop.br/pymocd/changelog/),
  [contributing](https://csilab.ufop.br/pymocd/contributing/) and
  [citing pymocd](https://csilab.ufop.br/pymocd/citing/)

### Citation

If you use pymocd in your research, please cite:

> Santos, G. O., Vieira, L. S., Rossetti, G., Ferreira, C. H. G., & Moreira,
> G. J. P. (2025). A high-performance evolutionary multiobjective community
> detection algorithm. *Social Network Analysis and Mining*, 15(1), 110.
> <https://doi.org/10.1007/s13278-025-01519-7>

```bibtex
@article{Santos2025,
  author    = {Santos, Guilherme O. and Vieira, Lucas S. and Rossetti, Giulio and Ferreira, Carlos H. G. and Moreira, Gladston J. P.},
  title     = {A high-performance evolutionary multiobjective community detection algorithm},
  journal   = {Social Network Analysis and Mining},
  year      = {2025},
  volume    = {15},
  number    = {1},
  pages     = {110},
  doi       = {10.1007/s13278-025-01519-7},
  url       = {https://doi.org/10.1007/s13278-025-01519-7},
  issn      = {1869-5469},
  date      = {2025-11-18}
}
```

For the re-implemented baselines, also cite their original papers; see
[Citing pymocd](https://csilab.ufop.br/pymocd/citing/).

pymocd is licensed under the GPL-3.0-or-later.
