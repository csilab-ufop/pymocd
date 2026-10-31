---
description: How to cite pymocd in research, with the BibTeX entry of the HP-MOCD paper and of every paper whose method pymocd re-implements.
---

# Citing pymocd

The paper to cite for pymocd is the HP-MOCD paper, Santos et al. (2025),
published in *Social Network Analysis and Mining*
([doi:10.1007/s13278-025-01519-7](https://doi.org/10.1007/s13278-025-01519-7)):
the repository README asks anyone who uses any of the library's algorithms in
research to cite it. If your work uses one of the re-implemented baselines,
also cite the paper that defines that method; the entries are listed
[below](#citing-a-re-implemented-baseline). This page also explains what to do
for RIMPSO, whose paper is not yet published, and for the software itself.

## The pymocd paper

!!! cite "How to cite"
    Santos, G. O., Vieira, L. S., Rossetti, G., Ferreira, C. H. G., &
    Moreira, G. J. P. (2025). A high-performance evolutionary multiobjective
    community detection algorithm. *Social Network Analysis and Mining*,
    15(1), 110.
    [https://doi.org/10.1007/s13278-025-01519-7](https://doi.org/10.1007/s13278-025-01519-7)

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

This is the entry from the repository README. Every field agrees with the
Crossref record of the DOI: the authors, title, journal, volume 15, issue 1,
the electronic ISSN and the online publication date of 18 November 2025. The
value in `pages` is the article number: Crossref records 110 as the article
number and gives no page range.

Some of the same authors have also published two short Portuguese-language
summaries of this paper, listed on the [Publications](publications.md#hp-mocd)
page. Both refer readers to the journal paper for the full results; cite the
journal paper.

### Report the version you used

Defaults, selection rules and detector names have changed between releases
(see the [Changelog](changelog.md)), so state the pymocd version next to the
citation. The installed version is reported by the package metadata. This site
was built with:

```python exec="on" source="above" session="citing" result="text"
from importlib.metadata import version

print(version("pymocd"))
```

Some detectors, HP-MOCD among them, can return a different partition on each
call; see
[Performance and reproducibility](guide/performance.md#determinism) for
which detectors are deterministic.

## Citing RIMPSO

RIMPSO, the library's recommended default detector, has not been published
yet: the repository cites it as *Santos, in preparation (2026)*, and there is
no DOI or BibTeX record for it. Until its paper appears, cite the pymocd paper
above, as the README asks for any of the library's algorithms, refer to the
method as RIMPSO, and give the pymocd version you used. The method itself is
described on the [RIMPSO](algorithms/rimpso.md) page.
<!-- TODO(docs): add the RIMPSO citation, DOI and doi.org BibTeX here once the paper is published; README.md and src/core/algorithms/README.md only say "Santos, in prep. (2026)". -->

## Citing a re-implemented baseline

If you use a re-implemented baseline, also cite its original paper. These
detectors are ports written from the papers, and each
[algorithm page](algorithms/index.md#re-implemented-baselines) lists where the
port deliberately diverges from its paper. When you report results obtained
with them, say that you used the pymocd implementation.

| Method | Entry points | Original paper | DOI |
|---|---|---|---|
| [CDRME](algorithms/cdrme.md) | [`cdrme`][pymocd.cdrme] | Dabaghi-Zarandi, Afkhami & Ashoori (2025), *Journal of Network and Computer Applications* 234, 104070 | [10.1016/j.jnca.2024.104070](https://doi.org/10.1016/j.jnca.2024.104070) |
| [MMCoMO](algorithms/mmcomo.md) | [`mmcomo`][pymocd.mmcomo], [`mmcomo_fronts`][pymocd.mmcomo_fronts] | Zhang, Yang, Yang & Zhang (2023), *IEEE Computational Intelligence Magazine* 18(3), 69–86 | [10.1109/MCI.2023.3277773](https://doi.org/10.1109/MCI.2023.3277773) |
| [NSGA-III-CCM](algorithms/ccm.md) | [`ccm`][pymocd.ccm], [`ccm_fronts`][pymocd.ccm_fronts] | Shaik, Ravi & Deb (2021), *SN Computer Science* 2(1), article 1 | [10.1007/s42979-020-00382-x](https://doi.org/10.1007/s42979-020-00382-x) |
| [NSGA-III-KRM](algorithms/krm.md) | [`krm`][pymocd.krm], [`krm_fronts`][pymocd.krm_fronts] | Shaik, Ravi & Deb (2021), same paper as NSGA-III-CCM | [10.1007/s42979-020-00382-x](https://doi.org/10.1007/s42979-020-00382-x) |
| [GDPSO](algorithms/gdpso.md) | [`gdpso`][pymocd.gdpso] | Cai, Gong, Ma, Ruan, Yuan & Jiao (2015), *Information Sciences* 316, 503–516 | [10.1016/j.ins.2014.09.041](https://doi.org/10.1016/j.ins.2014.09.041) |
| [Shi-MOCD](algorithms/shi-mocd.md) | [`mocd_q`][pymocd.mocd_q], [`mocd_d`][pymocd.mocd_d] | Shi, Yan, Cai & Wu (2012), *Applied Soft Computing* 12(2), 850–859 | [10.1016/j.asoc.2011.10.005](https://doi.org/10.1016/j.asoc.2011.10.005) |
| [MOGA-Net](algorithms/moga-net.md) | [`moga_net`][pymocd.moga_net], [`moga_net_fronts`][pymocd.moga_net_fronts] | Pizzuti (2012), *IEEE Transactions on Evolutionary Computation* 16(3), 418–430 | [10.1109/TEVC.2011.2161090](https://doi.org/10.1109/TEVC.2011.2161090) |

Volumes, issues, pages and article numbers follow the Crossref record of each
DOI. The NSGA-III-CCM and NSGA-III-KRM paper was published online in November
2020 and appeared in the 2021 volume, which is why its BibTeX entry below says
2020. The MOGA-Net module README also cites Pizzuti's earlier conference paper
at ICTAI 2009
([doi:10.1109/ICTAI.2009.58](https://doi.org/10.1109/ICTAI.2009.58)), which
defines the operators, and treats the 2012 journal paper as the authority for
the parameters.

The entries below are the publisher records returned by
[doi.org](https://www.doi.org/) content negotiation for each DOI, laid out one
field per line, with two changes. The GDPSO and MOGA-Net (2012) entries give
the month as `sep` and `jun`, because doi.org writes `Sept` and `June`, which
are not BibTeX month macros. The NSGA-III-CCM and NSGA-III-KRM entry adds a
`pages` field, which the doi.org record leaves out; its value is the article
number 1, from Crossref.

??? cite "BibTeX: CDRME (Dabaghi-Zarandi et al., 2025)"
    ```bibtex
    @article{Dabaghi_Zarandi_2025,
      title     = {Community Detection method based on Random walk and Multi objective Evolutionary algorithm in complex networks},
      volume    = {234},
      ISSN      = {1084-8045},
      url       = {http://dx.doi.org/10.1016/j.jnca.2024.104070},
      DOI       = {10.1016/j.jnca.2024.104070},
      journal   = {Journal of Network and Computer Applications},
      publisher = {Elsevier BV},
      author    = {Dabaghi-Zarandi, Fahimeh and Afkhami, Mohammad Mehdi and Ashoori, Mohammad Hossein},
      year      = {2025},
      month     = Feb,
      pages     = {104070}
    }
    ```

??? cite "BibTeX: MMCoMO (Zhang et al., 2023)"
    ```bibtex
    @article{Zhang_2023,
      title     = {A Macro-Micro Population-Based Co-Evolutionary Multi-Objective Algorithm for Community Detection in Complex Networks [Research Frontier]},
      volume    = {18},
      ISSN      = {1556-6048},
      url       = {http://dx.doi.org/10.1109/MCI.2023.3277773},
      DOI       = {10.1109/mci.2023.3277773},
      number    = {3},
      journal   = {IEEE Computational Intelligence Magazine},
      publisher = {Institute of Electrical and Electronics Engineers (IEEE)},
      author    = {Zhang, Lei and Yang, Haipeng and Yang, Shangshang and Zhang, Xingyi},
      year      = {2023},
      month     = Aug,
      pages     = {69–86}
    }
    ```

??? cite "BibTeX: NSGA-III-CCM and NSGA-III-KRM (Shaik et al., 2021)"
    ```bibtex
    @article{Shaik_2020,
      title     = {Evolutionary Multi-Objective Optimization Algorithm for Community Detection in Complex Social Networks},
      volume    = {2},
      ISSN      = {2661-8907},
      url       = {http://dx.doi.org/10.1007/s42979-020-00382-x},
      DOI       = {10.1007/s42979-020-00382-x},
      number    = {1},
      journal   = {SN Computer Science},
      publisher = {Springer Science and Business Media LLC},
      author    = {Shaik, Tanveerulhuq and Ravi, Vadlamani and Deb, Kalyanmoy},
      year      = {2020},
      month     = Nov,
      pages     = {1}
    }
    ```

??? cite "BibTeX: GDPSO (Cai et al., 2015)"
    ```bibtex
    @article{Cai_2015,
      title     = {Greedy discrete particle swarm optimization for large-scale social network clustering},
      volume    = {316},
      ISSN      = {0020-0255},
      url       = {http://dx.doi.org/10.1016/j.ins.2014.09.041},
      DOI       = {10.1016/j.ins.2014.09.041},
      journal   = {Information Sciences},
      publisher = {Elsevier BV},
      author    = {Cai, Qing and Gong, Maoguo and Ma, Lijia and Ruan, Shasha and Yuan, Fuyan and Jiao, Licheng},
      year      = {2015},
      month     = sep,
      pages     = {503–516}
    }
    ```

??? cite "BibTeX: Shi-MOCD (Shi et al., 2012)"
    ```bibtex
    @article{Shi_2012,
      title     = {Multi-objective community detection in complex networks},
      volume    = {12},
      ISSN      = {1568-4946},
      url       = {http://dx.doi.org/10.1016/j.asoc.2011.10.005},
      DOI       = {10.1016/j.asoc.2011.10.005},
      number    = {2},
      journal   = {Applied Soft Computing},
      publisher = {Elsevier BV},
      author    = {Shi, Chuan and Yan, Zhenyu and Cai, Yanan and Wu, Bin},
      year      = {2012},
      month     = Feb,
      pages     = {850–859}
    }
    ```

??? cite "BibTeX: MOGA-Net (Pizzuti, 2012)"
    ```bibtex
    @article{Pizzuti_2012,
      title     = {A Multiobjective Genetic Algorithm to Find Communities in Complex Networks},
      volume    = {16},
      ISSN      = {1941-0026},
      url       = {http://dx.doi.org/10.1109/TEVC.2011.2161090},
      DOI       = {10.1109/tevc.2011.2161090},
      number    = {3},
      journal   = {IEEE Transactions on Evolutionary Computation},
      publisher = {Institute of Electrical and Electronics Engineers (IEEE)},
      author    = {Pizzuti, Clara},
      year      = {2012},
      month     = jun,
      pages     = {418–430}
    }
    ```

??? cite "BibTeX: MOGA-Net operators (Pizzuti, 2009)"
    ```bibtex
    @inproceedings{Pizzuti_2009,
      title     = {A Multi-objective Genetic Algorithm for Community Detection in Networks},
      url       = {http://dx.doi.org/10.1109/ICTAI.2009.58},
      DOI       = {10.1109/ictai.2009.58},
      booktitle = {2009 21st IEEE International Conference on Tools with Artificial Intelligence},
      publisher = {IEEE},
      author    = {Pizzuti, Clara},
      year      = {2009},
      month     = Nov,
      pages     = {379–386}
    }
    ```

## Citing the software

The repository root contains a
[`CITATION.cff`](https://github.com/csilab-ufop/pymocd/blob/master/CITATION.cff)
file in version 1.2.0 of the
[Citation File Format](https://citation-file-format.github.io/). It describes
the software and names the pymocd paper [above](#the-pymocd-paper) as its
preferred citation. GitHub reads the file from the default branch and adds a
*Cite this repository* link to the right sidebar of the
[repository page](https://github.com/csilab-ufop/pymocd), which exports the
citation as APA or BibTeX. Because the file sets a preferred citation, that
link gives the pymocd paper, not the software.

No release has a DOI yet, so there is no archived, citable copy of the
software. Cite the pymocd paper. If a venue also asks for a software
reference, `CITATION.cff` records the following:
<!-- TODO(docs): archive releases with a DOI (e.g. through the Zenodo GitHub integration, which needs the maintainer's GitHub and Zenodo accounts), then add that DOI to CITATION.cff and give the DOI and its BibTeX here. -->

Title
:   pymocd

Author
:   Guilherme O. Santos, as in `pyproject.toml`

Version
:   6.1.0, released on 31 October 2026, the release the file was written
    for. If you used another version, give the one printed by
    `version("pymocd")` [above](#report-the-version-you-used).

Repository
:   [https://github.com/csilab-ufop/pymocd](https://github.com/csilab-ufop/pymocd)

Documentation
:   [https://csilab.ufop.br/pymocd/](https://csilab.ufop.br/pymocd/)

License
:   GPL-3.0-or-later

The package itself is published on PyPI as
[pymocd](https://pypi.org/project/pymocd/).

## See also

- [Publications](publications.md): the papers behind pymocd and the methods it
  re-implements.
- [Algorithms](algorithms/index.md): one page per detector, each ending with
  the full reference of its paper.
- [Changelog](changelog.md): what changed in each release.
