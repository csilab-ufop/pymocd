---
description: The papers that describe pymocd's own community detection methods, the published methods it re-implements, and how to report work that uses pymocd.
---

# Publications

This page lists the papers behind pymocd. It gives the published paper on
HP-MOCD, with two short summaries of it, and the publication status of
RIMPSO, the library's two own methods, then the papers whose methods pymocd
re-implements as baselines. A last section is reserved for published work
that uses pymocd.
For BibTeX entries and guidance on what to cite, see
[Citing pymocd](citing.md).

## pymocd's methods

### HP-MOCD

Santos, G. O., Vieira, L. S., Rossetti, G., Ferreira, C. H. G., & Moreira,
G. J. P. (2025). A high-performance evolutionary multiobjective community
detection algorithm. *Social Network Analysis and Mining*, 15(1), 110.
[doi:10.1007/s13278-025-01519-7](https://doi.org/10.1007/s13278-025-01519-7)

This is the paper the repository README asks users of any of the library's
algorithms to cite. It introduces HP-MOCD, an NSGA-II search over the
intra/inter decomposition of modularity, which pymocd
provides as [`hpmocd`][pymocd.hpmocd], [`hpmocd_fronts`][pymocd.hpmocd_fronts]
and the configurable [`HpMocd`][pymocd.HpMocd] class. The Crossref record of
the DOI lists the article under a CC BY 4.0 licence. A preprint of the same
paper is available as
[arXiv:2506.01752](https://arxiv.org/abs/2506.01752); the arXiv record links
to the published version.

Some of the same authors have since published two short Portuguese-language
summaries of the paper. Both present HP-MOCD and refer readers to the journal
paper for the full results, so cite the journal paper for the method.

- Santos, G. O., Ferreira, C. H. G., & Moreira, G. J. P. (2026). HP-MOCD: Um
  Algoritmo Evolucionário Multiobjetivo de Alto Desempenho para Detecção de
  Comunidades em Redes Complexas. *Revista Eletrônica de Iniciação Científica
  em Computação*, 24, 418–424.
  [doi:10.5753/reic.2026.8488](https://doi.org/10.5753/reic.2026.8488)
- Santos, G. O., Vieira, L. S., Ferreira, C. H. G., & Moreira, G. J. P.
  (2026). HP-MOCD: Viabilizando a Detecção de Comunidades Multiobjetivo em
  Redes de Larga Escala. In *Anais do XV Brazilian Workshop on Social Network
  Analysis and Mining (BraSNAM 2026)*, 262–268. Sociedade Brasileira de
  Computação.
  [doi:10.5753/brasnam.2026.22901](https://doi.org/10.5753/brasnam.2026.22901)

The first names pymocd as the library that provides HP-MOCD and gives
`github.com/oliveira-sh/pymocd`, now a fork of this repository, as the
location of the source code. Bibliographic details follow the Crossref record
of each DOI.

[HP-MOCD algorithm page](algorithms/hpmocd.md){ .md-button }
[BibTeX](citing.md#the-pymocd-paper){ .md-button }

### RIMPSO

RIMPSO, which its module README expands as *RIMPSO-CD: Resolution-Indexed
Memetic Particle Swarm Optimisation for Community Detection*, is the library's
recommended default detector, provided as [`rimpso`][pymocd.rimpso],
[`rimpso_fronts`][pymocd.rimpso_fronts] and
[`rimpso_select`][pymocd.rimpso_select]. Its paper is in preparation: the
repository cites it as *Santos, in preparation (2026)*, and no DOI exists
yet. Until it is published, the [RIMPSO](algorithms/rimpso.md) page and the
module README in `src/core/algorithms/rimpso/` are its full description.
<!-- TODO(docs): add the RIMPSO paper's full citation and DOI once it is published; README.md and src/core/algorithms/README.md only say "Santos, in prep. (2026)". -->

## Papers re-implemented in pymocd

pymocd re-implements the methods of the following papers as baselines. Each
algorithm page gives the method's objectives, its parameters and the
divergences of the port from its paper that the module README records.
Bibliographic details follow the Crossref record of each DOI.

| Method | Paper | Algorithm page |
|---|---|---|
| CDRME | Dabaghi-Zarandi, F., Afkhami, M. M., & Ashoori, M. H. (2025). Community Detection method based on Random walk and Multi objective Evolutionary algorithm in complex networks. *Journal of Network and Computer Applications*, 234, 104070. [doi:10.1016/j.jnca.2024.104070](https://doi.org/10.1016/j.jnca.2024.104070) | [CDRME](algorithms/cdrme.md) |
| MMCoMO | Zhang, L., Yang, H., Yang, S., & Zhang, X. (2023). A Macro-Micro Population-Based Co-Evolutionary Multi-Objective Algorithm for Community Detection in Complex Networks [Research Frontier]. *IEEE Computational Intelligence Magazine*, 18(3), 69–86. [doi:10.1109/MCI.2023.3277773](https://doi.org/10.1109/MCI.2023.3277773) | [MMCoMO](algorithms/mmcomo.md) |
| NSGA-III-CCM, NSGA-III-KRM | Shaik, T., Ravi, V., & Deb, K. (2021). Evolutionary Multi-Objective Optimization Algorithm for Community Detection in Complex Social Networks. *SN Computer Science*, 2(1), 1. [doi:10.1007/s42979-020-00382-x](https://doi.org/10.1007/s42979-020-00382-x) | [NSGA-III-CCM](algorithms/ccm.md), [NSGA-III-KRM](algorithms/krm.md) |
| GDPSO | Cai, Q., Gong, M., Ma, L., Ruan, S., Yuan, F., & Jiao, L. (2015). Greedy discrete particle swarm optimization for large-scale social network clustering. *Information Sciences*, 316, 503–516. [doi:10.1016/j.ins.2014.09.041](https://doi.org/10.1016/j.ins.2014.09.041) | [GDPSO](algorithms/gdpso.md) |
| Shi-MOCD (MOCD-Q, MOCD-D) | Shi, C., Yan, Z., Cai, Y., & Wu, B. (2012). Multi-objective community detection in complex networks. *Applied Soft Computing*, 12(2), 850–859. [doi:10.1016/j.asoc.2011.10.005](https://doi.org/10.1016/j.asoc.2011.10.005) | [Shi-MOCD](algorithms/shi-mocd.md) |
| MOGA-Net | Pizzuti, C. (2012). A Multiobjective Genetic Algorithm to Find Communities in Complex Networks. *IEEE Transactions on Evolutionary Computation*, 16(3), 418–430. [doi:10.1109/TEVC.2011.2161090](https://doi.org/10.1109/TEVC.2011.2161090) | [MOGA-Net](algorithms/moga-net.md) |

The Shaik et al. paper was published online in November 2020 and appeared in
volume 2 (2021) as article 1. For MOGA-Net, the module README also draws on
Pizzuti's earlier conference paper, "A Multi-objective Genetic Algorithm for
Community Detection in Networks", *21st IEEE International Conference on Tools
with Artificial Intelligence (ICTAI)*, 2009, pp. 379–386
([doi:10.1109/ICTAI.2009.58](https://doi.org/10.1109/ICTAI.2009.58)), which
defines the operators.

HP-MOCD's objective pair follows the formulation of Shi, Zhong, Yan, Cai &
Wu (2010), "A multi-objective approach for community detection in complex
network", *IEEE Congress on Evolutionary Computation*, pp. 1–8
([doi:10.1109/CEC.2010.5585987](https://doi.org/10.1109/CEC.2010.5585987)),
which is the Shi et al. paper the HP-MOCD paper cites. The module README takes
the equations from Shi et al.'s later journal paper, Shi et al. (2012),
Eqs. 3.5 and 3.6, which is the paper Shi-MOCD ports. The
[HP-MOCD](algorithms/hpmocd.md) page lists the works it builds on.

## Publications using pymocd

No published work by other authors that uses pymocd has been confirmed yet.
Of the works found that cite the HP-MOCD paper, those whose full text could
be read cite it without using pymocd.
<!-- TODO(docs): list peer-reviewed work by other authors that uses pymocd (full citation and DOI). Searched 2026-10-31: works citing doi:10.1007/s13278-025-01519-7 in OpenAlex (cites:W4414898022, 1 work), Semantic Scholar (7, one of them the authors' REIC summary) and Crossref (cited-by count 5, list not public); full-text search for "pymocd" and "HP-MOCD" in OpenAlex (besides the HP-MOCD paper itself, only the authors' REIC and BraSNAM summaries; "pymocd" also matched an unrelated chemistry review), Semantic Scholar keyword search (0; its full-text snippet search was rate-limited) and Europe PMC (0). Full text read, cited only: Wahid Ahmed & Hama Ali Faraj 2026 (doi:10.5121/ijaia.2026.17405), Amedapu & Velusamy 2026 (doi:10.3390/app16073329), de Oliveira et al. 2026 (doi:10.1145/3795766.3799768, arXiv:2604.18586 version read). Paywalled, not checked: Rahimi et al. 2026 (doi:10.1007/s00521-026-12599-7), Sun et al. 2026 (doi:10.1007/s41060-026-01128-3), Li et al. 2026 (doi:10.1016/j.ins.2026.123170). Not yet published: OHP-MOCD (github.com/aanasakhtar/ohp-mocd), a fork of pymocd whose README cites a Springer LNCS paper with no DOI; list it once it has one. -->

If you have published work that uses pymocd, please let us know by
[opening an issue](https://github.com/csilab-ufop/pymocd/issues) on GitHub
with the full citation and its DOI, so that it can be listed here.

## See also

- [Citing pymocd](citing.md): what to cite, with BibTeX entries for every
  paper on this page.
- [Algorithms](algorithms/index.md): the library's methods and the
  re-implemented baselines side by side.
