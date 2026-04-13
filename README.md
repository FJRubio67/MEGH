# MEGH: General Hazard Models for Clustered Survival Data (R package)

[![R package](https://img.shields.io/badge/language-R-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Open Access](https://img.shields.io/badge/Open%20Access-✓-green.svg)](https://doi.org/10.1177/09622802221102620)

<div align="center">
<img src="https://github.com/FJRubio67/MEGH/assets/33929387/9aa5d701-bbe9-4aea-a732-4f65905ba5e0" width="300" height="300" alt="MEGH model illustration"/>
</div>

## Overview

`MEGH` is an R package implementing a parametric class of **Mixed Effects General
Hazard (MEGH)** models for clustered survival data. It extends the General
Hazard (GH) framework, which nests the Proportional Hazards (PH),
Accelerated Failure Time (AFT), and Accelerated Hazards (AH) models as special
cases, to settings where individuals are grouped into clusters (e.g. patients
within hospitals, or subjects within treatment centres), by incorporating
**cluster-specific random effects** into the hazard structure.

The package is based on the methodology proposed in:

> Rubio, F.J. and Drikvandi, R. (2022). MEGH: A parametric class of general
> hazard models for clustered survival data. *Statistical Methods in Medical
> Research* **31**(8): 1603–1616.
> [DOI: 10.1177/09622802221102620](https://doi.org/10.1177/09622802221102620)
> *(Open access)* |
> [Preprint](https://drive.google.com/file/d/1YjHkOKYWK_4gZNt8kMAk6YqndI9mJ_aI/view?usp=sharing) |
> [Supplementary Material](https://drive.google.com/file/d/1A4V3eRCl23tinv7XLKAk4Fek6d0q0Uwe/view?usp=sharing)

## Installation

```r
# install.packages("devtools")
devtools::install_github("FJRubio67/MEGH")
library(MEGH)
```

## Main function

| Function | Description |
|---|---|
| `MEGHMLE` | Maximum likelihood estimation for MEGH models |

For full documentation: `?MEGHMLE`

## Supported hazard structures and random effect distributions

The MEGH framework supports the following cluster-level hazard models:

| Model | Abbreviation |
|---|---|
| General Hazard | GH |
| Proportional Hazards | PH |
| Accelerated Failure Time | AFT |
| Accelerated Hazards | AH |

Random effects can be specified at the cluster level to capture within-cluster
dependence in survival times.

## Tutorial

- [MEGH: Leukemia data set](https://rpubs.com/FJRubio/MEGHLeuk) — illustrative
  real-data application of MEGH models to leukaemia survival data

## Citation

If you use this package, please cite:

```bibtex
@article{rubio:2022,
  author  = {Rubio, F.J. and Drikvandi, R.},
  title   = {{MEGH}: A parametric class of general hazard models for clustered
             survival data},
  journal = {Statistical Methods in Medical Research},
  volume  = {31},
  number  = {8},
  pages   = {1603--1616},
  year    = {2022},
  doi     = {10.1177/09622802221102620}
}
```

## Related resources

- [HazReg](https://github.com/FJRubio67/HazReg) — GH models for overall
  survival (no random effects; required dependency)
- [ShortCourseParamSurvival](https://github.com/FJRubio67/ShortCourseParamSurvival) —
  short course covering the GH framework on which MEGH is built

## License

This package is licensed under the [MIT License](LICENSE).
