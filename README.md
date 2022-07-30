# MEGH R package

## MEGH: A parametric class of general hazard models for clustered survival data

This repository contains a real data application of the clustered survival models proposed in the paper (Open access):

 > Rubio, F.J. and Drikvandi, R. (2022). MEGH: A parametric class of general hazard models for clustered survival data. Statistical Methods in Medical Research, 31(8) 1603–1616. https://doi.org/10.1177%2F09622802221102620

[[Preprint](https://drive.google.com/file/d/1YjHkOKYWK_4gZNt8kMAk6YqndI9mJ_aI/view?usp=sharing)] [[Supplementary Material](https://drive.google.com/file/d/1A4V3eRCl23tinv7XLKAk4Fek6d0q0Uwe/view?usp=sharing)] [[Journal](https://doi.org/10.1177%2F09622802221102620)]

The models are fitted using the R package `MEGH`. To install the `MEGH` R package use:

```
library(devtools)
install_github("FJRubio67/MEGH")

library(MEGH)
```

An example of the use of this package using real data can be found at:

[MEGH: A parametric class of general hazard models for clustered survival, Leukemia data set](https://rpubs.com/FJRubio/MEGHLeuk)

See also: [GHSurv](https://github.com/FJRubio67/GHSurv), [HazReg](https://github.com/FJRubio67/HazReg)
