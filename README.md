# VBRobM
Variational Bayesian robust mediation analysis. Paper link:

Yu-Quan Wang, Da-Peng Shi, Christoph Scherber, Ben A. Woodcock, Yue-Qing Hu, Nian-Feng Wan. Understanding biodiversity effects on trophic interactions with a robust approach to path analysis. Cell Reports Sustainability, 2025. [https://www.cell.com/cell-reports-sustainability/fulltext/S2949-7906(25)00058-8](https://www.cell.com/cell-reports-sustainability/fulltext/S2949-7906(25)00058-8)

# Installation
To install the package, run the following command in R
```
devtools::install_github("YuquanW/VBRobM")
```

# Examples
The main function is `vbrobm`. You can choose one of `rREML`, `bcrREML` and `vbrob` to fit a linear mixed model. For example,
```
library(VBRobM)
fit <- vbrobm(vbrob(smd_herbivore ~ 1+(1|Code),
                    data = trophic,
                    vi = trophic$vi_herbivore,
                    maxit = 1000),
              vbrob(smd_plant ~ 1+smd_herbivore+(1|Code),
                    data = trophic,
                    vi = trophic$vi_plant,
                    maxit = 1000))
plot(fit)
```
