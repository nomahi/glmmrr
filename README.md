
# glmmrr package


## Logistic mixed-effects model analysis with pseudo-observations for estimating risk ratios in clustered binary data analysis

Logistic mixed-effects model has been widely used as a multilevel statistical model for analyzing clustered binary outcome data (e.g., longitudinal studies, cluster-randomised trials, and multi-center clinical trials). However, the resultant odds ratio estimator can only be interpreted as an approximation of the risk ratio estimator for low-frequency events; it cannot be directly interpreted as an effect measure. To overcome this issue, the modified Poisson regression analysis and its extention to GEE methodology has been widely applied in recent clinical and epidemiological studies, but these estimating equation-based semiparametric methods cannot be straightforwardly extended to mixed-effects models. Noma (2025) proposed a new method to provide consistent risk ratio estimator on a multilevel statistical model framework using logistic mixed-effects model incorporating pseudo-observations. The advantage of the new method is it is implementable using a standard statistical software of GLMM only through modifying dataset. This package involves computational functions for implementing the risk ratio estimation method through multilevel modelling framework.



## Installation

Please download "glmmrr_1.1-1.tar.gz" and install it by R menu: "packages" -> "Install package(s) from local files...".

Download: [please click this link](https://github.com/nomahi/glmmrr/raw/main/glmmrr_1.1-1.tar.gz)

Manual: [please click this link](https://github.com/nomahi/glmmrr/raw/main/glmmrr_1.1-1.pdf)

Example code: [please click this link](https://github.com/nomahi/glmmrr/raw/main/examplecode.r)
