
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SGDinference

<!-- badges: start -->

[![R-CMD-check](https://github.com/SGDinference-Lab/SGDinference/workflows/R-CMD-check/badge.svg)](https://github.com/SGDinference-Lab/SGDinference/actions)
[![codecov](https://codecov.io/gh/SGDinference-Lab/SGDinference/branch/master/graph/badge.svg?token=YTBY15IXEP)](https://codecov.io/gh/SGDinference-Lab/SGDinference)
<!-- badges: end -->

**SGDinference** is an R package that provides estimation and inference
methods for large-scale mean and quantile regression models via
stochastic (sub-)gradient descent (S-subGD) algorithms. The inference
procedure handles cross-sectional data sequentially:

1)  updating the parameter estimate with each incoming “new
    observation”,
2)  aggregating it as a Polyak-Ruppert average, and
3)  computing an asympotically pivotal statistic for inference through
    random scaling.

The methodology used in the SGDinference package is described in detail
in the following papers:

- Lee, S., Liao, Y., Seo, M.H. and Shin, Y., 2022. Fast and robust
  online inference with stochastic gradient descent via random scaling.
  In *Proceedings of the AAAI Conference on Artificial Intelligence*
  (Vol. 36, No. 7, pp. 7381-7389).
  <https://doi.org/10.1609/aaai.v36i7.20701>.

- Lee, S., Liao, Y., Seo, M.H. and Shin, Y., 2023. Fast Inference for
  Quantile Regression with Tens of Millions of Observations.
  arXiv:2209.14502 \[econ.EM\]
  <https://doi.org/10.48550/arXiv.2209.14502>.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools") # if you have not installed "devtools" package
devtools::install_github("SGDinference-Lab/SGDinference")
```
