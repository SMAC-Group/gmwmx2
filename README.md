# `gmwmx2` Overview <img src="man/figures/logo.png" align="right" style="width: 15%; height: 15%"/>

<!-- badges: start 

-->
![](https://img.shields.io/github/last-commit/SMAC-Group/gmwmx2) 
[<img src="https://s-a.github.io/license/img/agpl-3.0.svg" />](https://s-a.github.io/license/?license=agpl-3.0&fullname=Stephan%20Ahlf&year=2015&profile=https://github.com/s-a&projectUrl=https://github.com/s-a/license&projectName=License%20Demo "")
![R-CMD-check](https://github.com/SMAC-Group/gmwmx2/actions/workflows/R-CMD-check.yaml/badge.svg)
[![minimal R version](https://img.shields.io/badge/R%3E%3D-4.0.0-6666ff.svg)](https://cran.r-project.org/)
<!-- badges: end -->

The `gmwmx2` `R` package implements the Generalized Method of Wavelet Moments with Exogenous Inputs estimator (GMWMX) presented in [Voirol, L., Xu, H., Zhang, Y., Insolia, L., Molinari, R. and Guerrier, S. (2024)](https://arxiv.org/abs/2409.05160).
The GMWMX estimator is a computationally efficient estimator to estimate large scale regression problems with complex dependence structure in presence of missing data.
The `gmwmx2` `R` package  allows to estimate (i) functional/structural parameters, (ii) stochastic parameters describing the dependence structure and (iii) nuisance parameters of the missingness process of large regression models with dependent observations and missing data.
To illustrate the capability of the GMWMX estimator, the `gmwmx2` `R` package provides functions to download an plot Global Navigation Satellite System (GNSS) position time series from the [Nevada Geodetic Laboratory](https://geodesy.unr.edu/) and allow to estimate linear model with a specific dependence structure modeled by composite stochastic processes, allowing to estimate tectonic velocities and crustal uplift from GNSS position time series.

Find vignettes with detailed examples as well as the user's manual at the [package website](https://smac-group.github.io/gmwmx2/index.html).

Below are instructions on how to install and make use of the `gmwmx2` package.

## Installation Instructions

The `gmwmx2` package is currently only available on GitHub. You can install the `gmwmx2` package with:

``` r
# Install dependencies
install.packages(c("devtools"))

# Install/Update the package from GitHub
devtools::install_github("SMAC-Group/gmwmx2")

# Install the package with Vignettes/User Guides 
devtools::install_github("SMAC-Group/gmwmx2", build_vignettes = TRUE)
```

### External `R` libraries

The `gmwmx2` package relies on a limited number of external libraries, but notably on `Rcpp` and `RcppArmadillo` which require a `C++` compiler for installation, such as for example `gcc`.


## Note on `gmwmx2` vs `gmwmx`

The original [`gmwmx`](https://github.com/SMAC-Group/gmwmx) package was designed to compare estimated parameters obtained from the GMWMX with the ones obtained with the Maximum Likelihood Estimator (MLE) implemented in [Hector](https://teromovigo.com/product/hector/). 
This allowed for the replication of examples and simulations discussed in [Cucci, D. A., Voirol, L., Kermarrec, G., Montillet, J. P., and Guerrier, S. (2022)](https://doi.org/10.1007/s00190-023-01702-8).
However, as we advanced in the methodological and computational development of the GMWMX method, we sought a standalone implementation that did not include [Hector](https://teromovigo.com/hector/).
Additionally, many of the new computational techniques and structural improvements would have been challenging to incorporate into the previous `gmwmx` package. 
Therefore, we will now exclusively support and develop the `gmwmx2` package.

## Upcoming features

The `gmwmx2` package is currently in the early stages of development. While the supported features are stable, we have numerous additional methods and computational enhancements planned for gradual integration. These include:

- Computational optimization to improve speed
- Support for a wider range of stochastic models
- A computationally efficient model selection criterion for stochastic models


## License

This source code is released under is the GNU AFFERO GENERAL PUBLIC LICENSE (AGPL) v3.0. 

## References
Voirol, L., Xu, H., Zhang, Y., Insolia, L., Molinari, R., and Guerrier, S. (2024). Inference for Large Scale Regression Models with Dependent Errors. [doi:10.48550/arXiv.2409.05160](https://doi.org/10.48550/arXiv.2409.05160).

Guerrier, S., Skaloud, J., Stebler, Y. and Victoria-Feser, M.P., 2013. Wavelet-variance-based estimation for composite stochastic processes. Journal of the American Statistical Association, 108(503), pp.1021-1030. [doi:10.1080/01621459.2013.799920](https://doi.org/10.1080/01621459.2013.799920)
