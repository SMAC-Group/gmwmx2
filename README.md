# `gmwmx2` Overview <img src="man/figures/logo.png" align="right" style="width: 20%; height: 20%"/>

The `gmwmx2` `R` package implements the Generalized Method of Wavelet Moments with Exogenous Inputs estimator (GMWMX) presented in [Voirol, L., Xu, H., Zhang, Y., Insolia, L., Molinari, R. and Guerrier, S. (2024)](https://arxiv.org/abs/2409.05160) and provides functions to estimate times series models that can be expressed as linear models with correlated residuals in presence of missing data.
Moreover, the `gmwmx2` package provides tools to load an plot Global Navigation Satellite System (GNSS) data from the Nevada Geodetic Laboratory.

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


## License

This source code is released under is the GNU AFFERO GENERAL PUBLIC LICENSE (AGPL) v3.0. 

## References
Voirol, L., Xu, H., Zhang, Y., Insolia, L., Molinari, R., and Guerrier, S. (2024). Inference for Large Scale Regression Models with Dependent Errors. [doi:10.48550/arXiv.2409.05160](https://doi.org/10.48550/arXiv.2409.05160).

