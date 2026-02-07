# Estimate stochastic model parameters by matching wavelet variance

Fits a `time_series_model` or `sum_model` to data by minimizing the
weighted distance between empirical and theoretical wavelet variance.
Optimization is performed in real (unconstrained) space and transformed
to the model's parameter domain internally.

## Usage

``` r
gmwm2(x, model, omega = NULL, method = "L-BFGS-B", control = list(), ...)
```

## Arguments

- x:

  Numeric vector, or a `generated_time_series` /
  `generated_composite_model_time_series` object (its `series` is used).

- model:

  A `time_series_model` or `sum_model`.

- omega:

  Optional weighting matrix. If `NULL`, a default based on the empirical
  WV confidence intervals is used.

- method:

  Optimization method passed to
  [`stats::optim`](https://rdrr.io/r/stats/optim.html).

- control:

  Optional list of control parameters for
  [`stats::optim`](https://rdrr.io/r/stats/optim.html).

- ...:

  Additional arguments passed to
  [`stats::optim`](https://rdrr.io/r/stats/optim.html).

## Value

An object of class `gmwm2_fit` with elements: `theta_hat` (real space),
`theta_domain` (constrained space), `model`, `empirical_wvar`,
`theoretical_wvar`, `optim`, and `n`.

## Details

The default weighting matrix is diagonal with entries proportional to
the inverse squared width of the empirical WV confidence intervals.
Provide `omega` to use a custom weighting (e.g., from a theoretical
covariance).

## Examples

``` r
model <- wn(sigma2 = 1) + ar1(phi = 0.8, sigma2 = 0.5)
x <- generate(model, n = 1000, seed = 123)
plot(x)

fit <- gmwm2(x, model = wn()+ar1())
fit
#> GMWM fit
#> 
#> Stochastic model
#>   Sum of 2 processes
#> 
#>   [1] White Noise
#>       Parameters : sigma2
#> 
#>   [2] AR(1)
#>       Parameters : phi, sigma2
#> 
#> Initial parameters
#>   1) White Noise: sigma2 = 2.424
#>   2) AR(1): phi = -0.1953, sigma2 = 2.424
#> 
#> Estimated parameters
#>   1) White Noise: sigma2 = 0.9147
#>   2) AR(1): phi = 0.7744, sigma2 = 0.5982
#> 
#> Optimization
#>   Convergence : converged (code 0)
#>   Iterations  : 38
#>   Loss        : 0.03692
```
