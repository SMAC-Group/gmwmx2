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
m <- wn(sigma2 = 1) + ar1(phi = 0.8, sigma2 = 0.5)
x <- generate(m, n = 500, seed = 123)
plot(x)

plot(wv::wvar(x$series))

fit <- gmwm2(x, m)
fit$theta_domain
#> $`White Noise_1`
#>   sigma2 
#> 1.082831 
#> 
#> $`AR(1)_2`
#>       phi    sigma2 
#> 0.7764063 0.5182533 
#> 
```
