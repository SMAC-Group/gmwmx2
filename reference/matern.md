# Matérn process

Constructs a `time_series_model` for a Matérn covariance process with
variance `sigma2`, range `lambda`, and smoothness `alpha`. The
autocovariance is \\ \gamma(h) = \mathrm{cov}(X_t, X\_{t+h}) = \frac{2
\sigma^2}{\Gamma(\alpha-1 / 2) 2^{\alpha-1 / 2}}\|\lambda h\|^{\alpha-1
/ 2} \mathcal{K}\_{\|\alpha-1 / 2\|}(\| \lambda h\|)\\ where
\\\mathcal{K}\_\omega(x)\\ is the modified Bessel function of the second
kind of order \\\omega\\.

## Usage

``` r
matern(sigma2 = NULL, lambda = NULL, alpha = NULL)
```

## Arguments

- sigma2:

  Marginal variance (\> 0).

- lambda:

  Range/scale parameter (\> 0).

- alpha:

  Smoothness parameter (\> 1/2).

## Value

A `time_series_model` object.

## Examples

``` r
mod <- matern(sigma2 = 1, lambda = 0.2, alpha = 1.0)
mod
#> Stochastic process
#>   Model      : Matérn 
#>   Parameters : sigma2 =     1, lambda =   0.2, alpha =     1 
```
