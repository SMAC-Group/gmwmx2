# Theoretical wavelet variance from model parameters

Computes the theoretical WV implied by a model and parameter vector.

## Usage

``` r
get_theoretical_wv(theta, model, n, wv_obj = NULL, tau = NULL, prep = NULL)
```

## Arguments

- theta:

  Real-valued parameter vector (unconstrained).

- model:

  A `time_series_model` or `sum_model`.

- n:

  Length of autocovariance to compute.

- wv_obj:

  Optional [`wv::wvar`](https://rdrr.io/pkg/wv/man/wvar.html) object
  (uses its scales).

- tau:

  Optional vector of scales; overrides `wv_obj` if provided.

- prep:

  Optional output from `prepare_optim_layout` for sum models.

## Value

Numeric vector of theoretical wavelet variances.
