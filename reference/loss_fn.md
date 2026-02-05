# Loss function for GMWM2 optimization (internal)

Computes the weighted squared error between empirical wavelet variance
and theoretical wavelet variance implied by a model and parameter
vector.

## Usage

``` r
loss_fn(theta, model, n, prep, wv_obj, omega = NULL)
```

## Arguments

- theta:

  Real-valued parameter vector.

- model:

  A `time_series_model` or `sum_model`.

- n:

  Length of autocovariance to compute.

- prep:

  Output from `prepare_optim_layout`.

- wv_obj:

  A [`wv::wvar`](https://rdrr.io/pkg/wv/man/wvar.html) object.

- omega:

  Optional weighting matrix. If `NULL`, uses inverse CI width.

## Value

Scalar objective value.
