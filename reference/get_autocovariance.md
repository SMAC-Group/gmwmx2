# Compute autocovariance for a model (internal)

If `theta` is `NULL`, uses domain parameters stored in the model. If
`theta` is provided, it is treated as unconstrained parameters in real
space and mapped to the domain via the model's transformation.

## Usage

``` r
get_autocovariance(object, n, theta = NULL, prep = NULL, ...)
```

## Arguments

- object:

  A `time_series_model` or `sum_model`.

- n:

  Length of autocovariance vector.

- theta:

  Optional real-valued parameter vector for optimization.

- prep:

  Optional output from `prepare_optim_layout` (required for sum models
  with `theta`).

- ...:

  Passed to methods.

## Value

Numeric vector of autocovariances of length `n`.
