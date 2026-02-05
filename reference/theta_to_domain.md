# Convert optimization parameters to domain parameters

Maps a real-valued optimization vector `theta` into the constrained
parameter space used by each model component.

## Usage

``` r
theta_to_domain(model, theta, prep = NULL)
```

## Arguments

- model:

  A `time_series_model` or `sum_model`.

- theta:

  Numeric vector in real space (typically from `optim`).

- prep:

  Optional output from `prepare_optim_layout` for sum models.

## Value

Named numeric vector (single model) or a named list (sum model).
