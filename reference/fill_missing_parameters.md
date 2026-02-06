# Fill missing model parameters using initial-parameter functions

Ensures any NULL / missing parameters are populated from the model's
`get_initial_parameters_function(signal)` while preserving any
user-provided parameters (assumed to be in the domain).

## Usage

``` r
fill_missing_parameters(model, signal)
```

## Arguments

- model:

  A `time_series_model` or `sum_model`.

- signal:

  Numeric vector used to derive initial parameters.

## Value

Model with all parameters populated.
