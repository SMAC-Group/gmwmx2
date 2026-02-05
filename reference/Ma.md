# Matérn autocovariance

Computes the Matérn correlation at lags `x` for smoothness `alpha`. This
is used internally to build the full autocovariance vector.

## Usage

``` r
Ma(x, alpha)
```

## Arguments

- x:

  Numeric vector of lags (non-negative).

- alpha:

  Smoothness parameter (\> 1/2).

## Value

Numeric vector of correlations for each lag in `x`.
