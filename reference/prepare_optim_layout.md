# Prepare optimization layout for a model

Creates a layout to map a flat parameter vector `theta` (real space)
into component parameters for a `sum_model`.

## Usage

``` r
prepare_optim_layout(model)
```

## Arguments

- model:

  A `time_series_model` or `sum_model`.

## Value

A list with `kind`, `theta0`, and (for sums) a `layout`.

## Examples

``` r
mod <- wn(1) + pl(kappa = 0.5, sigma2 = 2)
```
