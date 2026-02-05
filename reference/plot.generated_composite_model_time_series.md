# Plot a generated composite time series

Produces stacked line plots for each component and the sum for a
`generated_composite_model_time_series` object.

## Usage

``` r
# S3 method for class 'generated_composite_model_time_series'
plot(x, ...)
```

## Arguments

- x:

  A `generated_composite_model_time_series`.

- ...:

  Additional arguments passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly returns `x`.

## Examples

``` r
m2 <- wn(sigma2 = 1) + ar1(phi = 0.8, sigma2 = 0.5)
y2 <- generate(m2, n = 200, seed = 123)
plot(y2)
```
