# Plot a generated time series

Produces a single line plot for a `generated_time_series` object.

## Usage

``` r
# S3 method for class 'generated_time_series'
plot(x, ...)
```

## Arguments

- x:

  A `generated_time_series`.

- ...:

  Additional arguments passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

Invisibly returns `x`.

## Examples

``` r
m1 <- wn(sigma2 = 1)
y1 <- generate(m1, n = 200, seed = 123)
plot(y1)
```
