# Plot a `fit_gnss_ts_ngl` object

Plot a `fit_gnss_ts_ngl` object

## Usage

``` r
# S3 method for class 'fit_gnss_ts_ngl'
plot(x, ...)
```

## Arguments

- x:

  A `fit_gnss_ts_ngl` object.

- ...:

  Additional graphical parameters.

## Value

No return value. Plot a `fit_gnss_ts_ngl` object.

## Examples

``` r
x <- download_station_ngl("0AMB")
fit_N <- gmwmx2(x, n_seasonal = 2, component = "N")
plot(fit_N)
```
