# Plot a `gnss_ts_ngl` object

Plot a `gnss_ts_ngl` object

## Usage

``` r
# S3 method for class 'gnss_ts_ngl'
plot(x, component = NULL, ...)
```

## Arguments

- x:

  A `gnss_ts_ngl` object.

- component:

  A `string` with value either "N", "E" or "V" that specify which
  component to plot (Northing, Easting or Vertical).

- ...:

  Additional graphical parameters.

## Value

No return value. Plot a `gnss_ts_ngl` object.

## Examples

``` r
station_1LSU <- download_station_ngl("1LSU")
plot(station_1LSU)

plot(station_1LSU, component = "N")

plot(station_1LSU, component = "E")

plot(station_1LSU, component = "V")
```
