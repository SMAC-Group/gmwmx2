# Get plot colors for generated time series

Returns a vector of plot colors of length `n`. If `n` exceeds the size
of `gmwmx2_plot_colors`, a larger palette is created via interpolation.

## Usage

``` r
.gmwmx2_get_plot_colors(n)
```

## Arguments

- n:

  Number of colors to return.

## Value

Character vector of hex colors.

## See also

gmwmx2_plot_colors
