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

## Examples

``` r
.gmwmx2_get_plot_colors(3)
#> Error in .gmwmx2_get_plot_colors(3): could not find function ".gmwmx2_get_plot_colors"
.gmwmx2_get_plot_colors(20)
#> Error in .gmwmx2_get_plot_colors(20): could not find function ".gmwmx2_get_plot_colors"
```
