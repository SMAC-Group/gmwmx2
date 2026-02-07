# Plot method for gmwm2_fit

Plots empirical wavelet variance with the fitted theoretical curve.

## Usage

``` r
# S3 method for class 'gmwm2_fit'
plot(
  x,
  show_ci = TRUE,
  col_emp = "black",
  col_theo = "darkorange",
  col_ci = "#e6f7fb",
  lwd = 2,
  pch_emp = 16,
  pch_theo = 21,
  cex_theo = 1.4,
  legend_pos = "auto",
  ...
)
```

## Arguments

- x:

  A `gmwm2_fit` object.

- show_ci:

  Logical; if TRUE and available, show empirical CI bars.

- col_emp:

  Color for empirical WV points/line.

- col_theo:

  Color for theoretical WV line.

- lwd:

  Line width for theoretical curve.

- ...:

  Additional arguments passed to
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html).

- pch:

  Plotting character for empirical points.

## Value

The input object, invisibly.
