# Data generation

This vignette shows how to build stochastic models, combine them with
`+`, generate data, and plot the results.

``` r
library(gmwmx2)
set.seed(123)
```

## Available models

You can build a single stochastic process with any of the following
model constructors:

- [`wn()`](../reference/wn.md) (white noise)
- [`ar1()`](../reference/ar1.md) (AR(1))
- [`pl()`](../reference/pl.md) (power-law)
- [`matern()`](../reference/matern.md) (Matérn)
- [`rw()`](../reference/rw.md) (random walk)
- [`flicker()`](../reference/flicker.md) (flicker / 1/f)

Each constructor returns a `time_series_model` object that can be
plotted or combined with others using `+`.

## Single model generation and plotting

``` r
model_single <- ar1(phi = 0.9, sigma2 = 1)
series_single <- generate(model_single, n = 500)
plot(series_single)
```

![Simulated AR(1) time
series.](data_generation_files/figure-html/single-model-plot-1.png)

The generated object is a `generated_time_series` with a numeric series
in `$series`:

``` r
head(series_single$series)
```

    ## [1] 0.3362423 1.6150311 1.1883829 1.6127387 1.0371249 0.4571655

## Composite models (sum of processes)

Use `+` to build a composite model from multiple stochastic processes.
The result is a `sum_model` that can be passed to
[`generate()`](../reference/generate.md).

``` r
model_comp <- wn(sigma2 = 1) + ar1(phi = 0.8, sigma2 = 0.5) + pl(kappa = 0.3, sigma2 = 2)
series_comp <- generate(model_comp, n = 500)
plot(series_comp)
```

![Simulated composite time series with white noise, AR(1), and power-law
components.](data_generation_files/figure-html/composite-model-plot-1.png)

The composite output is a `generated_composite_model_time_series` with:

- `series` (the total sum)
- `components` (a list of each component series)
- `model` (names of each component)

``` r
names(series_comp)
```

    ## [1] "series"     "components" "n"          "model"      "parameters"

## Combining and plotting

The [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method for
a composite model stacks each component vertically and adds the sum at
the top. For single models,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the time
series with the model name in the title.
