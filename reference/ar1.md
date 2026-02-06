# AR(1) process

Constructs a `time_series_model` for a stationary AR(1) process with
parameter `phi` and innovation variance `sigma2`. The model is \\X_t =
\phi X\_{t-1} + \varepsilon_t, \\ \varepsilon_t
\stackrel{\text{i.i.d.}}{\sim} N(0, \sigma^2)\\. The autocovariance is
\\ \gamma(h) = \mathrm{cov}(X_t, X\_{t+h}) = \frac{\sigma^2}{1 -
\phi^2}\\\phi^{\lvert h \rvert} \\.

## Usage

``` r
ar1(phi = NULL, sigma2 = NULL)
```

## Arguments

- phi:

  AR(1) coefficient in (-1, 1).

- sigma2:

  Innovation variance (\> 0).

## Value

A `time_series_model` object.

## Examples

``` r
mod <- ar1(phi = 0.8, sigma2 = 1)
mod
#> Stochastic process
#>   Model      : AR(1) 
#>   Parameters : phi =   0.8, sigma2 =     1 
```
