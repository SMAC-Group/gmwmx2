# White noise process

Constructs a `time_series_model` for white noise with variance `sigma2`.
The process is defined as \\X_t \stackrel{\text{i.i.d.}}{\sim} N(0,
\sigma^2)\\ with autocovariance \\\gamma(h) = \mathrm{cov}(X_t,
X\_{t+h}) = \sigma^2 \mathbf{1}\\h=0\\\\

## Usage

``` r
wn(sigma2 = NULL)
```

## Arguments

- sigma2:

  Innovation variance (\> 0).

## Value

A `time_series_model` object.

## Examples

``` r
mod <- wn(sigma2 = 1)
mod
#> Stochastic process
#>   Model      : White Noise 
#>   Parameters : sigma2 =     1 
```
