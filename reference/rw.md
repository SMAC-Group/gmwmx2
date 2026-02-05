# Random walk process

Constructs a `time_series_model` for a random walk with innovation
variance `sigma2`. The autocovariance returned is the mean of the
diagonal and super-diagonals of the covariance matrix. The model is
\\X_t = X\_{t-1} + \varepsilon_t, \\ \varepsilon_t
\stackrel{\text{i.i.d.}}{\sim} N(0, \sigma^2)\\.

## Usage

``` r
rw(sigma2 = NULL)
```

## Arguments

- sigma2:

  Innovation variance (\> 0).

## Value

A `time_series_model` object.

## Examples

``` r
mod <- rw(sigma2 = 1)
mod
#> Stochastic process
#>   Model      : Random Walk 
#>   Parameters : sigma2 =     1 
```
