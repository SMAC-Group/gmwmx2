# Sum of stochastic models (internal)

Creates a composite model representing a sum of independent
`time_series_model` components. This is an internal constructor used by
the `+` methods.

## Usage

``` r
sum_model(models)
```

## Arguments

- models:

  A list of `time_series_model` objects.

## Value

A `sum_model` object.

## Examples

``` r
mod <- pl(kappa = 0.5, sigma2 = 2) + wn(sigma2 = 1)
print(mod)
#> Stochastic model: sum of 2 processes
#> 
#>  [1] Stationary PowerLaw
#>      Parameters : kappa =   0.5, sigma2 =     2 
#> 
#>  [2] White Noise
#>      Parameters : sigma2 =     1 
#> 
```
