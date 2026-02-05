# Add two model objects

Combines `time_series_model` and/or `sum_model` into a `sum_model`.

## Usage

``` r
# S3 method for class 'time_series_model'
e1 + e2
```

## Arguments

- e1:

  Left operand.

- e2:

  Right operand.

## Value

A `sum_model`.

## Examples

``` r
m1 <- wn(sigma2 = 1)
m2 <- ar1(phi = 0.8, sigma2 = 0.5)
model <- m1 + m2
model
#> Stochastic model: sum of 2 processes
#> 
#>  [1] White Noise
#>      Parameters : sigma2 =     1 
#> 
#>  [2] AR(1)
#>      Parameters : phi =   0.8, sigma2 =   0.5 
#> 
```
