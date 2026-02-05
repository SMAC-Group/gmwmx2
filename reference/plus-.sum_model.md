# Add to a sum_model

Add to a sum_model

## Usage

``` r
# S3 method for class 'sum_model'
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
m3 <- pl(kappa = 0.3, sigma2 = 2)
model <- (m1 + m2) + m3
```
