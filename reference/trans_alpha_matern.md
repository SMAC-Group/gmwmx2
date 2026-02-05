# Transform Matérn alpha from real line to domain

Ensures Matérn smoothness parameter `alpha` is \> 1/2 by mapping real
values with `exp(x) + 1/2`.

## Usage

``` r
trans_alpha_matern(x)
```

## Arguments

- x:

  Numeric vector on the real line.

## Value

Numeric vector with values \> 1/2.
