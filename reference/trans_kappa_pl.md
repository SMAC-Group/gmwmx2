# Transform power-law kappa from real line to domain

Maps unconstrained real values to the stationary domain for the
power-law parameter `kappa` using a `tanh` transform.

## Usage

``` r
trans_kappa_pl(x)
```

## Arguments

- x:

  Numeric vector on the real line.

## Value

Numeric vector of transformed values in (-1, 1).
