# Print method for gmwm2_fit

Print method for gmwm2_fit

## Usage

``` r
# S3 method for class 'gmwm2_fit'
print(x, digits = 4, ...)
```

## Arguments

- x:

  A `gmwm2_fit` object.

- digits:

  Significant digits for printing.

- ...:

  Unused.

## Value

The input object, invisibly.

## Examples

``` r
model <- wn(sigma2 = 1) + ar1(phi = 0.8, sigma2 = 0.5)
x <- generate(model, n = 1000, seed = 123)
plot(x)

fit <- gmwm2(x, model = wn()+ar1())
fit
#> GMWM fit
#> 
#> Stochastic model
#>   Sum of 2 processes
#> 
#>   [1] White Noise
#>       Parameters : sigma2
#> 
#>   [2] AR(1)
#>       Parameters : phi, sigma2
#> 
#> Initial parameters
#>   1) White Noise: sigma2 = 2.424
#>   2) AR(1): phi = 0.8972, sigma2 = 2.424
#> 
#> Estimated parameters
#>   1) White Noise: sigma2 = 0.9147
#>   2) AR(1): phi = 0.7744, sigma2 = 0.5982
#> 
#> Optimization
#>   Convergence : converged (code 0)
#>   Iterations  : 40
#>   Loss        : 0.03692
```
