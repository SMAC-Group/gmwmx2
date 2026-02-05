# Stationary Power-Law process

Constructs a `time_series_model` representing a stationary power-law
process with parameters `kappa` and `sigma2`. In the frequency domain, a
power-law process is often described by a spectrum \\P(f) = P_0
f^{\kappa}\\ (Bos et al., 2008), where \\f\\ is the frequency, \\P_0\\
is a constant and \\\kappa\\ is the spectral index. Note that we use the
convention that the power spectral density satisfies \\P(f) \propto
\|f\|^{\kappa}\\, where \\\kappa \> -1\\ ensures second-order
stationarity. This corresponds to the alternative notation \\P(f)
\propto \|f\|^{-\alpha}\\ with \\\alpha = -\kappa\\. The autocovariance
used here (Hosking, 1981) is \\\gamma(0) = \sigma^{2}
\frac{\Gamma(1+\kappa)}{\Gamma\left(1+\kappa/2\right)^2}\\, and for \\h
\> 0\\ \\\gamma(h) =\mathrm{cov}(X_t, X\_{t+h}) = \frac{-\kappa/2 + h -
1}{\kappa/2 + h}\\\gamma(h-1)\\.

## Usage

``` r
pl(kappa = NULL, sigma2 = NULL)
```

## Arguments

- kappa:

  Power-law parameter in (-1, 1). Use `inv_trans_kappa_pl` for
  unconstrained optimization.

- sigma2:

  Process variance (\> 0).

## Value

A `time_series_model` object.

## References

Bos MS, Fernandes RMS, Williams SDP, Bastos L (2008). "Fast error
analysis of continuous GPS observations." *Journal of Geodesy*, 82,
157-166.

Hosking JRM (1981). "Fractional differencing." *Biometrika*, 68(1),
165-176.

## Examples

``` r
mod <- pl(kappa = 0.5, sigma2 = 2)
mod
#> Stochastic process
#>   Model      : Stationary PowerLaw 
#>   Parameters : kappa =   0.5, sigma2 =     2 
```
