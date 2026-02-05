# Flicker noise process

Constructs a `time_series_model` for flicker noise with variance
`sigma2`. The process has spectral density \\S(f) \propto
\frac{1}{\|f\|}\\. Hence, \\\kappa = -1\\ (Bos et al., 2008). The
process is non-stationary and its covariance matrix is assumed to be
given by \$\$ \mathbf C = \sigma^2 \mathbf U^\top \mathbf U, \$\$ where
\\\mathbf U \in \mathbb{R}^{N \times N}\\ is an upper-triangular
Toeplitz matrix with entries \$\$ U\_{i,j} = \begin{cases} h\_{j-i}, & j
\ge i, \\ 0, & j \< i, \end{cases} \qquad i,j = 1, \ldots, N. \$\$ The
coefficients \\\\h_i\\\_{i \ge 0}\\ define a causal linear filter and
are given recursively by \$\$ h_0 = 1, \qquad h_i = \left(i -
\frac{\kappa}{2} - 1\right)\frac{h\_{i-1}}{i}, \quad i \> 0. \$\$

## Usage

``` r
flicker(sigma2 = NULL)
```

## Arguments

- sigma2:

  Innovation variance (\> 0).

## Value

A `time_series_model` object.

## References

Bos MS, Fernandes RMS, Williams SDP, Bastos L (2008). "Fast error
analysis of continuous GPS observations." *Journal of Geodesy*, 82,
157-166.

## Examples

``` r
mod <- flicker(sigma2 = 1)
mod
#> Stochastic process
#>   Model      : Flicker 
#>   Parameters : sigma2 =     1 
```
