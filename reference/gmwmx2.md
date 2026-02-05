# Estimate a trajectory model for a `gnss_ts_ngl` object considering a white noise plus colored noise as the stochastic model for the residuals and model missingness with a Markov process using the GMWMX estimator.

Estimate a trajectory model for a `gnss_ts_ngl` object considering a
white noise plus colored noise as the stochastic model for the residuals
and model missingness with a Markov process using the GMWMX estimator.

## Usage

``` r
gmwmx2(
  x,
  n_seasonal = 2,
  vec_earthquakes_relaxation_time = NULL,
  component = "N",
  toeplitz_approx_var_cov_wv = TRUE,
  stochastic_model = "wn + fl"
)
```

## Arguments

- x:

  A `gnss_ts_ngl` object.

- n_seasonal:

  An `integer` specifying the number of seasonal signals in the time
  series. "1" specify only one annual periodic signal and "2"specify an
  annual and a semiannual periodic signal.

- vec_earthquakes_relaxation_time:

  A `vector` specifying the relaxation time for each earthquakes
  indicated for the time series.

- component:

  A `string` with value either "N", "E" or "V" that specify which
  component to estimate (Northing, Easting or Vertical).

- toeplitz_approx_var_cov_wv:

  A `boolean` that specify if the variance of the wavelet variance
  should be computed based on a toeplitz approximation of the variance
  covariance matrix of the residuals.

- stochastic_model:

  A `string` that specify the stochastic model considered for the
  residuals. Either "wn + fl" for white noise and flicker/pink noise or
  "wn + pl" for white noise and stationary power-law noise.

## Examples

``` r
x <- download_station_ngl("CHML")
fit <- gmwmx2(x, n_seasonal = 2, component = "N")
```
