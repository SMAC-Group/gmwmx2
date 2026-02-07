# Package index

## Composite stochastic models and simulation

- [`wn()`](wn.md) : White noise process
- [`ar1()`](ar1.md) : AR(1) process
- [`pl()`](pl.md) : Stationary Power-Law process
- [`matern()`](matern.md) : Matern process
- [`rw()`](rw.md) : Random walk process
- [`flicker()`](flicker.md) : Flicker noise process
- [`` `+`( ``*`<time_series_model>`*`)`](plus-.time_series_model.md) :
  Add two model objects
- [`` `+`( ``*`<sum_model>`*`)`](plus-.sum_model.md) : Add to a
  sum_model
- [`generate()`](generate.md) : Generate a time series from a model
- [`plot(`*`<generated_time_series>`*`)`](plot.generated_time_series.md)
  : Plot a generated time series
- [`plot(`*`<generated_composite_model_time_series>`*`)`](plot.generated_composite_model_time_series.md)
  : Plot a generated composite time series

## Composite stochastic models estimation

- [`gmwm2()`](gmwm2.md) : Estimate stochastic model parameters by
  matching wavelet variance
- [`print(`*`<gmwm2_fit>`*`)`](print.gmwm2_fit.md) : Print method for
  gmwm2_fit
- [`plot(`*`<gmwm2_fit>`*`)`](plot.gmwm2_fit.md) : Plot method for
  gmwm2_fit

## Load and plot NGL data

- [`download_station_ngl()`](download_station_ngl.md) : Download GNSS
  position time series and steps reference from the Nevada Geodetic
  Laboratory with IGS14 or IGS20 reference frame.

- [`download_all_stations_ngl()`](download_all_stations_ngl.md) :
  Download all stations name and location from the Nevada Geodetic
  Laboratory

- [`download_estimated_velocities_ngl()`](download_estimated_velocities_ngl.md)
  : Download estimated velocities provided by the Nevada Geodetic
  Laboratory for all stations.

- [`plot(`*`<gnss_ts_ngl>`*`)`](plot.gnss_ts_ngl.md) :

  Plot a `gnss_ts_ngl` object

## Estimate a model

- [`gmwmx2()`](gmwmx2.md) :

  Estimate a trajectory model for a `gnss_ts_ngl` object considering a
  white noise plus colored noise as the stochastic model for the
  residuals and model missingness with a Markov process using the GMWMX
  estimator.

- [`summary(`*`<fit_gnss_ts_ngl>`*`)`](summary.fit_gnss_ts_ngl.md) :

  Extract estimated parameters from a `fit_gnss_ts_ngl`

- [`plot(`*`<fit_gnss_ts_ngl>`*`)`](plot.fit_gnss_ts_ngl.md) :

  Plot a `fit_gnss_ts_ngl` object

## Data

- [`df_estimated_velocities_gmwmx`](df_estimated_velocities_gmwmx.md) :
  Estimated northward and eastward velocity and their standard deviation
  using the GMWMX estimator
