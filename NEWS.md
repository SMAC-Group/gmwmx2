# gmwmx2 version 0.0.5

- Added a composable stochastic-model interface. Models are now built with
  `wn()`, `ar1()`, `rw()`, `pl()`, `flicker()`, and `matern()`, and independent
  components can be combined with `+`.
- Added `gmwm2()` for fitting stochastic models directly from a time series by
  matching empirical and theoretical wavelet variances. The returned
  `gmwm2_fit` objects include print and plot methods.
- Reworked `gmwmx2()` as an S3 generic. It now supports both the GNSS workflow
  for `gnss_ts_ngl` objects and a generic regression interface
  `gmwmx2(X, y, model, ...)` for arbitrary design matrices and response
  vectors, including responses with missing observations.
- Changed GNSS model specification from string labels such as
  `stochastic_model = "wn + pl"` to model objects such as
  `model = wn() + pl()` or `model = wn() + flicker()`. GNSS fit objects now use
  class `gmwmx2_fit_gnss_ts_ngl` and have updated print and plot methods.
- Added simulation utilities through `generate()` methods for single
  stochastic models, composite stochastic models, and two-state Markov
  missingness models created with `markov_two_states()`.
- Added infrastructure for model parameter initialization, constrained
  parameter transformations, autocovariance evaluation, theoretical wavelet
  variance computation, and model-implied covariance matrices.
- Added random-walk covariance C++ routines and registered them through
  `RcppExports`.
- Expanded documentation and vignettes for data generation, composite
  stochastic model estimation, generic dependent-error regression, regression
  with missing observations, and GNSS time-series estimation.
- Updated package metadata for the new modeling stack, including the `longmemo`
  dependency and pkgdown website URL.

# gmwmx2 version 0.0.4

- updated functions `download_all_stations_ngl`, `download_estimated_velocities_ngl` and `download_station_ngl` so that it exit gracefully if NGL server is not reachable.

# gmwmx2 version 0.0.3

- Solve minor bug when creating design matrix if jumps or earthquakes are indicated after the last observation.

# gmwmx2 version 0.0.2

- updated functions `download_all_stations_ngl`, `download_estimated_velocities_ngl` and `download_station_ngl` with new address of the NGL: https://geodesy.unr.edu and using `httr2`.

# gmwmx2 version 0.0.1

first version
