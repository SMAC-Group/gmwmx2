# gmwmx2 version 0.0.5

- Added a new composable stochastic-model framework with `time_series_model` / `sum_model` objects and component constructors: `wn()`, `ar1()`, `rw()`, `pl()`, `flicker()`, and `matern()`.
- Added `gmwm2()` to fit stochastic models directly from wavelet variance, with new helpers for parameter layout, transformations, theoretical WV computation, and model printing/plotting.
- Added `gmwmx2_new()` as a generic GMWMX interface for arbitrary `X` and `y` (in addition to GNSS-specific workflows), including explicit paths for complete and incomplete series (`gmwmx2_new_no_missing()` / `gmwmx2_new_with_missing()`).
- Extended simulation/generation utilities with `generate()` methods for single and composite models plus missingness support (`markov_two_states()` and related plotting methods).
- Added random-walk covariance C++ back-end routines (`get_sigma_mat_rw()`, `get_mean_diagonal_super_diagonals_cov_mat_rw_cpp()`) and integrated them through `RcppExports`.
- Expanded documentation with many new man pages and three new vignettes: `data_generation.Rmd`, `estimate_composite_stochastic_models.Rmd`, and `gmwmx2_new.Rmd`.
- Updated package metadata and dependencies for the new modeling stack (notably adding `longmemo` and pkgdown URL metadata updates).

# gmwmx2 version 0.0.4

- updated functions `download_all_stations_ngl`, `download_estimated_velocities_ngl` and `download_station_ngl` so that it exit gracefully if NGL server is not reachable.

# gmwmx2 version 0.0.3

- Solve minor bug when creating design matrix if jumps or earthquakes are indicated after the last observation.

# gmwmx2 version 0.0.2

- updated functions `download_all_stations_ngl`, `download_estimated_velocities_ngl` and `download_station_ngl` with new address of the NGL: https://geodesy.unr.edu and using `httr2`.

# gmwmx2 version 0.0.1

first version
