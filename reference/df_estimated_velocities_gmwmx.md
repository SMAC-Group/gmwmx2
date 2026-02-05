# Estimated northward and eastward velocity and their standard deviation using the GMWMX estimator

Estimated northward and eastward velocity and standard deviation for a
subset of 1202 GNSS station with more than 10 years of daily data.

## Usage

``` r
df_estimated_velocities_gmwmx
```

## Format

A data frame with 1202 rows and 12 variables:

- station_name:

  Name of the GNSS station.

- estimated_trend_N:

  Estimated northward velocity trend (in meters per day).

- std_estimated_trend_N:

  Standard deviation of the estimated northward velocity trend.

- estimated_trend_E:

  Estimated eastward velocity trend (in meters per day).

- std_estimated_trend_E:

  Standard deviation of the estimated eastward velocity trend.

- length_signal:

  Length of the signal (in days).

- estimated_trend_N_scaled:

  Scaled estimated northward velocity trend (multiplying by 365.25 for
  yearly values).

- std_estimated_trend_N_scaled:

  Scaled standard deviation of the estimated northward velocity trend.

- estimated_trend_E_scaled:

  Scaled estimated eastward velocity trend (multiplying by 365.25 for
  yearly values).

- std_estimated_trend_E_scaled:

  Scaled standard deviation of the estimated eastward velocity trend.

- latitude:

  Latitude of the GNSS station.

- longitude:

  Longitude of the GNSS station.
