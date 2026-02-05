# Changelog

## gmwmx2 version 0.0.4

CRAN release: 2026-01-07

- updated functions `download_all_stations_ngl`,
  `download_estimated_velocities_ngl` and `download_station_ngl` so that
  it exit gracefully if NGL server is not reachable.

## gmwmx2 version 0.0.3

CRAN release: 2025-08-19

- Solve minor bug when creating design matrix if jumps or earthquakes
  are indicated after the last observation.

## gmwmx2 version 0.0.2

CRAN release: 2025-04-10

- updated functions `download_all_stations_ngl`,
  `download_estimated_velocities_ngl` and `download_station_ngl` with
  new address of the NGL: <https://geodesy.unr.edu> and using `httr2`.

## gmwmx2 version 0.0.1

CRAN release: 2025-03-10

first version
