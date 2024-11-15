#'
#' Estimated northward and eastward velocity and their standard deviation using the GMWMX estimator
#'
#' @description Estimated northward and eastward velocity and standard deviation for a subset of 1202 GNSS station with more than 10 years of daily data.
#'
#' @format A data frame with 1202 rows and 12 variables:
#' \describe{
#'   \item{station_name}{Name of the GNSS station.}
#'   \item{estimated_trend_N}{Estimated northward velocity trend (in meters per day).}
#'   \item{std_estimated_trend_N}{Standard deviation of the estimated northward velocity trend.}
#'   \item{estimated_trend_E}{Estimated eastward velocity trend (in meters per day).}
#'   \item{std_estimated_trend_E}{Standard deviation of the estimated eastward velocity trend.}
#'   \item{length_signal}{Length of the signal (in days).}
#'   \item{estimated_trend_N_scaled}{Scaled estimated northward velocity trend (multiplying by 365.25 for yearly values).}
#'   \item{std_estimated_trend_N_scaled}{Scaled standard deviation of the estimated northward velocity trend.}
#'   \item{estimated_trend_E_scaled}{Scaled estimated eastward velocity trend (multiplying by 365.25 for yearly values).}
#'   \item{std_estimated_trend_E_scaled}{Scaled standard deviation of the estimated eastward velocity trend.}
#'   \item{latitude}{Latitude of the GNSS station.}
#'   \item{longitude}{Longitude of the GNSS station.}
#' }
#'
"df_estimated_velocities_gmwmx"
#'
