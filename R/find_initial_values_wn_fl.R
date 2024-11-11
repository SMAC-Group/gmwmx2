#' @importFrom wv wn_to_wv
find_initial_values_wn_fl <- function(wv_emp) {
  sigma2_wn_start <- 2 * wv_emp$variance[1]
  wv_after_removing_wn <- wv_emp$variance - wv::wn_to_wv(sigma2 = sigma2_wn_start, tau = wv_emp$scales)
  n <- tail(wv_emp$scales, 1)
  mean_autocov_flicker <- vec_mean_autocov_powerlaw(kappa = -1, n)
  unit_wv_flicker <- autocovariance_to_wv(acf = mean_autocov_flicker, tau = wv_emp$scales)
  sigma2_flicker_start <- tail(wv_after_removing_wn, 1) / tail(unit_wv_flicker, 1)
  return(abs(c(sigma2_wn_start, sigma2_flicker_start)))
}
