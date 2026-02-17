#' Loss function for GMWMX without missing value
#'
#' Computes the weighted squared error between empirical wavelet variance
#' and theoretical wavelet variance implied by a model and parameter vector.
#'
#' @param theta Real-valued parameter vector.
#' @param model A `time_series_model` or `sum_model`.
#' @param n Length of autocovariance to compute.
#' @param prep Output from `prepare_optim_layout`.
#' @param wv_obj A `wv::wvar` object.
#' @param omega Optional weighting matrix. If `NULL`, uses inverse CI width.
#' @return Scalar objective value.
#' @keywords internal
loss_fn_gmwmx_no_missing <- function(theta, model, n, prep, wv_obj, quantities_D, omega = NULL) {
  # compute autocovariance from theta
  autocov_vec <- get_autocovariance(object = model, n = n, theta = theta, prep = prep)

  # retreat autocovariance to take account of the fact that we are using the residuals of a regression, not the original series
  vec_mean_autocov_eps_hat <- compute_all_mean_diag_fast_w_linear_interp_only_required_cpp(
    mat_D_q_term_1 = quantities_D$mat_D_q_term_1,
    mat_D_q_term_2 = quantities_D$mat_D_q_term_2,
    sum_on_sub_diag_of_D = quantities_D$sum_on_sub_diag_of_D,
    vec_autocov = autocov_vec, approx_type = "3"
  )

  # compute wv from autocovariance
  theo_wv <- autocovariance_to_wv(vec_mean_autocov_eps_hat, tau = wv_obj$scales)

  # define omega
  nu_hat <- wv_obj$variance
  if (is.null(omega)) {
    omega <- diag(1 / (wv_obj$ci_low - wv_obj$ci_high)^2)
  }

  difference <- nu_hat - theo_wv
  objective <- as.numeric(t(difference) %*% omega %*% difference)


  return(objective)
  # compute loss
}




#' Loss function for GMWMX with missing value
#'
#' Computes the weighted squared error between empirical wavelet variance
#' and theoretical wavelet variance implied by a model and parameter vector.
#'
#' @param theta Real-valued parameter vector.
#' @param model A `time_series_model` or `sum_model`.
#' @param n Length of autocovariance to compute.
#' @param prep Output from `prepare_optim_layout`.
#' @param wv_obj A `wv::wvar` object.
#' @param omega Optional weighting matrix. If `NULL`, uses inverse CI width.
#' @return Scalar objective value.
#' @keywords internal
loss_fn_gmwmx_with_missing <- function(theta, model, n, prep, wv_obj, quantities_D, vec_autocov_omega, pstar_hat, omega = NULL) {
  # compute autocovariance from theta
  autocov_vec <- get_autocovariance(object = model, n = n, theta = theta, prep = prep)

  # retreat autocovariance to take account of the fact that we are using the residuals of a regression, not the original series
  vec_mean_autocov_eps_hat <- compute_all_mean_diag_fast_w_linear_interp_only_required_cpp(
    mat_D_q_term_1 = quantities_D$mat_D_q_term_1,
    mat_D_q_term_2 = quantities_D$mat_D_q_term_2,
    sum_on_sub_diag_of_D = quantities_D$sum_on_sub_diag_of_D,
    vec_autocov = autocov_vec, approx_type = "3"
  )

  # retreat autocovariance of epsilon hat with the missing data mechanism
  vec_mean_per_diag_w_missing <- vec_mean_autocov_eps_hat * (vec_autocov_omega + pstar_hat^2)

  # compute wv from autocovariance
  theo_wv <- autocovariance_to_wv(vec_mean_per_diag_w_missing, tau = wv_obj$scales)

  # define omega
  nu_hat <- wv_obj$variance
  if (is.null(omega)) {
    omega <- diag(1 / (wv_obj$ci_low - wv_obj$ci_high)^2)
  }

  # compute difference between empirical and theoretical wv, and compute objective
  difference <- nu_hat - theo_wv

  # compute objective
  objective <- as.numeric(t(difference) %*% omega %*% difference)



  return(objective)
}




