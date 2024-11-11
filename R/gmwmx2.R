create_X_matrix = function (all_mjd_index,
                            jumps,
                            n_seasonal,
                            vec_earthquakes_index_mjd,
                            vec_earthquakes_relaxation_time) {

  # x = download_station_ngl("CHML")
  # # create all jumps by combining jumps due to equipment change and jumps due to earthquakes
  # jumps = c(x$df_equipment_software_changes$modified_julian_date,
  #           x$df_earthquakes$modified_julian_date )
  #
  # # if multiple jumps  to prevent not invertible matrix
  # jumps = unique(jumps)
  # vec_earthquakes_index_mjd = c(x$df_earthquakes$modified_julian_date)
  # # if multiple earthquakes  to prevent not invertible matrix
  # vec_earthquakes_index_mjd = unique(vec_earthquakes_index_mjd)
  # vec_earthquakes_relaxation_time = NULL
  #
  # all_mjd_index = seq(head(x$df_position$modified_julian_day, 1), tail(x$df_position$modified_julian_day, 1))



  # ensure that vec_earthquakes_relaxation_time is the same length as vec_earthquakes_index_mjd if not NULL
  if(!is.null(vec_earthquakes_relaxation_time)){
    if(length(vec_earthquakes_index_mjd) != length(vec_earthquakes_relaxation_time)){
      stop("Vector of relaxation time not of same length as vector of MJD index of earthquakes")
    }
  }


  # set to default 365.25 days if no earthquakes relaxation time provided
  if(is.null(vec_earthquakes_relaxation_time)){
    vec_earthquakes_relaxation_time = rep(365.25, length(vec_earthquakes_index_mjd))
  }


  # create empty matrix number of columns bias+trend+2*nbr sinusoidal + nbr jumps + number earthquakes
  X = matrix(0, nrow = length(all_mjd_index), ncol = 2 + 2 * n_seasonal +
               length(jumps) + length(vec_earthquakes_index_mjd))

  # add bias, intercept
  X[, 1] = 1

  # add component for trend
  reference_time =  0.5*(all_mjd_index[1]+tail(all_mjd_index,1))
  X[,2] = all_mjd_index - reference_time


  # add seasonal
  if (n_seasonal > 0) {
    for (i in 1:n_seasonal) {
      omega_i = (i/365.25) * 2 * pi
      X[, 2 + (i - 1) * 2 + 1] = sin((all_mjd_index) *  omega_i)
      X[, 2 + (i - 1) * 2 + 2] = cos((all_mjd_index) *  omega_i)
    }
  }

  # add offsets
  if (!is.null(jumps)) {
    for (i in 1:length(jumps)) {
      it = min(which(all_mjd_index > jumps[i] - 1e-06))
      X[, 2 + 2 * n_seasonal + i] = c(rep(0, it - 1), rep(1,
                                                          length(all_mjd_index) - it + 1))
    }
  }


  # exponential decay function for post seismic relaxation
  if(!is.null(vec_earthquakes_index_mjd)){
    for(i in seq_along(vec_earthquakes_index_mjd)){
      tau_i = vec_earthquakes_relaxation_time[i]
      earthquake_mjd_i = vec_earthquakes_index_mjd[i]
      # create vector
      decay_values <- ifelse(all_mjd_index > earthquake_mjd_i,
                             1 - exp(-(all_mjd_index - earthquake_mjd_i) / tau_i),
                             0)

      # create column in matrix
      X[, 2 + 2 * n_seasonal + length(jumps) + i] = decay_values
    }
  }


  # # Slow slip events (tanh)
  # if(!is.null(vec_tanh_mid_point)){
  #   for(i in seq_along(vec_tanh_mid_point)){
  #     # i = 1
  #     t_k = vec_tanh_mid_point[i]
  #     T_k = vec_tanh_length[i]
  #     X[, 2 + 2 * n_seasonal + length(jumps) + i] = 0.5 * (tanh((t_nogap - t_k)/T_k) -1)
  #   }
  # }
  rownames(X) = all_mjd_index
  return(X)

}


# define optimization function
objective_function_wn_flicker_w_missing <- function(theta, wv_obj, n, quantities_D, approx_type, vec_autocov_omega, pstar_hat, no_missing=T, omega = NULL) {
  theta_t <- vector(mode = "numeric", length = 2)
  theta_t[1] <- exp(theta[1]) # sigma2 wn
  theta_t[2] <- exp(theta[2]) # sigma2 fl

  vec_mean_autocov <- vec_mean_autocov_powerlaw(kappa = -1, n) *  theta_t[2]
  vec_mean_autocov[1] <- vec_mean_autocov[1] +   theta_t[1]

  # approx with linear interpolation on errors
  vec_mean_autocov_eps_hat <- compute_all_mean_diag_fast_w_linear_interp_only_required_cpp(
    mat_D_q_term_1 = quantities_D$mat_D_q_term_1,
    mat_D_q_term_2 = quantities_D$mat_D_q_term_2,
    sum_on_sub_diag_of_D = quantities_D$sum_on_sub_diag_of_D,
    vec_autocov = vec_mean_autocov, approx_type = approx_type
  )
  if(!no_missing){
    vec_mean_per_diag_w_missing = vec_mean_autocov_eps_hat * ( vec_autocov_omega + pstar_hat^2)
    theo_wv <- autocovariance_to_wv(vec_mean_per_diag_w_missing, tau = wv_obj$scales)

  }else{
    theo_wv <- autocovariance_to_wv(vec_mean_autocov_eps_hat, tau = wv_obj$scales)

  }
  nu_hat <- wv_obj$variance
  if (is.null(omega)) {
    omega <- diag(1 / (wv_obj$ci_low - wv_obj$ci_high)^2)
  }

  difference <- nu_hat - theo_wv
  objective <- t(difference) %*% omega %*% difference
  return(objective)
}



#' Estimate a linear model with White noise and Flicker noise for the residuals in presence of missing data using the GMWMX estimator.
#' @param x A \code{gnss_ts_ngl} object.
#' @param n_seasonal An \code{integer} specifying the number of seasonal signals in the time series. "1" specify only one annual periodic signal and "2"specify an annual and a semiannual periodic signal.
#' @param vec_earthquakes_relaxation_time A \code{vecor} specifying the relaxation time for each earthquakes indicated for the time series.
#' @param component A \code{string} with value either "N", "E" or "V" that specify which component to estimate (Northing, Easting or Vertical).
#' @importFrom wv wvar
#' @importFrom dplyr between
#' @importFrom Matrix solve
#' @importFrom MASS ginv
#' @importFrom stats .lm.fit optim toeplitz
#' @export
gmwmx2 = function(x, n_seasonal=2, vec_earthquakes_relaxation_time = NULL, component ="N"){

  # x = download_station_ngl("CHML")
  # plot(x, component = "N")
  # vec_earthquakes_relaxation_time = NULL
  # component ="N"
  # n_seasonal=2




  # check that component is either N, E or V
  if(!component %in% c("N", "E", "V")){
    stop("Specified component should be either 'N', 'E' or 'V'")
  }

  # create full index
  all_mjd_index = seq(head(x$df_position$modified_julian_day,1), tail(x$df_position$modified_julian_day, 1),by =1)

  # create all jumps by combining jumps due to equipment change and jumps due to earthquakes
  jumps = c(x$df_equipment_software_changes$modified_julian_date,
            x$df_earthquakes$modified_julian_date )

  # if multiple jumps  to prevent not invertible matrix
  jumps = unique(jumps)
  vec_earthquakes_index_mjd = c(x$df_earthquakes$modified_julian_date)

  # if multiple earthquakes  to prevent not invertible matrix
  vec_earthquakes_index_mjd = unique(vec_earthquakes_index_mjd)

  if(length(jumps) == 0){
    jumps=NULL
  }

  # create matrix
  X = create_X_matrix(all_mjd_index = all_mjd_index,
                      jumps = jumps,
                      n_seasonal = n_seasonal,
                      vec_earthquakes_index_mjd = vec_earthquakes_index_mjd,
                      vec_earthquakes_relaxation_time = vec_earthquakes_relaxation_time )

  # obtain X_sub
  id_X_sub = which(rownames(X) %in% x$df_position$modified_julian_day)
  X_sub = X[id_X_sub,]

  # Extract Y given specified component
  if(component == "N"){
    y = x$df_position$northings_fractional_portion
  }else if (component == "E"){
    y = x$df_position$eastings_fractional_portion
  }else if(component == "V"){
    y = x$df_position$vertical_fractional_portion
  }

  # obtain beta hat
  beta_hat <- .lm.fit(y = y , x = X_sub)$coefficients

  # # plot signal and estimated model
  # plot(x=X_sub[,2], y=y, type="l")
  # lines(x=X_sub[, 2], y= X_sub%*% beta_hat, col="red")
  # # obtain observed residuals
  eps_hat_sub <- y - X_sub %*% beta_hat

  # create
  eps_hat_filled <- vector(mode = "numeric", length = length(all_mjd_index))

  # fill in observed residuals
  eps_hat_filled[id_X_sub] = eps_hat_sub

  # compute empirical wv variance on the "filled" vector of residuals
  wv_emp_eps_hat_filled <- wv::wvar(eps_hat_filled)

  # plot(wv_emp_eps_hat_filled)
  # define vector of observed index
  vec_omega = as.numeric(all_mjd_index %in% x$df_position$modified_julian_day )

  # estimate parameter of Markov process missingness process using the MLE
  p_hat <- estimate_p1_p2_mle_cpp(vec_omega)

  pstar_hat <- p_hat[2] / (p_hat[1] + p_hat[2])

  # get vec autocovariance theo omega
  vec_autocov_omega <- create_vec_theo_autocov_omega_cpp(p1 = p_hat[1], p2 = p_hat[2], length(all_mjd_index))

  XtX = t(X) %*% X
  inv_XtX = Matrix::solve(XtX) # to check later, why does when we have multiple earthquake XtX becomes not invertible
  H <- X %*% inv_XtX %*% t(X)
  D <- diag(length(all_mjd_index)) - H

  # precompute quantities on D
  quantities_D <- pre_compute_quantities_on_D_only_required_smarter_cpp(D, approx_type = "3")

  # define missingness case
  if (all(vec_omega == 1)) {
    no_missing <- T
  } else {
    no_missing <- F
  }


  # define gamma init
  gamma_init = c(log(8e-7), log(8e-9)) # arbitrary values around the ones observed in practice for these kind of data

  # fit gmwmx on empirical wv
  res_gmwmx <- optim(
    par = gamma_init,
    fn = objective_function_wn_flicker_w_missing,
    wv_obj = wv_emp_eps_hat_filled,
    n = length(vec_omega),
    quantities_D = quantities_D,
    approx_type = "3",
    vec_autocov_omega = vec_autocov_omega,
    pstar_hat = pstar_hat,
    no_missing = no_missing
  )

  # res_gmwmx
  gamma_hat_1 = exp(res_gmwmx$par)

  # construct sigma matrix of white noise + ficker
  var_cov_mat_wn = gamma_hat_1[1] * diag(length(vec_omega))
  var_cov_mat_flicker = var_cov_powerlaw_cpp(sigma2 = gamma_hat_1[2], kappa = -1, n = length(vec_omega))
  var_cov_mat_epsilon = var_cov_mat_wn+var_cov_mat_flicker

  # get var cov missingness process
  var_cov_omega = toeplitz(as.vector(vec_autocov_omega))

  # define variance covariance of residuals with missing
  var_cov_eps_hat_w_missing = var_cov_mat_epsilon * ( var_cov_omega + pstar_hat^2)

  # compute variance of wavelet variance computed on epsilon (negleting the fact that we compute the wv on the "estimated" residuals)
  var_cov_nu_hat = compute_cov_wv_cpp_approx_faster(Sigma_X = var_cov_eps_hat_w_missing)

  # get inverse
  inv_var_cov_nu_hat =  Matrix::solve(var_cov_nu_hat)

  # re estimate gamma with optimal weighting matrix
  res_gmwmx_2 <- optim(
    par = res_gmwmx$par,
    fn = objective_function_wn_flicker_w_missing,
    wv_obj = wv_emp_eps_hat_filled,
    n = length(vec_omega),
    quantities_D = quantities_D,
    approx_type = "3",
    vec_autocov_omega = vec_autocov_omega,
    omega = inv_var_cov_nu_hat,
    pstar_hat = pstar_hat,
    no_missing = no_missing
  )

  gamma_hat_2 = exp(res_gmwmx_2$par)

  # get theo wv with approximation obtained by last fit
  vec_mean_autocov <- vec_mean_autocov_powerlaw(kappa = -1, length(vec_omega)) *  gamma_hat_2[2]
  vec_mean_autocov[1] <- vec_mean_autocov[1] +   gamma_hat_2[1]

  # approx with linear interpolation on errors
  vec_mean_autocov_eps_hat <- compute_all_mean_diag_fast_w_linear_interp_only_required_cpp(
    mat_D_q_term_1 = quantities_D$mat_D_q_term_1,
    mat_D_q_term_2 = quantities_D$mat_D_q_term_2,
    sum_on_sub_diag_of_D = quantities_D$sum_on_sub_diag_of_D,
    vec_autocov = vec_mean_autocov, approx_type = "3"
  )
  if(!no_missing){
    vec_mean_per_diag_w_missing = vec_mean_autocov_eps_hat * ( vec_autocov_omega + pstar_hat^2)
    theo_wv <- autocovariance_to_wv(vec_mean_per_diag_w_missing, tau = wv_emp_eps_hat_filled$scales)

  }else{
    theo_wv <- autocovariance_to_wv(vec_mean_autocov_eps_hat, tau = wv_emp_eps_hat_filled$scales)

  }



  # Compute variance covariance of beta hat
  if (no_missing) {
    # XtX_inv <- Matrix::solve(t(X) %*% X)
    var_cov_beta_hat = inv_XtX %*% t(X) %*% var_cov_mat_epsilon  %*% X %*% inv_XtX
  } else {
    # XtX_inv <- Matrix::solve(t(X) %*% X)
    var_cov_beta_hat <- pstar_hat^(-2) * inv_XtX %*% t(X) %*% ((var_cov_omega + pstar_hat^2) * var_cov_mat_epsilon) %*% X %*% inv_XtX
  }

  std_beta_hat_gmwmx_3 <- sqrt(diag(var_cov_beta_hat))

  ret = list("beta_hat" = beta_hat,
             "std_beta_hat" = std_beta_hat_gmwmx_3,
             "gamma_hat" = gamma_hat_2,
             "vartheta_hat" = p_hat,
             "component" = component,
             "design_matrix_X" = X_sub,
             "y" = y,
             "empirical_wvar" = wv_emp_eps_hat_filled,
             "theoretical_wvar" =theo_wv,
             "df_position" = x$df_position,
             "df_earthquakes" = x$df_earthquakes,
             "df_equipment_software_changes" = x$df_equipment_software_changes)

  class(ret) = "fit_gnss_ts_ngl"

  return(ret)


}



#' Extract estimated parameters from a \code{fit_gnss_ts_ngl}
#' @param object A \code{fit_gnss_ts_ngl} object.
#' @param ... Additional parameters.
#' @export
summary.fit_gnss_ts_ngl <- function(object, ...) {
  # Print header
  cat("Summary of Estimated Model\n")
  cat("-------------------------------------------------------------\n")
  cat("Functional parameters\n")
  cat("-------------------------------------------------------------\n")

  cat(" Estimate       Std_Deviation     95% CI Lower    95% CI Upper\n")
  cat("-------------------------------------------------------------\n")

  for (i in seq_along(object$beta_hat)) {
    lower_ci <- object$beta_hat[i] - 1.96 * object$std_beta_hat[i]
    upper_ci <- object$beta_hat[i] + 1.96 * object$std_beta_hat[i]

    # Print values with 8 decimal places
    cat(sprintf("%14.8f %16.8f %16.8f %16.8f\n",
                object$beta_hat[i], object$std_beta_hat[i], lower_ci, upper_ci))
  }

  cat("-------------------------------------------------------------\n")
  cat("Stochastic parameters\n")
  cat("-------------------------------------------------------------\n")
  cat(sprintf(" White Noise Variance  : %14.8f\n", object$gamma_hat[1]))
  cat(sprintf(" Flicker Noise Variance: %14.8f\n", object$gamma_hat[2]))

  cat("-------------------------------------------------------------\n")
}




#' Plot a \code{fit_gnss_ts_ngl} object
#' @param x A \code{fit_gnss_ts_ngl} object.
#' @param ... Additional graphical parameters.
#' @export
#' @return No return value. Plot a \code{fit_gnss_ts_ngl} object.
plot.fit_gnss_ts_ngl = function(x, ...){
  #
  #
  # library(gmwmx2)
  # station_data = download_station_ngl("0AMB")
  # x = gmwmx2(station_data, n_seasonal = 2, component = "N")

  # Save the current graphical parameters
  old_par <- par(no.readonly = TRUE)

  # set parameters for layout
  mat_layout = matrix(c(1,2, 3),ncol=1, nrow=3)
  layout(mat_layout, heights = c(.1,1, 1))
  par(mar=c(0,0,0,0))
  plot.new()
  legend("center", horiz = T,
         legend = c("NA", "Equipment/Software change", "Earthquake", "Estimated fit"),
         col = c("grey60", "blue", "darkorange", "red"),
         pch = c(15, NA, NA, NA),
         pt.cex=c(2, NA, NA, NA),
         # x.intersp = 0.8,
         text.width = c(.1,.3,.2,.1)  ,
         lty=c(NA, 1,1,1), bty="n")
  par(mar=c(4, 4.1, 2, 2.1))


  component = x$component
  if(component == "E"){
    axis_name = "Easting (m)"
  }else if(component == "N"){
    axis_name = "Northing (m)"
  }else if(component == "V"){
    axis_name = "Vertical (m)"
  }

  # plot data
  plot(x$design_matrix_X[,2], x$y, type="l", las=1,ylab=axis_name, xlab="MJD")

  # add estimated fit
  lines(x=x$design_matrix_X[,2], y= x$design_matrix_X%*% x$beta_hat, col="red")

  # compute NA over the time series
  all_mjd = seq(head(x$df_position$modified_julian_day,1),
                tail(x$df_position$modified_julian_day,1))
  missing_mjd = all_mjd[which(!all_mjd %in% x$df_position$modified_julian_day)]

  # missing data
  for(i in seq_along(missing_mjd)){
    abline(v = missing_mjd[i], col="grey60")
  }

  # add equipment change
  for(i in seq((dim(x$df_equipment_software_changes)[1]))){
    abline(v = x$df_equipment_software_changes$modified_julian_date, col="blue")
  }

  # add earthquake
  for(i in seq((dim(x$df_earthquakes)[1]))){
    abline(v = x$df_earthquakes$modified_julian_date, col="darkorange")
  }
  box()

  plot(x$empirical_wvar)
  # add theoretical implied wv
  lines(x=x$empirical_wvar$scales, y=x$theoretical_wvar, col="red")
  box()

  # Restore the original graphical parameters
  par(old_par)

  # Reset to a single plot layout
  layout(1)
}
