create_X_matrix <- function(all_mjd_index,
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
  if (!is.null(vec_earthquakes_relaxation_time)) {
    if (length(vec_earthquakes_index_mjd) != length(vec_earthquakes_relaxation_time)) {
      stop("Vector of relaxation time not of same length as vector of MJD index of earthquakes")
    }
  }


  # set to default 365.25 days if no earthquakes relaxation time provided
  if (is.null(vec_earthquakes_relaxation_time)) {
    vec_earthquakes_relaxation_time <- rep(365.25, length(vec_earthquakes_index_mjd))
  }


  # create empty matrix number of columns bias+trend+2*nbr sinusoidal + nbr jumps + number earthquakes
  X <- matrix(0, nrow = length(all_mjd_index), ncol = 2 + 2 * n_seasonal +
    length(jumps) + length(vec_earthquakes_index_mjd))

  # add bias, intercept
  X[, 1] <- 1

  # add component for trend and scale with respect to middle of time axis
  reference_time <- 0.5 * (all_mjd_index[1] + tail(all_mjd_index, 1))
  X[, 2] <- all_mjd_index - reference_time

  # add seasonal
  if (n_seasonal > 0) {
    for (i in 1:n_seasonal) {
      omega_i <- (i / 365.25) * 2 * pi
      X[, 2 + (i - 1) * 2 + 1] <- sin((all_mjd_index) * omega_i)
      X[, 2 + (i - 1) * 2 + 2] <- cos((all_mjd_index) * omega_i)
    }
  }

  # add offsets
  if (!is.null(jumps)) {
    for (i in 1:length(jumps)) {
      it <- min(which(all_mjd_index > jumps[i] - 1e-06))
      X[, 2 + 2 * n_seasonal + i] <- c(rep(0, it - 1), rep(
        1,
        length(all_mjd_index) - it + 1
      ))
    }
  }

  # exponential decay function for post seismic relaxation
  if (!is.null(vec_earthquakes_index_mjd)) {
    for (i in seq_along(vec_earthquakes_index_mjd)) {
      tau_i <- vec_earthquakes_relaxation_time[i]
      earthquake_mjd_i <- vec_earthquakes_index_mjd[i]
      # create vector
      decay_values <- ifelse(all_mjd_index >= earthquake_mjd_i,
        1 - exp(-(all_mjd_index - earthquake_mjd_i) / tau_i),
        0
      )

      # create column in matrix
      X[, 2 + 2 * n_seasonal + length(jumps) + i] <- decay_values
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
  rownames(X) <- all_mjd_index
  return(X)
}



# transformation function for spectral index of power law noise
trans_kappa_pl <- function(x) {
  exp(x) - 1
}


# inverse of transformation function for spectral index of power law noise
inv_trans_kappa_pl <- function(x) log(x + 1)


# define optimization function for stochastic model White noise + Flicker model or White noise + stationary powerlaw
# In order to perform the optimization efficiently,
# the objective function construct a fast approximation of the theoretical wavelet variance of the vector of missing and observed estimated residuals.
objective_function_w_missing <- function(theta, wv_obj, n, quantities_D, approx_type, vec_autocov_omega, pstar_hat, no_missing = T, omega = NULL, stochastic_model) {
  if (stochastic_model == "wn + fl") {
    theta_t <- vector(mode = "numeric", length = 2)
    theta_t[1] <- exp(theta[1]) # sigma2 wn
    theta_t[2] <- exp(theta[2]) # sigma2 fl

    vec_mean_autocov <- vec_mean_autocov_powerlaw(kappa = -1, n) * theta_t[2]
    vec_mean_autocov[1] <- vec_mean_autocov[1] + theta_t[1]
  } else if (stochastic_model == "wn + pl") {
    theta_t <- vector(mode = "numeric", length = 3)
    theta_t[1] <- exp(theta[1]) # sigma2 wn
    theta_t[2] <- trans_kappa_pl(theta[2]) # sigma2 wn
    theta_t[3] <- exp(theta[3]) # sigma2 pl

    vec_mean_autocov <- powerlaw_autocovariance(kappa = theta_t[2], sigma2 = theta_t[3], n = n)
    vec_mean_autocov[1] <- vec_mean_autocov[1] + theta_t[1]
  }




  # approx with linear interpolation on errors
  vec_mean_autocov_eps_hat <- compute_all_mean_diag_fast_w_linear_interp_only_required_cpp(
    mat_D_q_term_1 = quantities_D$mat_D_q_term_1,
    mat_D_q_term_2 = quantities_D$mat_D_q_term_2,
    sum_on_sub_diag_of_D = quantities_D$sum_on_sub_diag_of_D,
    vec_autocov = vec_mean_autocov, approx_type = approx_type
  )
  if (!no_missing) {
    vec_mean_per_diag_w_missing <- vec_mean_autocov_eps_hat * (vec_autocov_omega + pstar_hat^2)
    theo_wv <- autocovariance_to_wv(vec_mean_per_diag_w_missing, tau = wv_obj$scales)
  } else {
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



#' Estimate a trajectory model for a \code{gnss_ts_ngl} object considering a white noise plus colored noise as the stochastic model for the residuals and model missingness with a Markov process using the GMWMX estimator.
#' @param x A \code{gnss_ts_ngl} object.
#' @param n_seasonal An \code{integer} specifying the number of seasonal signals in the time series. "1" specify only one annual periodic signal and "2"specify an annual and a semiannual periodic signal.
#' @param vec_earthquakes_relaxation_time A \code{vector} specifying the relaxation time for each earthquakes indicated for the time series.
#' @param component A \code{string} with value either "N", "E" or "V" that specify which component to estimate (Northing, Easting or Vertical).
#' @param toeplitz_approx_var_cov_wv A \code{boolean} that specify if the variance of the wavelet variance should be computed based on a toeplitz approximation of the variance covariance matrix of the residuals.
#' @param stochastic_model A \code{string} that specify the stochastic model considered for the residuals. Either "wn + fl" for white noise and flicker/pink noise or "wn + pl" for white noise and stationary power-law noise.
#' @importFrom wv wvar
#' @importFrom dplyr between
#' @importFrom Matrix solve
#' @importFrom stats .lm.fit optim toeplitz
#' @examples
#' x <- download_station_ngl("CHML")
#' fit <- gmwmx2(x, n_seasonal = 2, component = "N")
#' @export
gmwmx2 <- function(x, n_seasonal = 2, vec_earthquakes_relaxation_time = NULL, component = "N", toeplitz_approx_var_cov_wv = TRUE, stochastic_model = "wn + fl") {
  # x = download_station_ngl("STPS")
  # plot(x)
  # # vec_earthquakes_relaxation_time <- NULL
  # component <- "N"
  # n_seasonal <- 2
  # toeplitz_approx_var_cov_wv=TRUE
  # stochastic_model = "wn + pl"
  # all_station = download_all_stations_ngl()
  # all_station[200]
  # # x = download_station_ngl("0ABN") # nice to show diff between pl and flicker
  # x=download_station_ngl("0KIS")
  # plot(x)
  # fit1 = gmwmx2(x = x, stochastic_model = "wn + fl", component = "E")
  # plot(fit1)
  # fit2 = gmwmx2(x = x, stochastic_model = "wn + pl", component = "E")
  # plot(fit2)

  # Check class
  if (!inherits(x, "gnss_ts_ngl")) {
    stop("Argument `x` should be a `gnss_ts_ngl` object")
  }

  # check stochastic model
  if (!stochastic_model %in% c("wn + pl", "wn + fl")) {
    stop("Argument `stochastic_model` should be either `wn + fl` or `wn + pl`")
  }

  # check that component is either N, E or V
  if (!component %in% c("N", "E", "V")) {
    stop("Argument `component` should take either value `N` or `E` or `V`")
  }

  # check that n_seasonal is either 1 or 2
  if (!n_seasonal %in% c(1, 2)) {
    stop("Argument `n_seasonal` should take either value `1` or `2`")
  }


  # Record start time
  start_time_gmwmx <- Sys.time()

  # create full index
  all_mjd_index <- seq(head(x$df_position$modified_julian_day, 1), tail(x$df_position$modified_julian_day, 1), by = 1)

  # create all jumps by combining jumps due to equipment change and jumps due to earthquakes
  jumps <- c(
    x$df_equipment_software_changes$modified_julian_date,
    x$df_earthquakes$modified_julian_date
  )

  # if multiple jumps  to prevent not invertible matrix
  jumps <- unique(jumps)
  vec_earthquakes_index_mjd <- c(x$df_earthquakes$modified_julian_date)

  # if multiple earthquakes to prevent not invertible matrix
  vec_earthquakes_index_mjd <- unique(vec_earthquakes_index_mjd)

  if (length(jumps) == 0) {
    jumps <- NULL
  }

  # ensure that no jumps or earthquake mjd are specified after the last date of recorded signal
  last_mjd_signal = tail(all_mjd_index,1)
  id_jumps_to_remove = which(jumps > last_mjd_signal)
  id_earthquake_index_to_remove = which(vec_earthquakes_index_mjd > last_mjd_signal)

  # Remove the identified indices if they exist
  if (length(id_jumps_to_remove) > 0) {
    jumps <- jumps[-id_jumps_to_remove]
  }

  if (length(id_earthquake_index_to_remove) > 0) {
    vec_earthquakes_index_mjd <- vec_earthquakes_index_mjd[-id_earthquake_index_to_remove]
  }

  # create design matrix
  X <- create_X_matrix(
    all_mjd_index = all_mjd_index,
    jumps = jumps,
    n_seasonal = n_seasonal,
    vec_earthquakes_index_mjd = vec_earthquakes_index_mjd,
    vec_earthquakes_relaxation_time = vec_earthquakes_relaxation_time
  )

  # obtain X_sub
  id_X_sub <- which(rownames(X) %in% x$df_position$modified_julian_day)
  X_sub <- X[id_X_sub, ]

  # Extract Y given specified component
  if (component == "N") {
    y <- x$df_position$northings_fractional_portion
  } else if (component == "E") {
    y <- x$df_position$eastings_fractional_portion
  } else if (component == "V") {
    y <- x$df_position$vertical_fractional_portion
  }

  # obtain beta hat
  beta_hat <- .lm.fit(y = y, x = X_sub)$coefficients

  # create vector of name of parameters
  if (n_seasonal == 1) {
    names_beta_hat <- c(
      "Intercept", "Trend", "Sin (Annual)", "Cos (Annual)",
      if (length(jumps) > 0) paste0("Jump: ", jumps),
      if (length(vec_earthquakes_index_mjd) > 0) paste0("Earthquake: ", vec_earthquakes_index_mjd)
    )
  } else if (n_seasonal == 2) {
    names_beta_hat <- c(
      "Intercept", "Trend", "Sin (Annual)", "Cos (Annual)",
      "Sin (Semi-Annual)", "Cos (Semi-Annual)",
      if (length(jumps) > 0) paste0("Jump: MJD ", jumps),
      if (length(vec_earthquakes_index_mjd) > 0) paste0("Earthquake: MJD ", vec_earthquakes_index_mjd)
    )
  }

  # assign names beta hat
  names(beta_hat) <- names_beta_hat

  # # plot signal and estimated model
  # plot(x=rownames(X_sub), y=y, type="l")
  # lines(x=rownames(X_sub), y= X_sub%*% beta_hat, col="red")

  # # obtain observed residuals
  eps_hat_sub <- y - X_sub %*% beta_hat

  # create vector that contains either estimated residual or 0 if data point is not observed
  eps_hat_filled <- vector(mode = "numeric", length = length(all_mjd_index))

  # fill in observed residuals when we observe data
  eps_hat_filled[id_X_sub] <- eps_hat_sub

  # compute empirical wv variance on the "filled" vector of residuals
  wv_emp_eps_hat_filled <- wv::wvar(eps_hat_filled)

  # plot(wv_emp_eps_hat_filled)
  # define vector of observed index
  vec_omega <- as.numeric(all_mjd_index %in% x$df_position$modified_julian_day)

  # estimate parameter of Markov process missingness process using the MLE
  p_hat <- estimate_p1_p2_mle_cpp(vec_omega)

  # define pstar hat (expecation of missingness process)
  pstar_hat <- p_hat[2] / (p_hat[1] + p_hat[2])

  # get vec autocovariance theo omega
  vec_autocov_omega <- create_vec_theo_autocov_omega_cpp(p1 = p_hat[1], p2 = p_hat[2], length(all_mjd_index))

  # compute (X^TX)^{-1} using QR decomposition for numerical stability
  XtX <- t(X) %*% X
  qr_decomp <- qr(X)
  R <- qr.R(qr_decomp)
  R_inv <- Matrix::solve(R)
  inv_XtX <- R_inv %*% t(R_inv)


  # inv_XtX <- Matrix::solve(XtX)

  # compute hat matrix
  H <- X %*% inv_XtX %*% t(X)
  D <- diag(length(all_mjd_index)) - H

  # pre-compute quantities on D=(I-H), with H = X(X^TX)^{-1}X^T to later use in optimization function
  quantities_D <- pre_compute_quantities_on_D_only_required_smarter_cpp(D, approx_type = "3")

  # define missingness case
  if (all(vec_omega == 1)) {
    no_missing <- TRUE
  } else {
    no_missing <- FALSE
  }

  # define gamma init based o nspecified stochastic model
  if (stochastic_model == "wn + fl") {
    gamma_init_not_transformed <- find_initial_values_wn_fl(wv_emp_eps_hat_filled)
    # log transform variance of WN and variance of flicker noise
    gamma_init <- log(gamma_init_not_transformed)
  } else if (stochastic_model == "wn + pl") {
    # get initial values
    gamma_init_not_transformed <- find_initial_values_wn_pl(signal = eps_hat_sub, wv_emp = wv_emp_eps_hat_filled)
    # log transform variance of WN, transform kappa and log transform variance of powerlaw noise
    gamma_init <- c(
      log(gamma_init_not_transformed[1]),
      inv_trans_kappa_pl(gamma_init_not_transformed[2]),
      log(gamma_init_not_transformed[3])
    )
  }
  # fit gmwmx on empirical wv
  res_gmwmx <- optim(
    par = gamma_init,
    fn = objective_function_w_missing,
    wv_obj = wv_emp_eps_hat_filled,
    n = length(vec_omega),
    quantities_D = quantities_D,
    approx_type = "3",
    vec_autocov_omega = vec_autocov_omega,
    pstar_hat = pstar_hat,
    no_missing = no_missing,
    stochastic_model = stochastic_model
  )

  # extract estimated stochastic parameters
  if (stochastic_model == "wn + fl") {
    gamma_hat_1 <- exp(res_gmwmx$par)
  } else if (stochastic_model == "wn + pl") {
    gamma_hat_1 <- c(
      exp(res_gmwmx$par[1]),
      trans_kappa_pl(res_gmwmx$par[2]),
      exp(res_gmwmx$par[3])
    )
  }

  if (toeplitz_approx_var_cov_wv) {
    # define max J
    max_J <- wv_emp_eps_hat_filled$J

    # minimum n required for computing var of wv
    min_n_to_get_var_wv <- length(vec_omega) + sum_of_powers_of_2(1, max_J - 1) + 2

    # get autocov of stochastic process (true residuals)
    if (stochastic_model == "wn + fl") {
      vec_mean_autocov <- vec_mean_autocov_powerlaw(-1, length(vec_omega)) * gamma_hat_1[2]
      vec_mean_autocov[1] <- vec_mean_autocov[1] + gamma_hat_1[1]
    } else if (stochastic_model == "wn + pl") {
      vec_mean_autocov <- powerlaw_autocovariance(kappa = gamma_hat_1[2], sigma2 = gamma_hat_1[3], n = length(vec_omega))
      vec_mean_autocov[1] <- vec_mean_autocov[1] + gamma_hat_1[1]
    }


    # get mean autocov of (I-H)Sigma
    vec_mean_autocov_eps_hat <- compute_all_mean_diag_fast_w_linear_interp_only_required_cpp(
      mat_D_q_term_1 = quantities_D$mat_D_q_term_1,
      mat_D_q_term_2 = quantities_D$mat_D_q_term_2,
      sum_on_sub_diag_of_D = quantities_D$sum_on_sub_diag_of_D,
      vec_autocov = vec_mean_autocov,
      approx_type = "3"
    )

    # get autocov of process times markov process
    if (!no_missing) {
      autocov_wn_fl_times_omega <- vec_mean_autocov_eps_hat * (vec_autocov_omega + pstar_hat^2)
    } else {
      autocov_wn_fl_times_omega <- vec_mean_autocov_eps_hat
    }

    # fill in zero after this entry to compute autocovariance of wavelet coefficients of all scales
    # autocov of wavelet coef will be wrong for some entries but the one used will be correct
    autocov_wn_fl_times_omega_w_zeroes <- vector(mode = "numeric", length = min_n_to_get_var_wv)
    autocov_wn_fl_times_omega_w_zeroes[1:(length(vec_omega))] <- autocov_wn_fl_times_omega

    # get var cov missingness process
    var_cov_omega <- fast_toeplitz_matrix_from_vector_cpp(as.vector(vec_autocov_omega))

    # get variance of wv based on this process
    Sigma_wv <- get_theo_cov_matrix_wvar_cpp(n = length(vec_omega), autocov_vec_X = autocov_wn_fl_times_omega_w_zeroes)
    inv_var_cov_nu_hat <- Matrix::solve(Sigma_wv)

    # Build required matrices to compute later variance covariance of beta hat
    if (stochastic_model == "wn + fl") {
      var_cov_mat_wn <- gamma_hat_1[1] * diag(length(vec_omega))
      var_cov_mat_flicker <- var_cov_powerlaw_cpp(sigma2 = gamma_hat_1[2], kappa = -1, n = length(vec_omega))
      var_cov_mat_epsilon <- var_cov_mat_wn + var_cov_mat_flicker
    } else if (stochastic_model == "wn + pl") {
      vec_autocov_stationary_powerlaw <- powerlaw_autocovariance(kappa = gamma_hat_1[2], sigma2 = gamma_hat_1[3], n = length(vec_omega))
      vec_autocov_stationary_powerlaw[1] <- gamma_hat_1[1] + vec_autocov_stationary_powerlaw[1]
      var_cov_mat_epsilon <- fast_toeplitz_matrix_from_vector_cpp(as.vector(vec_autocov_stationary_powerlaw))
    }
  } else {
    # Compute the variance covariance of the Wavelet variance when the process is not stationary and therefore using the whole variance covariance matrix of the process on which it is computed
    # construct sigma matrix of white noise + ficker

    # Obtain variance covariance of epsilon
    if (stochastic_model == "wn + fl") {
      var_cov_mat_wn <- gamma_hat_1[1] * diag(length(vec_omega))
      var_cov_mat_flicker <- var_cov_powerlaw_cpp(sigma2 = gamma_hat_1[2], kappa = -1, n = length(vec_omega))
      var_cov_mat_epsilon <- var_cov_mat_wn + var_cov_mat_flicker
    } else if (stochastic_model == "wn + pl") {
      vec_autocov_stationary_powerlaw <- powerlaw_autocovariance(kappa = gamma_hat_1[2], sigma2 = gamma_hat_1[3], n = length(vec_omega))
      vec_autocov_stationary_powerlaw[1] <- gamma_hat_1[1] + vec_autocov_stationary_powerlaw[1]
      var_cov_mat_epsilon <- fast_toeplitz_matrix_from_vector_cpp(as.vector(vec_autocov_stationary_powerlaw))
    }

    # get var cov missingness process
    var_cov_omega <- fast_toeplitz_matrix_from_vector_cpp(as.vector(vec_autocov_omega))

    # define variance covariance of residuals with missing
    var_cov_eps_hat_w_missing <- var_cov_mat_epsilon * (var_cov_omega + pstar_hat^2)

    # compute variance of wavelet variance computed on epsilon (negleting the fact that we compute the wv on the "estimated" residuals)
    var_cov_nu_hat <- compute_cov_wv_cpp_approx_faster(Sigma_X = var_cov_eps_hat_w_missing)

    # get inverse
    inv_var_cov_nu_hat <- Matrix::solve(var_cov_nu_hat)
  }


  # re estimate gamma with optimal weighting matrix
  res_gmwmx_2 <- optim(
    par = res_gmwmx$par,
    fn = objective_function_w_missing,
    wv_obj = wv_emp_eps_hat_filled,
    n = length(vec_omega),
    quantities_D = quantities_D,
    approx_type = "3",
    vec_autocov_omega = vec_autocov_omega,
    omega = inv_var_cov_nu_hat,
    pstar_hat = pstar_hat,
    no_missing = no_missing,
    stochastic_model = stochastic_model
  )

  if (stochastic_model == "wn + fl") {
    gamma_hat_2 <- exp(res_gmwmx_2$par)
  } else if (stochastic_model == "wn + pl") {
    gamma_hat_2 <- c(
      exp(res_gmwmx_2$par[1]),
      trans_kappa_pl(res_gmwmx_2$par[2]),
      exp(res_gmwmx_2$par[3])
    )
  }

  # get theo wv with approximation obtained by last fit
  if (stochastic_model == "wn + fl") {
    vec_mean_autocov <- vec_mean_autocov_powerlaw(kappa = -1, length(vec_omega)) * gamma_hat_2[2]
    vec_mean_autocov[1] <- vec_mean_autocov[1] + gamma_hat_2[1]
  } else if (stochastic_model == "wn + pl") {
    vec_mean_autocov <- powerlaw_autocovariance(kappa = gamma_hat_2[2], sigma2 = gamma_hat_2[3], n = length(vec_omega))
    vec_mean_autocov[1] <- vec_mean_autocov[1] + gamma_hat_2[1]
    vec_mean_autocov <- as.vector(vec_mean_autocov)
  }


  # compute vector of average of autocovariance of Sigma matrix of epsilon hat using approximation
  vec_mean_autocov_eps_hat <- compute_all_mean_diag_fast_w_linear_interp_only_required_cpp(
    mat_D_q_term_1 = quantities_D$mat_D_q_term_1,
    mat_D_q_term_2 = quantities_D$mat_D_q_term_2,
    sum_on_sub_diag_of_D = quantities_D$sum_on_sub_diag_of_D,
    vec_autocov = vec_mean_autocov, approx_type = "3"
  )

  # compute implied theoretical wavelet variance at estimated parameters
  if (!no_missing) {
    vec_mean_per_diag_w_missing <- vec_mean_autocov_eps_hat * (vec_autocov_omega + pstar_hat^2)
    theo_wv <- autocovariance_to_wv(vec_mean_per_diag_w_missing, tau = wv_emp_eps_hat_filled$scales)
  } else {
    theo_wv <- autocovariance_to_wv(vec_mean_autocov_eps_hat, tau = wv_emp_eps_hat_filled$scales)
  }



  # Compute variance covariance of beta hat
  if (no_missing) {
    var_cov_beta_hat <- inv_XtX %*% t(X) %*% var_cov_mat_epsilon %*% X %*% inv_XtX
  } else {
    var_cov_beta_hat <- pstar_hat^(-2) * inv_XtX %*% t(X) %*% ((var_cov_omega + pstar_hat^2) * var_cov_mat_epsilon) %*% X %*% inv_XtX
  }

  std_beta_hat_gmwmx <- sqrt(diag(var_cov_beta_hat))

  # Record end time
  end_time_gmwmx <- Sys.time()

  # get time gmwmx 1
  time_gmwmx <- as.numeric(difftime(end_time_gmwmx, start_time_gmwmx, units = "secs"))

  ret <- list(
    "beta_hat" = beta_hat,
    "std_beta_hat" = std_beta_hat_gmwmx,
    "gamma_hat" = gamma_hat_2,
    "vartheta_hat" = p_hat,
    "component" = component,
    "design_matrix_X" = X_sub,
    "y" = y,
    "empirical_wvar" = wv_emp_eps_hat_filled,
    "theoretical_wvar" = theo_wv,
    "df_position" = x$df_position,
    "df_earthquakes" = x$df_earthquakes,
    "df_equipment_software_changes" = x$df_equipment_software_changes,
    "p_hat" = p_hat,
    "p_star_hat" = pstar_hat,
    "stochastic_model" = stochastic_model,
    "running_time" = time_gmwmx
  )

  class(ret) <- "fit_gnss_ts_ngl"

  return(ret)
}



#' Extract estimated parameters from a \code{fit_gnss_ts_ngl}
#' @param object A \code{fit_gnss_ts_ngl} object.
#' @param scale_parameters A \code{boolean} indicating whether or not to scale estimated parameters so that the returned estimated trend is provided in m/year instead of m/day. Default is FALSE.
#' @param ... Additional parameters.
#' @examples
#' x <- download_station_ngl("P820")
#' fit1 <- gmwmx2(x, n_seasonal = 2, component = "N", stochastic_model = "wn + pl")
#' summary(fit1)
#' summary(fit1, scale_parameters = TRUE)
#' fit2 <- gmwmx2(x, n_seasonal = 2, component = "N", stochastic_model = "wn + fl")
#' summary(fit2)
#'
#' @export
summary.fit_gnss_ts_ngl <- function(object, scale_parameters = FALSE, ...) {
  # Print header
  cat("Summary of Estimated Model\n")
  cat("-------------------------------------------------------------\n")
  cat("Functional parameters\n")
  cat("-------------------------------------------------------------\n")

  cat("Parameter                  Estimate  Std_Deviation  95% CI Lower  95% CI Upper\n")
  cat("-------------------------------------------------------------\n")


  # define normal confidence interval for parameters of the trajectory model (betas), optionnal scaling
  for (i in seq_along(object$beta_hat)) {
    if (scale_parameters) {
      lower_ci <- object$beta_hat[i] * 365.25 - qnorm(.975) * object$std_beta_hat[i] * 365.25
      upper_ci <- object$beta_hat[i] * 365.25 + qnorm(.975) * object$std_beta_hat[i] * 365.25
    } else {
      lower_ci <- object$beta_hat[i] - qnorm(.975) * object$std_beta_hat[i]
      upper_ci <- object$beta_hat[i] + qnorm(.975) * object$std_beta_hat[i]
    }


    # Print values with 8 decimal places, left-align the names, and align other columns
    if (scale_parameters) {
      cat(sprintf(
        "%-22s %12.8f %12.8f %12.8f %12.8f\n",
        names(object$beta_hat)[i],
        object$beta_hat[i] * 365.25, object$std_beta_hat[i] * 365.25,
        lower_ci, upper_ci
      ))
    } else {
      cat(sprintf(
        "%-22s %12.8f %12.8f %12.8f %12.8f\n",
        names(object$beta_hat)[i],
        object$beta_hat[i], object$std_beta_hat[i],
        lower_ci, upper_ci
      ))
    }
  }

  cat("-------------------------------------------------------------\n")
  cat("Stochastic parameters\n")
  cat("-------------------------------------------------------------\n")

  if (object$stochastic_model == "wn + fl") {
    cat(sprintf(" White Noise Variance  : %14.8f\n", object$gamma_hat[1]))
    cat(sprintf(" Flicker Noise Variance: %14.8f\n", object$gamma_hat[2]))
  } else if (object$stochastic_model == "wn + pl") {
    cat(sprintf(" White Noise Variance  : %14.8f\n", object$gamma_hat[1]))
    cat(sprintf(" Stationary powerlaw Spectral index: %14.8f\n", object$gamma_hat[2]))
    cat(sprintf(" Stationary powerlaw Variance: %14.8f\n", object$gamma_hat[3]))
  }



  cat("-------------------------------------------------------------\n")
  cat("Missingness parameters\n")
  cat("-------------------------------------------------------------\n")
  cat(sprintf(" P(Z_{i+1} = 0 | Z_{i} = 1): %.8f\n", object$p_hat[1]))
  cat(sprintf(" P(Z_{i+1} = 1 | Z_{i} = 0): %.8f\n", object$p_hat[2]))
  cat(sprintf(" \\hat{E[Z]}: %.8f\n", object$p_star_hat))
  cat("-------------------------------------------------------------\n")
  cat(paste0("Running time: ", round(object$running_time, 2), " seconds", "\n"))
  cat("-------------------------------------------------------------\n")
}




#' Plot a \code{fit_gnss_ts_ngl} object
#' @importFrom graphics axis grid mtext polygon
#' @importFrom stats qnorm
#' @param x A \code{fit_gnss_ts_ngl} object.
#' @param ... Additional graphical parameters.
#' @export
#' @examples
#' x <- download_station_ngl("0AMB")
#' fit_N <- gmwmx2(x, n_seasonal = 2, component = "N")
#' plot(fit_N)
#' fit_E <- gmwmx2(x, n_seasonal = 2, component = "E")
#' plot(fit_E)
#' @return No return value. Plot a \code{fit_gnss_ts_ngl} object.
plot.fit_gnss_ts_ngl <- function(x, ...) {
  #
  #
  # library(gmwmx2)
  # station_data = download_station_ngl("0AMB")
  # x = gmwmx2(station_data, n_seasonal = 2, component = "E")

  # Save the current graphical parameters
  old_par <- par(no.readonly = TRUE)

  # set parameters for layout
  mat_layout <- matrix(c(1, 2, 3, 4), ncol = 1, nrow = 4)
  layout(mat_layout, heights = c(.1, .1, 1, 1))
  par(mar = c(0, 0, 0, 0))
  plot.new()
  if (x$component == "N") {
    axis_name <- "Northing (m)"
  } else if (x$component == "E") {
    axis_name <- "Easting (m)"
  } else if (x$component == "V") {
    axis_name <- "Vertical (m)"
  }


  if (x$stochastic_model == "wn + pl") {
    name_stoch <- "White noise and Powerlaw"
  } else if (x$stochastic_model == "wn + fl") {
    name_stoch <- "White noise and Flicker"
  }
  text_station <- c(paste0(
    "Station ", unique(x$df_position$station_name),
    ": ",
    axis_name, " | ",
    name_stoch
  ))
  legend("center",
    horiz = T,
    legend = text_station,
    bty = "n", cex = 1.3
  )
  plot.new()
  legend("center",
    horiz = T,
    legend = c("NA", "Equipment/Software change", "Earthquake", "Estimated fit"),
    col = c("grey60", "blue", "darkorange", "red"),
    pch = c(15, NA, NA, NA),
    pt.cex = c(2, NA, NA, NA),
    # x.intersp = 0.8,
    text.width = c(.1, .3, .2, .1),
    lty = c(NA, 1, 1, 1), bty = "n"
  )
  par(mar = c(3.5, 4.5, 2, 2.1))


  component <- x$component
  if (component == "E") {
    axis_name <- "Easting (m)"
  } else if (component == "N") {
    axis_name <- "Northing (m)"
  } else if (component == "V") {
    axis_name <- "Vertical (m)"
  }

  # plot data
  plot(x = rownames(x$design_matrix_X), y = x$y, type = "l", las = 1, ylab = "", xlab = "")
  mtext(side = 1, "Modified Julian Date", line = 3)
  mtext(side = 2, axis_name, line = 3.2)

  grid(col = "grey90", lty = 2)
  # add signal
  lines(x = rownames(x$design_matrix_X), y = x$y)

  # add estimated fit
  lines(x = rownames(x$design_matrix_X), y = x$design_matrix_X %*% x$beta_hat, col = "red")

  # compute NA over the time series
  all_mjd <- seq(
    head(x$df_position$modified_julian_day, 1),
    tail(x$df_position$modified_julian_day, 1)
  )
  missing_mjd <- all_mjd[which(!all_mjd %in% x$df_position$modified_julian_day)]

  # missing data
  for (i in seq_along(missing_mjd)) {
    abline(v = missing_mjd[i], col = "grey60")
  }

  # add equipment change
  for (i in seq((dim(x$df_equipment_software_changes)[1]))) {
    abline(v = x$df_equipment_software_changes$modified_julian_date, col = "blue")
  }

  # add earthquake
  for (i in seq((dim(x$df_earthquakes)[1]))) {
    abline(v = x$df_earthquakes$modified_julian_date, col = "darkorange")
  }
  box()


  # plot empirical WV and theorical implied wv
  yl <- c(min(x$empirical_wvar$ci_low), max(x$empirical_wvar$ci_high))

  plot(NA,
    ylim = yl,
    xlim = c(min((x$empirical_wvar$scales)), max(x$empirical_wvar$scales)),
    main = "",
    ylab = "",
    xlab = "",
    log = "xy",
    xaxt = "n",
    yaxt = "n"
  )
  mtext(side = 1, text = "Scales", line = 2.3)
  mtext(side = 2, text = "Wavelet Variance", line = 3.2)

  # Create a vector of labels in the desired format
  exponents <- log2(x$empirical_wvar$scales)
  labels <- sapply(exponents, function(x) as.expression(bquote(2^.(x))))

  axis(1, at = x$empirical_wvar$scales, labels = labels)

  polygon(
    x = c(x$empirical_wvar$scales, rev(x$empirical_wvar$scales)),
    y = c(x$empirical_wvar$ci_low, rev(x$empirical_wvar$ci_high)),
    col = "#ccf1f8",
    border = NA
  )

  ylab <- floor(log10(yl[1])):ceiling(log10(yl[2]))
  axis(2, at = 10^ylab, labels = parse(text = sprintf("10^%.0f", ylab)), las = 1)

  # add grid
  abline(h = 10^ylab, col = "grey90", lt = 2)
  for (i in seq_along(exponents)) {
    abline(v = 2^exponents[i], col = "grey90", lt = 2)
  }

  lines(x$empirical_wvar$scales, x$empirical_wvar$variance, type = "b", col = "black", pch = 16)
  lines(x = x$empirical_wvar$scales, y = x$theoretical_wvar, type = "b", col = "red", pch = 21, cex = 1.4)
  # add empirical WV

  legend(
    "bottomleft",
    legend = c("Empirical WV", "Estimated Theoretical WV"),
    col = c("black", "red"),
    pch = c(16, 21),
    pt.cex = c(1, 1.4),
    lty = c(1),
    horiz = FALSE,
    bty = "n", bg = "transparent"
  )

  box()

  # Restore the original graphical parameters
  par(old_par)

  # Reset to a single plot layout
  layout(1)
}
