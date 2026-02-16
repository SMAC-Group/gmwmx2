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

  # prevent agains optimization problem diverging due to numerical issues, if objective is not finite, return a large number
  # if(!is.finite(objective)) return(1e30)

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

  difference <- nu_hat - theo_wv
  objective <- as.numeric(t(difference) %*% omega %*% difference)

  # prevent agains optimization problem diverging due to numerical issues, if objective is not finite, return a large number
  # if(!is.finite(objective)) return(1e30)


  return(objective)

}
















#' Variance-covariance matrix implied by a model
#'
#' Constructs the variance-covariance matrix for a model using its
#' `get_variance_covariance_matrix_signal` method. Supports both single
#' `time_series_model` objects and composite `sum_model` objects.
#'
#' @param model A `time_series_model` or `sum_model`.
#' @param n Length of the signal.
#' @param theta Optional parameter vector (already in domain space of the parameters).
#' @param prep Optional output from `prepare_optim_layout(model)` (unused for now).
#' @return Variance-covariance matrix of dimension `n x n`.
#' @keywords internal
get_variance_covariance_matrix_model <- function(model, n, theta = NULL, prep = NULL, ...) {
  UseMethod("get_variance_covariance_matrix_model")
}

#' @keywords internal
get_variance_covariance_matrix_model.time_series_model <- function(model, n, theta = NULL, prep = NULL, ...) {

  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n <= 0L) {
    stop("`n` must be a positive integer.", call. = FALSE)
  }
  # if theta not provided, then use model$parameters, otherwise use theta, but set names to theta according to the order of parameters in model$parameters
  if(is.null(theta)){
    pars <- model$parameters
    if (is.null(pars) || any(is.null(pars))) {
      stop("Model parameters must be set (not NULL).", call. = FALSE)
    }
    do.call(model$get_variance_covariance_matrix_signal, c(as.list(pars), list(n = n)))
    # else evaluate get_variance_covariance_matrix_signal with theta, but set names to theta according to the order of parameters in model$parameters
  }else{
    # set names to theta according to the order of parameters in model$parameters
    names(theta) = names(model$parameters)
    do.call(model$get_variance_covariance_matrix_signal, c(as.list(theta), list(n)))
  }
}

#' @keywords internal
get_variance_covariance_matrix_model.sum_model <- function(model, n, theta = NULL, prep = NULL, ...) {


  #------------------------
  # model = ar1() +rw()
  # n = 10
  # theta = c(.5,1,1)
  # prep = prepare_optim_layout(ar1() +rw())
  #------------------------

  # check on n
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n <= 0L) {
    stop("`n` must be a positive integer.", call. = FALSE)
  }

  # if theta is not provided, ensure all model parameters are set
  if (is.null(theta)) {
    missing_params <- vapply(model$models, function(m) {
      pars <- m$parameters
      is.null(pars) || any(is.null(pars))
    }, logical(1L))
    if (any(missing_params)) {
      stop("Model parameters must be set (not NULL).", call. = FALSE)
    }
  }

  # if theta not provided, then evaluate get_variance_covariance_matrix_signal for each model in model$models with their respective parameters, and sum the resulting covariance matrices
  if (is.null(theta)) {
    cov_list <- lapply(model$models, function(m) {
      pars <- m$parameters
      if (is.null(pars) || any(is.null(pars))) {
        stop("Model parameters must be set (not NULL).", call. = FALSE)
      }
      do.call(m$get_variance_covariance_matrix_signal, c(as.list(pars), list(n = n)))
    })
    cov_mat = Reduce(`+`, cov_list)
    return(cov_mat)
  }else{

    # throw an error if prep is null, because we need it to know how to extract the parameters for each component from theta
    if(is.null(prep)){
     stop("`prep` must be provided when `theta` is provided, as it contains the layout information to extract parameters for each component from `theta`.", call. = FALSE)
    }

    # theta is provided and prep is provided so we construct the covariance matrix for each model in model$models with the respective parameters from theta, and sum the resulting covariance matrices
    # theta is assumed to be already in the domain of the parameters
    number_of_components = length(prep$layout)
    cov_mat = matrix(0, nrow = n, ncol = n)
    for(component in seq(number_of_components)){

      # extract index of the parameters in theta vector
      index_param_in_theta_for_component = prep$layout[[component]]$idx
      # extract parameter names for this component
      param_names_for_component = prep$layout[[component]]$pnames
      # extract from theta and give name
      theta_component = theta[index_param_in_theta_for_component]
      names(theta_component) = param_names_for_component
      # Evaluate this component's autocovariance in DOMAIN
      cov_mat_component <- do.call(model$models[[component]]$get_variance_covariance_matrix_signal, c(as.list(theta_component), list(n = n)))

      # Add it to the sum
      cov_mat <- cov_mat + cov_mat_component

    }
     return(cov_mat)
  }
}


# # do some test
# mat1  = get_variance_covariance_matrix_model(model = ar1(.5, 1) +rw(1), n = 10)
# mat2  =get_variance_covariance_matrix_model(model = ar1() +rw(), n = 10, theta = c(.5,1,1), prep = prepare_optim_layout(ar1() +rw()))
# all.equal(mat1, mat2)
# mat3 = get_variance_covariance_matrix_model(model = ar1(phi = .5,2) +flicker(1), n = 10)
# mat4 = get_variance_covariance_matrix_model(model = ar1() +flicker(), n = 10, theta = c(.5,2,1), prep = prepare_optim_layout(ar1() +flicker()))
# all.equal(mat3, mat4)
# # verify some combination to be sure
# mat5  = get_variance_covariance_matrix_model(model = ar1(.5, 1) +rw(1), n = 10)
# sigma2 = 1
# phi = .5
# n=10
# vec_autocov = sigma2 * (phi^(0:(n - 1))) / (1 - phi^2)
# mat6 = fast_toeplitz_matrix_from_vector_cpp(vec_autocov) + get_sigma_mat_rw(n, 1)
# all.equal(mat5, mat6)










#' GMWMX estimator
#'
#' Draft interface for a future `gmwmx2()` that can be called either with a
#' \code{gnss_ts_ngl} object (current workflow) or with a generic design matrix
#' and response vector.
#'
#' @param X Design matrix for a generic regression interface.
#' @param y Response vector for a generic regression interface.
#' @param model Stochastic model specification.
#' @param omega Optional weighting matrix. If `NULL`, uses inverse CI width.
#' @param method Optimization method passed to `stats::optim`.
#' @param control Control list passed to `stats::optim`.
#' @param ... Reserved for future extensions.
#' @return A fitted model object (to be defined).
#' @keywords internal
gmwmx2_new_no_missing <- function(X = NULL, y = NULL, model = NULL, omega = NULL, method = "L-BFGS-B", control = list(), ...) {

  # record start time
  start_time <- Sys.time()

  # Placeholder for the actual implementation
  if (is.null(X) || is.null(y)) {
    stop("Both X and y must be provided for the generic regression interface.")
  }

  # check model
  if (is.null(model) || (!inherits(model, "time_series_model") && !inherits(model, "sum_model"))) {
    stop("`model` must be a 'time_series_model' or 'sum_model'.", call. = FALSE)
  }

  # check that y length matches number of rows in X
  if (length(y) != nrow(X)) {
    stop("`y` length must match the number of rows in `X`.", call. = FALSE)
  }

  # check that there are no NA in X, stop if so
  if (any(is.na(X))) {
    stop("Design matrix X contains NA values. Please remove or impute missing data.")
  }

  # check if there are no NA in y, stop if so
  if (any(is.na(y))) {
    stop("Response vector y contains NA values. Please remove or impute missing data.")
  }


  # get dimension of X and y
  n = nrow(X)
  p = ncol(X)

  # obtain beta hat
  beta_hat <- .lm.fit(y = y, x = X)$coefficients

  # obtain epsilon hat
  eps_hat <- y - X %*% beta_hat

  # compute (X^TX)^{-1} using QR decomposition for numerical stability
  X_transpose = t(X)
  XtX <- X_transpose %*% X
  qr_decomp <- qr(X)
  R <- qr.R(qr_decomp)
  R_inv <- Matrix::solve(R)
  inv_XtX <- R_inv %*% t(R_inv)

  # compute hat matrix
  H <- X %*% inv_XtX %*% X_transpose
  D <- diag(n) - H

  # pre-compute quantities on D=(I-H), with H = X(X^TX)^{-1}X^T to later use in optimization function
  quantities_D <- pre_compute_quantities_on_D_only_required_smarter_cpp(D, approx_type = "3")

  # fill missing parameters in model
  model <- fill_missing_parameters(model, signal = eps_hat)

  # prepare optim layout
  prep <- prepare_optim_layout(model)

  # compute empirical wv on estimated residuals
  wv_emp <- wv::wvar(eps_hat)

  # perform optimization to estimate stochastic parameters
  res <- optim(
    par = prep$theta0,
    fn = loss_fn_gmwmx_no_missing,
    model = model,
    n = n,
    prep = prep,
    quantities_D = quantities_D,
    wv_obj = wv_emp,
    omega = omega,
    method = method,
    control = control,
    ...

  )

  # transform estimated parameters to domain
  theta_domain <- theta_to_domain(model, res$par, prep = prep)
  theta_init_domain <- theta_to_domain(model, prep$theta0, prep = prep)

  # flatten theta domain to a vector
  theta_domain_vec <- unlist(theta_domain)

  # get variance covariance matrix of epsilon hat with model at estimated parameters
  variance_covariance_mat_epsilon <- get_variance_covariance_matrix_model(model, n, theta = theta_domain_vec, prep = prep)

  # construct variance covariance of beta hat with model at estimated parameters
  variance_covariance_beta_hat = inv_XtX %*% X_transpose %*% variance_covariance_mat_epsilon %*% X %*% inv_XtX

  std_beta_hat = sqrt(diag(variance_covariance_beta_hat))

  # construct output
  out = list(
    beta_hat = beta_hat,
    std_beta_hat = std_beta_hat,
    theta_domain = theta_domain,
    convergence = res$convergence,
    value = res$value,
    model = model,
    run_time_sec = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  )

  # assign class to output
  class(out) = "gmwmx2_fit"
  return(out)
}

#' Print method for gmwmx2_fit
#'
#' Displays a table of regression coefficients with standard errors and
#' summarizes the fitted stochastic model with estimated parameters.
#'
#' @param x A `gmwmx2_fit` object.
#' @param digits Significant digits to display.
#' @param ... Passed to print methods.
#' @return The input object, invisibly.
#' @export
print.gmwmx2_fit <- function(x, digits = 4, ...) {
  cat("GMWMX fit\n")

  # regression table
  coef_names <- names(x$beta_hat)
  if (is.null(coef_names) || any(coef_names == "")) {
    coef_names <- paste0("beta", seq_along(x$beta_hat))
  }
  coef_tab <- data.frame(
    Estimate = as.numeric(x$beta_hat),
    Std.Error = as.numeric(x$std_beta_hat),
    row.names = coef_names
  )
  print(coef_tab, digits = digits, ...)

  if (!is.null(x$missing_params) && !is.null(x$missing_prop)) {
    cat("\nMissingness model\n")
    cat("  Proportion missing :", formatC(x$missing_prop, digits = digits, format = "f"), "\n")
    cat("  p1                 :", formatC(x$missing_params$p1, digits = digits, format = "f"), "\n")
    cat("  p2                 :", formatC(x$missing_params$p2, digits = digits, format = "f"), "\n")
    cat("  p*                 :", formatC(x$missing_params$pstar, digits = digits, format = "f"), "\n")
  }

  cat("\nStochastic model\n")
  if (inherits(x$model, "time_series_model")) {
    cat("  Model      :", x$model$model, "\n")
    pars <- if (!is.null(x$theta_domain)) x$theta_domain else x$model$parameters
    cat("  Parameters :", format_params(pars, digits = digits), "\n")
  } else if (inherits(x$model, "sum_model")) {
    n <- length(x$model$models)
    cat("  Sum of", n, "processes\n")
    for (i in seq_along(x$model$models)) {
      m <- x$model$models[[i]]
      pars <- m$parameters
      if (!is.null(x$theta_domain) && is.list(x$theta_domain) && length(x$theta_domain) >= i) {
        pars <- x$theta_domain[[i]]
      }
      cat(sprintf("  [%d] %s\n", i, m$model))
      cat("       Estimated parameters :", format_params(pars, digits = digits), "\n")
    }
  }

  if (!is.null(x$run_time_sec)) {
    cat("\nRuntime (seconds)\n")
    cat("  Total              :", formatC(x$run_time_sec, digits = digits, format = "f"), "\n")
  }

  invisible(x)
}

#' Print method for gmwmx2_fit_gnss_ts_ngl
#'
#' Displays regression coefficients with standard errors and confidence
#' intervals, along with the fitted stochastic and missingness models.
#'
#' @param x A `gmwmx2_fit_gnss_ts_ngl` object.
#' @param digits Significant digits to display.
#' @param scale_parameters Logical; if TRUE, scale trend-related parameters to
#'   m/year (multiplies by 365.25).
#' @param ... Passed to print methods.
#' @return The input object, invisibly.
#' @export
print.gmwmx2_fit_gnss_ts_ngl <- function(x, digits = 4, ...) {
  cat("GMWMX fit (GNSS)\n")

  # regression table with CI
  coef_names <- names(x$beta_hat)
  if (is.null(coef_names) || any(coef_names == "")) {
    coef_names <- paste0("beta", seq_along(x$beta_hat))
  }


  est <- as.numeric(x$beta_hat)
  se <- as.numeric(x$std_beta_hat)
  ci_low <- est - qnorm(0.975) * se
  ci_high <- est + qnorm(0.975) * se

  coef_tab <- data.frame(
    Estimate = est,
    Std.Error = se,
    CI.Lower = ci_low,
    CI.Upper = ci_high,
    row.names = coef_names
  )
  print(coef_tab, digits = digits, ...)

  if (!is.null(x$missing_params) && !is.null(x$missing_prop)) {
    cat("\nMissingness model\n")
    cat("  Proportion missing :", formatC(x$missing_prop, digits = digits, format = "f"), "\n")
    cat("  p1                 :", formatC(x$missing_params$p1, digits = digits, format = "f"), "\n")
    cat("  p2                 :", formatC(x$missing_params$p2, digits = digits, format = "f"), "\n")
    cat("  p*                 :", formatC(x$missing_params$pstar, digits = digits, format = "f"), "\n")
  }

  cat("\nStochastic model\n")
  if (inherits(x$model, "time_series_model")) {
    cat("  Model      :", x$model$model, "\n")
    pars <- if (!is.null(x$theta_domain)) x$theta_domain else x$model$parameters
    cat("  Parameters :", format_params(pars, digits = digits), "\n")
  } else if (inherits(x$model, "sum_model")) {
    n <- length(x$model$models)
    cat("  Sum of", n, "processes\n")
    for (i in seq_along(x$model$models)) {
      m <- x$model$models[[i]]
      pars <- m$parameters
      if (!is.null(x$theta_domain) && is.list(x$theta_domain) && length(x$theta_domain) >= i) {
        pars <- x$theta_domain[[i]]
      }
      cat(sprintf("  [%d] %s\n", i, m$model))
      cat("       Estimated parameters :", format_params(pars, digits = digits), "\n")
    }
  }

  if (!is.null(x$run_time_sec)) {
    cat("\nRuntime (seconds)\n")
    cat("  Total              :", formatC(x$run_time_sec, digits = digits, format = "f"), "\n")
  }

  invisible(x)
}


#' Plot a \code{gmwmx2_fit_gnss_ts_ngl} object
#'
#'
#' @param x A \code{gmwmx2_fit_gnss_ts_ngl} object.
#' @param ... Additional graphical parameters.
#' @return No return value. Plot a \code{gmwmx2_fit_gnss_ts_ngl} object.
#' @export
plot.gmwmx2_fit_gnss_ts_ngl <- function(x, ...) {

  # Save the current graphical parameters
  old_par <- par(no.readonly = TRUE)



  plot.new()
  title(main = "gmwmx2_fit_gnss_ts_ngl plot (draft)")
  invisible(x)

  # Restore the original graphical parameters
  par(old_par)

  # Reset to a single plot layout
  layout(1)

}




#' GMWMX estimator with missing
#'
#' Draft interface for a future `gmwmx2()` that can be called either with a
#' \code{gnss_ts_ngl} object (current workflow) or with a generic design matrix
#' and response vector.
#'
#' @param X Design matrix for a generic regression interface.
#' @param y response vector for a generic regression interface.
#' @param model Stochastic model specification.
#' @param omega Optional weighting matrix. If `NULL`, uses inverse CI width.
#' @param method Optimization method passed to `stats::optim`.
#' @param control Control list passed to `stats::optim`.
#' @param ... Reserved for future extensions.
#' @return A fitted model object (to be defined).
#' @keywords internal
gmwmx2_new_with_missing <- function(X = NULL, y = NULL, model = NULL, omega = NULL, method = "L-BFGS-B", control = list(), ...) {

  # record start time
  start_time <- Sys.time()



  # Placeholder for the actual implementation
  if (is.null(X) || is.null(y)) {
    stop("Both X and y must be provided for the generic regression interface.")
  }

  # check model
  if (is.null(model) || (!inherits(model, "time_series_model") && !inherits(model, "sum_model"))) {
    stop("`model` must be a 'time_series_model' or 'sum_model'.", call. = FALSE)
  }

  # check that y length matches number of rows in X
  if (length(y) != nrow(X)) {
    stop("`y` length must match the number of rows in `X`.", call. = FALSE)
  }

  # check that there are no NA in X, stop if so
  if (any(is.na(X))) {
    stop("Design matrix X contains NA values. Please remove or impute missing data.")
  }

  # require missing values in y for this implementation
  if (!any(is.na(y))) {
    stop("`y` contains no missing values. Use `gmwmx2_new_no_missing()`.", call. = FALSE)
  }

  # get dimension of X and y
  n = nrow(X)
  p = ncol(X)

  # identify missing observation
  vec_is_present = !is.na(y)
  vec_is_present = as.numeric(vec_is_present)
  id_non_missing = which(!is.na(y))

  # obtain X sub and y sub
  X_sub = X[id_non_missing, ]
  y_sub = y[id_non_missing]

  # obtain beta hat
  beta_hat <- .lm.fit(y = y_sub, x = X_sub)$coefficients

  # obtain residuals when we have data
  eps_hat_sub = y_sub - X_sub %*% beta_hat

  # fill in vector with zero when we have missing values
  eps_hat_filled <- vector(mode = "numeric", length = n)

  # fill in residuals where we have data, otherwise zero
  eps_hat_filled[id_non_missing] = eps_hat_sub

  # compute empirical wv on this vector filled with zero
  wv_emp = wv::wvar(eps_hat_filled)

  # estimate missing data mechanism parameters
  p_hat = estimate_p1_p2_mle_cpp(vec_is_present)

  # define pstar hat (expecation of missingness process)
  pstar_hat <- p_hat[2] / (p_hat[1] + p_hat[2])
  missing_prop <- mean(is.na(y))

  # get vec autocovariance theo omega
  vec_autocov_omega <- create_vec_theo_autocov_omega_cpp(p1 = p_hat[1], p2 = p_hat[2], n)

  # compute (X^TX)^{-1} using QR decomposition for numerical stability
  X_transpose = t(X)
  XtX <- X_transpose %*% X
  qr_decomp <- qr(X)
  R <- qr.R(qr_decomp)
  R_inv <- Matrix::solve(R)
  inv_XtX <- R_inv %*% t(R_inv)

  # compute hat matrix
  H <- X %*% inv_XtX %*% X_transpose
  D <- diag(n) - H

  # pre-compute quantities on D=(I-H), with H = X(X^TX)^{-1}X^T to later use in optimization function
  quantities_D <- pre_compute_quantities_on_D_only_required_smarter_cpp(D, approx_type = "3")

  # fill missing parameters in model
  model <- fill_missing_parameters(model, signal = eps_hat_filled)

  # prepare optim layout
  prep <- prepare_optim_layout(model)

  # perform optimization to estimate stochastic parameters
  res <- optim(
    par = prep$theta0,
    fn = loss_fn_gmwmx_with_missing,
    model = model,
    n = n,
    prep = prep,
    quantities_D = quantities_D,
    vec_autocov_omega = vec_autocov_omega,
    pstar_hat = pstar_hat,
    wv_obj = wv_emp,
    omega = omega,
    method = method,
    control = control,
    ...

  )


  # transform estimated parameters to domain
  theta_domain <- theta_to_domain(model, res$par, prep = prep)
  theta_init_domain <- theta_to_domain(model, prep$theta0, prep = prep)

  # flatten theta domain to a vector
  theta_domain_vec <- unlist(theta_domain)

  # get variance covariance matrix of epsilon hat with model at estimated parameters
  variance_covariance_mat_epsilon <- get_variance_covariance_matrix_model(model, n, theta = theta_domain_vec, prep = prep)

  # compute variance covariance of Z (missingness process)
  var_cov_omega <- fast_toeplitz_matrix_from_vector_cpp(as.vector(vec_autocov_omega))

  # compute variance covariance of beta hat
  variance_covariance_beta_hat <- pstar_hat^(-2) * inv_XtX %*% X_transpose %*% ((var_cov_omega + pstar_hat^2) * variance_covariance_mat_epsilon) %*% X %*% inv_XtX

  # get std of beta hat
  std_beta_hat = sqrt(diag(variance_covariance_beta_hat))


  # construct output
  out = list(
    beta_hat = beta_hat,
    std_beta_hat = std_beta_hat,
    theta_domain = theta_domain,
    convergence = res$convergence,
    value = res$value,
    model = model,
    missing_params = list(p1 = p_hat[1], p2 = p_hat[2], pstar = pstar_hat),
    missing_prop = missing_prop,
    run_time_sec = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  )

  # assign class to output
  class(out) = "gmwmx2_fit"
  return(out)
}


#' GMWMX estimator
#'
#' Dispatches either to the generic regression interface (design matrix + response)
#' or to a \code{gnss_ts_ngl} workflow.
#'
#' @param X Either a design matrix (generic regression interface) or a
#'   \code{gnss_ts_ngl} object (GNSS time-series interface).
#' @param ... Additional arguments forwarded to the selected method.
#' @return A fitted model object.
#' @export
gmwmx2_new <- function(X, ...) {
  UseMethod("gmwmx2_new")
}

#' GMWMX estimator
#'
#' Convenience wrapper that selects the missing or non-missing implementation
#' based on the presence of `NA` values in `y`.
#'
#' @param X Design matrix for a generic regression interface.
#' @param y Response vector for a generic regression interface.
#' @param model Stochastic model specification.
#' @param omega Optional weighting matrix. If `NULL`, uses inverse CI width.
#' @param method Optimization method passed to `stats::optim`.
#' @param control Control list passed to `stats::optim`.
#' @param ... Reserved for future extensions.
#' @return A fitted model object.
#' @rdname gmwmx2_new
#' @export
gmwmx2_new.default <- function(X, y, model, omega = NULL, method = "L-BFGS-B", control = list(), ...) {
  if (is.null(X) || is.null(y) || is.null(model)) {
    stop("`X`, `y`, and `model` must be provided.", call. = FALSE)
  }

  # check model
  if (!inherits(model, "time_series_model") && !inherits(model, "sum_model")) {
    stop("`model` must be a 'time_series_model' or 'sum_model'.", call. = FALSE)
  }

  # check that y length matches number of rows in X
  if (length(y) != nrow(X)) {
    stop("`y` length must match the number of rows in `X`.", call. = FALSE)
  }

  if (all(is.na(y))) {
    stop("`y` contains only NA values; cannot fit model.", call. = FALSE)
  }

  if (any(is.na(y))) {
    return(gmwmx2_new_with_missing(
      X = X,
      y = y,
      model = model,
      omega = omega,
      method = method,
      control = control,
      ...
    ))
  }

  gmwmx2_new_no_missing(
    X = X,
    y = y,
    model = model,
    omega = omega,
    method = method,
    control = control,
    ...
  )
}

#' GMWMX estimator for a \code{gnss_ts_ngl} object
#'
#'
#' @param X A \code{gnss_ts_ngl} object (GNSS time-series interface).
#' @param n_seasonal Number of seasonal signals.
#' @param vec_earthquakes_relaxation_time Relaxation time for each earthquake.
#' @param component Component to estimate ("N", "E", or "V").
#' @param model Stochastic model specification.
#' @param ... Reserved for future extensions.
#' @return A fitted model object.
#' @rdname gmwmx2_new
#' @export
gmwmx2_new.gnss_ts_ngl <- function(
  X,
  n_seasonal = 2,
  vec_earthquakes_relaxation_time = NULL,
  component = NULL,
  model = NULL,
  ...
) {


  # #------------------------------------------------------------
  #
  # # with missing values
  # X = gmwmx2::download_station_ngl("1LSU")
  #
  #
  # # without missing values
  # X = gmwmx2::download_station_ngl("4VIR")
  # plot(X)
  #
  # model =wn() +pl()
  # vec_earthquakes_relaxation_time = NULL
  # n_seasonal = 2
  # component = "N"
  # method = "L-BFGS-B"
  # control = list()
  # omega = NULL
  #
  # #------------------------------------------------------------




  # if some required objects are NULL
  if (is.null(X) || is.null(n_seasonal) || is.null(component) || is.null(model)) {
    stop("`X`, `n_seasonal`, `component`, and `model` must be provided.", call. = FALSE)
  }

  # check stochastic model
  if (!inherits(model, "time_series_model") && !inherits(model, "sum_model")) {
    stop("`model` must be a 'time_series_model' or 'sum_model'.", call. = FALSE)
  }

  # Check class
  if (!inherits(X, "gnss_ts_ngl")) {
    stop("Argument `x` should be a `gnss_ts_ngl` object")
  }

  # check that component is either N, E or V
  if (!component %in% c("N", "E", "V")) {
    stop("Argument `component` should take either value `N` or `E` or `V`")
  }

  # check that n_seasonal is either 1 or 2
  if (!n_seasonal %in% c(1, 2)) {
    stop("Argument `n_seasonal` should take either value `1` or `2`")
  }


  #------------------------------------------------------------ START

  # get start time
  start_time <- Sys.time()


  # create full index of MJD values
  all_mjd_index <- seq(head(X$df_position$modified_julian_day, 1), tail(X$df_position$modified_julian_day, 1), by = 1)

  # set n
  n = length(all_mjd_index)

  # create all jumps by combining jumps due to equipment change and jumps due to earthquakes
  jumps <- c(
    X$df_equipment_software_changes$modified_julian_date,
    X$df_earthquakes$modified_julian_date
  )

  # if multiple jumps  to prevent not invertible matrix
  jumps <- unique(jumps)

  # set times where earthquakes happen
  vec_earthquakes_index_mjd <- c(X$df_earthquakes$modified_julian_date)

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

  # if after removing jumps or earthquake that are indicated after the last date, set to NULL
  if(length(vec_earthquakes_index_mjd) == 0){
    vec_earthquakes_index_mjd = NULL
  }

  # if after removing jumps or earthquake that are indicated after the last date, set to NULL
  if(length(jumps) == 0){
    jumps = NULL
  }


  # create design matrix
  X_mat <- create_X_matrix(
    all_mjd_index = all_mjd_index,
    jumps = jumps,
    n_seasonal = n_seasonal,
    vec_earthquakes_index_mjd = vec_earthquakes_index_mjd,
    vec_earthquakes_relaxation_time = vec_earthquakes_relaxation_time
  )

  # Extract Y given specified component
  if (component == "N") {
    y_sub <- X$df_position$northings_fractional_portion
  } else if (component == "E") {
    y_sub <- X$df_position$eastings_fractional_portion
  } else if (component == "V") {
    y_sub <- X$df_position$vertical_fractional_portion
  }

  # identify which of all_mjd index is present and subset X and y
  id_mjd_present = which(all_mjd_index %in% X$df_position$modified_julian_day)
  vec_presence = as.numeric(all_mjd_index %in% X$df_position$modified_julian_day)

  # at this step we can change the logic depending if we have missing observation or not on the signal
  # with missing observations
  if(!all(vec_presence == 1)){

    # obtain X sub
    X_sub = X_mat[id_mjd_present, ]

    # obtain beta hat
    beta_hat <- .lm.fit(y = y_sub, x = X_sub)$coefficients

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
    names(beta_hat) = names_beta_hat

    # obtain residuals when we have data
    eps_hat_sub = y_sub - X_sub %*% beta_hat

    # fill in vector with zero when we have missing values
    eps_hat_filled <- vector(mode = "numeric", length = n)

    # fill in residuals where we have data, otherwise zero
    eps_hat_filled[id_mjd_present] = eps_hat_sub

    # compute empirical wv on this vector filled with zero
    wv_emp = wv::wvar(eps_hat_filled)

    # estimate missing data mechanism parameters
    p_hat = estimate_p1_p2_mle_cpp(vec_presence)

    # define pstar hat (expecation of missingness process)
    pstar_hat <- p_hat[2] / (p_hat[1] + p_hat[2])

    # get vec autocovariance theo omega
    vec_autocov_omega <- create_vec_theo_autocov_omega_cpp(p1 = p_hat[1], p2 = p_hat[2], n)

    # compute (X^TX)^{-1} using QR decomposition for numerical stability
    X_transpose = t(X_mat)
    XtX <- X_transpose %*% X_mat
    qr_decomp <- qr(X_mat)
    R <- qr.R(qr_decomp)
    R_inv <- Matrix::solve(R)
    inv_XtX <- R_inv %*% t(R_inv)

    # compute hat matrix
    H <- X_mat %*% inv_XtX %*% X_transpose
    D <- diag(n) - H

    # pre-compute quantities on D=(I-H), with H = X(X^TX)^{-1}X^T to later use in optimization function
    quantities_D <- pre_compute_quantities_on_D_only_required_smarter_cpp(D, approx_type = "3")

    # fill missing parameters in model
    model <- fill_missing_parameters(model, signal = eps_hat_filled)

    # prepare optim layout
    prep <- prepare_optim_layout(model)

    # perform optimization to estimate stochastic parameters
    res <- optim(
      par = prep$theta0,
      fn = loss_fn_gmwmx_with_missing,
      model = model,
      n = n,
      prep = prep,
      quantities_D = quantities_D,
      vec_autocov_omega = vec_autocov_omega,
      pstar_hat = pstar_hat,
      wv_obj = wv_emp,
      omega = omega,
      method = method,
      control = control,
      ...
    )

    # transform estimated parameters to domain
    theta_domain <- theta_to_domain(model, res$par, prep = prep)
    theta_init_domain <- theta_to_domain(model, prep$theta0, prep = prep)

    # flatten theta domain to a vector
    theta_domain_vec <- unlist(theta_domain)

    # get variance covariance matrix of epsilon hat with model at estimated parameters
    variance_covariance_mat_epsilon <- get_variance_covariance_matrix_model(model, n, theta = theta_domain_vec, prep = prep)

    # # compute variance covariance of Z (missingness process)
    var_cov_omega <- fast_toeplitz_matrix_from_vector_cpp(as.vector(vec_autocov_omega))
    #
    # # compute variance covariance of beta hat
    variance_covariance_beta_hat <- pstar_hat^(-2) * inv_XtX %*% X_transpose %*% ((var_cov_omega + pstar_hat^2) * variance_covariance_mat_epsilon) %*% X_mat %*% inv_XtX
    #
    # # get std of beta hat
    std_beta_hat = sqrt(diag(variance_covariance_beta_hat))
    #
    #

    missing_prop <- mean(is.na(y))


    # # construct output
    out = list(
      beta_hat = beta_hat,
      std_beta_hat = std_beta_hat,
      theta_domain = theta_domain,
      convergence = res$convergence,
      value = res$value,
      model = model,
      missing_params = list(p1 = p_hat[1], p2 = p_hat[2], pstar = pstar_hat),
      missing_prop = missing_prop,
      run_time_sec = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    )

  # without missing observations (all mjd are in y)
  }else{

    # obtain beta hat
    beta_hat <- .lm.fit(y = y_sub, x = X_mat)$coefficients

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
    names(beta_hat) = names_beta_hat

    # obtain epsilon hat
    eps_hat <- y_sub - X_mat %*% beta_hat

    # compute (X^TX)^{-1} using QR decomposition for numerical stability
    X_transpose = t(X_mat)
    XtX <- X_transpose %*% X_mat
    qr_decomp <- qr(X_mat)
    R <- qr.R(qr_decomp)
    R_inv <- Matrix::solve(R)
    inv_XtX <- R_inv %*% t(R_inv)

    # compute hat matrix
    H <- X_mat %*% inv_XtX %*% X_transpose
    D <- diag(n) - H

    # pre-compute quantities on D=(I-H), with H = X(X^TX)^{-1}X^T to later use in optimization function
    quantities_D <- pre_compute_quantities_on_D_only_required_smarter_cpp(D, approx_type = "3")

    # fill missing parameters in model
    model <- fill_missing_parameters(model, signal = eps_hat)

    # prepare optim layout
    prep <- prepare_optim_layout(model)

    # compute empirical wv on estimated residuals
    wv_emp <- wv::wvar(eps_hat)

    # perform optimization to estimate stochastic parameters
    res <- optim(
      par = prep$theta0,
      fn = loss_fn_gmwmx_no_missing,
      model = model,
      n = n,
      prep = prep,
      quantities_D = quantities_D,
      wv_obj = wv_emp,
      omega = omega,
      method = method,
      control = control,
      ...

    )


    # transform estimated parameters to domain
    theta_domain <- theta_to_domain(model, res$par, prep = prep)
    theta_init_domain <- theta_to_domain(model, prep$theta0, prep = prep)

    # flatten theta domain to a vector
    theta_domain_vec <- unlist(theta_domain)

    # get variance covariance matrix of epsilon hat with model at estimated parameters
    variance_covariance_mat_epsilon <- get_variance_covariance_matrix_model(model, n, theta = theta_domain_vec, prep = prep)

    # construct variance covariance of beta hat with model at estimated parameters
    variance_covariance_beta_hat = inv_XtX %*% X_transpose %*% variance_covariance_mat_epsilon %*% X_mat %*% inv_XtX

    std_beta_hat = sqrt(diag(variance_covariance_beta_hat))

    # construct output
    out = list(
      beta_hat = beta_hat,
      std_beta_hat = std_beta_hat,
      theta_domain = theta_domain,
      convergence = res$convergence,
      value = res$value,
      model = model,
      run_time_sec = as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    )
  }

  # assign class to out
  class(out) = "gmwmx2_fit_gnss_ts_ngl"

  return(out)

}
