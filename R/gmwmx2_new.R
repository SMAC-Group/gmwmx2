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
  objective <- t(difference) %*% omega %*% difference
  return(objective)
  # compute loss
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
#' @param X Optional design matrix for a generic regression interface.
#' @param y Optional response vector for a generic regression interface.
#' @param model Optional stochastic model specification.
#' @param omega Optional weighting matrix. If `NULL`, uses inverse CI width.
#' @param method Optimization method passed to `stats::optim`.
#' @param control Control list passed to `stats::optim`.
#' @param ... Reserved for future extensions.
#' @return A fitted model object (to be defined).
#' @keywords internal
gmwmx2_new_no_missing <- function(X = NULL, y = NULL, model = NULL, omega = NULL, method = "L-BFGS-B", control = list(), ...) {
  #-------------------------------------------
  # n=5000
  # X =matrix(NA, nrow=n, ncol=2)
  # X[,1] = 1
  # X[,2] = 1:n
  # beta = c(1, .2)
  # eps = generate(ar1(phi=0.8, sigma2=20) + wn(20), n=n)$series
  # # plot(wv::wvar(eps))
  # y = X %*% beta + eps
  # plot(X[,2], y, type='l')
  # method = "L-BFGS-B"
  # control = list()
  # omega =NULL
  # model = ar1()+wn()
  #----------------------------------------------


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
    model = model
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

  invisible(x)
}









#
# n = 10000
# X = matrix(NA, nrow=n, ncol=4)
# # intercept
# X[,1] = 1
# # trend
# X[,2] = 1:n
# # add a sin signal
# omega_1 <- (1 / 365.25) * 2 * pi
# X[, 3] <- sin((1:n) * omega_1)
# X[, 4] <- cos((1:n) * omega_1)
# beta = c(1, .2, 3,4)
# yy = X%*% beta
# plot(X[,2], yy, type='l')
# eps = generate(wn(1) + pl(kappa = -.9, sigma2 = 1), n=n, seed = (123 + b))$series
# y = X %*% beta + eps
# fit = gmwmx2_new_no_missing(X = X, y = y, model = wn() + pl() )
# fit


#
# # do a little check, do not remove, to use later for a vignette
# n = 1000
# X = matrix(NA, nrow=n, ncol=4)
# # intercept
# X[,1] = 1
# # trend
# X[,2] = 1:n
# # add a sin signal
# omega_1 <- (1 / 365.25) * 2 * pi
# X[, 3] <- sin((1:n) * omega_1)
# X[, 4] <- cos((1:n) * omega_1)
# beta = c(1, .2, 3,4)
# yy = X%*% beta
# plot(X[,2], yy, type='l')
# B = 500
# mat_res = matrix(NA, nrow=B, ncol=19)
# for(b in seq(B)){
#   eps = generate(ar1(phi=0.95, sigma2=20) + wn(20), n=n, seed = (123 + b))$series
#   # plot(wv::wvar(eps))
#   y = X %*% beta + eps
#   fit = gmwmx2_new_no_missing(X = X, y = y, model = wn() + ar1() )
#   # mispecified model assuming white noise as the stochastic model
#   fit2 = lm(y~X[,2] + X[,3] + X[,4])
#
#   mat_res[b, ] = c(fit$beta_hat, fit$std_beta_hat,
#                    summary(fit2)$coefficients[,1],
#                    summary(fit2)$coefficients[,2],
#                    fit$theta_domain$`AR(1)_2`,
#                    fit$theta_domain$`White Noise_1`)
#   cat("Iteration ", b, " completed.\n")
# }
#
# # compute empirical coverage
# mat_res_df = as.data.frame(mat_res)
# colnames(mat_res_df) = c("gmwmx_beta0_hat", "gmwmx_beta1_hat",
#                          "gmwmx_beta2_hat", "gmwmx_beta3_hat",
#                          "gmwmx_std_beta0_hat", "gmwmx_std_beta1_hat",
#                          "gmwmx_std_beta2_hat", "gmwmx_std_beta3_hat",
#                          "lm_beta0_hat", "lm_beta1_hat", "lm_beta2_hat", "lm_beta3_hat",
#                          "lm_std_beta0_hat", "lm_std_beta1_hat", "lm_std_beta2_hat", "lm_std_beta3_hat",
#                          "phi_ar1","sigma_2_ar1" ,"sigma_2_wn")
# zval = qnorm(0.975)
# mat_res_df$upper_ci_gmwmx_beta0 = mat_res_df$gmwmx_beta0_hat + zval * mat_res_df$gmwmx_std_beta0_hat
# mat_res_df$lower_ci_gmwmx_beta0 = mat_res_df$gmwmx_beta0_hat - zval * mat_res_df$gmwmx_std_beta0_hat
# mat_res_df$upper_ci_gmwmx_beta1 = mat_res_df$gmwmx_beta1_hat + zval * mat_res_df$gmwmx_std_beta1_hat
# mat_res_df$lower_ci_gmwmx_beta1 = mat_res_df$gmwmx_beta1_hat - zval * mat_res_df$gmwmx_std_beta1_hat
# # empirical coverage of gmwmx beta
# dplyr::between(rep(1, 500), mat_res_df$lower_ci_gmwmx_beta0, mat_res_df$upper_ci_gmwmx_beta0) %>% mean()
# dplyr::between(rep(0.2, 500), mat_res_df$lower_ci_gmwmx_beta1, mat_res_df$upper_ci_gmwmx_beta1) %>% mean()
#
# # do the same for lm beta
# mat_res_df$upper_ci_lm_beta0 = mat_res_df$lm_beta0_hat + zval * mat_res_df$lm_std_beta0_hat
# mat_res_df$lower_ci_lm_beta0 = mat_res_df$lm_beta0_hat - zval * mat_res_df$lm_std_beta0_hat
# mat_res_df$upper_ci_lm_beta1 = mat_res_df$lm_beta1_hat + zval * mat_res_df$lm_std_beta1_hat
# mat_res_df$lower_ci_lm_beta1 = mat_res_df$lm_beta1_hat - zval * mat_res_df$lm_std_beta1_hat
# dplyr::between(rep(1, 500), mat_res_df$lower_ci_lm_beta0, mat_res_df$upper_ci_lm_beta0) %>% mean()
# dplyr::between(rep(0.2, 500), mat_res_df$lower_ci_lm_beta1, mat_res_df$upper_ci_lm_beta1) %>% mean()





#' GMWMX estimator with missing
#'
#' Draft interface for a future `gmwmx2()` that can be called either with a
#' \code{gnss_ts_ngl} object (current workflow) or with a generic design matrix
#' and response vector.
#'
#' @param X Optional design matrix for a generic regression interface.
#' @param y Optional response vector for a generic regression interface.
#' @param model Optional stochastic model specification.
#' @param omega Optional weighting matrix. If `NULL`, uses inverse CI width.
#' @param method Optimization method passed to `stats::optim`.
#' @param control Control list passed to `stats::optim`.
#' @param ... Reserved for future extensions.
#' @return A fitted model object (to be defined).
#' @keywords internal
gmwmx2_new_with_missing <- function(X = NULL, y = NULL, model = NULL, omega = NULL, method = "L-BFGS-B", control = list(), ...) {
  #-------------------------------------------
  # n=5000
  # X =matrix(NA, nrow=n, ncol=2)
  # X[,1] = 1
  # X[,2] = 1:n
  # beta = c(1, .2)
  # eps = generate(ar1(phi=0.8, sigma2=20) + wn(20), n=n)$series
  # plot(wv::wvar(eps))
  # y = X %*% beta + eps
  # plot(X[,2], y, type='l')
  # method = "L-BFGS-B"
  # control = list()
  # omega =NULL
  # model = ar1()+wn()
  #----------------------------------------------


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
    model = model
  )

  # assign class to output
  class(out) = "gmwmx2_fit"
  return(out)
}





