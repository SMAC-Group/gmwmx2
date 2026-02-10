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
  n=1000
  X =matrix(NA, nrow=n, ncol=2)
  X[,1] = 1
  X[,2] = 1:n
  beta = c(1, .2)
  y = X %*% beta + generate(ar1(phi=0.8, sigma2=20) + wn(20), n=n)$series
  plot(X[,2], y, type='l')
  method = "L-BFGS-B"
  control = list()
#----------------------------------------------

  # get dimension of X and y
  n = nrow(X)
  p = ncol(X)

  # Placeholder for the actual implementation
  if (is.null(X) || is.null(y)) {
    stop("Both X and y must be provided for the generic regression interface.")
  }

  # check that there are no NA in X, stop if so
  if (any(is.na(X))) {
    stop("Design matrix X contains NA values. Please remove or impute missing data.")
  }

  # check if there are no NA in y, stop if so
  if (any(is.na(y))) {
    stop("Response vector y contains NA values. Please remove or impute missing data.")
  }

  # obtain beta hat
  beta_hat <- .lm.fit(y = y, x = X)$coefficients

  # obtain epsilon hat
  eps_hat <- y - X %*% beta_hat

  # compute (X^TX)^{-1} using QR decomposition for numerical stability
  XtX <- t(X) %*% X
  qr_decomp <- qr(X)
  R <- qr.R(qr_decomp)
  R_inv <- Matrix::solve(R)
  inv_XtX <- R_inv %*% t(R_inv)

  # compute hat matrix
  H <- X %*% inv_XtX %*% t(X)
  D <- diag(n) - H

  # pre-compute quantities on D=(I-H), with H = X(X^TX)^{-1}X^T to later use in optimization function
  quantities_D <- pre_compute_quantities_on_D_only_required_smarter_cpp(D, approx_type = "3")

  # fill missing parameters in model
  model <- fill_missing_parameters(model, signal = eps_hat)

  # prepare optim layout
  prep <- prepare_optim_layout(model)

  # compute empirical wv on estimated residuals
  wv_emp <- wv::wvar(eps_hat)

  # res <- optim(
  #   par = prep$theta0,
  #   fn = loss_fn_gmwmx_no_missing,
  #   model = model,
  #   n = n,
  #   prep = prep,
  #   wv_obj = wv_emp,
  #   omega = omega,
  #   method = method,
  #   control = control,
  #   ...
  # )

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
    control = control

  )


  theta_domain <- theta_to_domain(model, res$par, prep = prep)
  theta_init_domain <- theta_to_domain(model, prep$theta0, prep = prep)





  # fit a stochastic model to the residuals








}

#' Variance-covariance matrix implied by a model
#'
#' Constructs the variance-covariance matrix for a model using its
#' `get_variance_covariance_matrix_signal` method. Supports both single
#' `time_series_model` objects and composite `sum_model` objects.
#'
#' @param model A `time_series_model` or `sum_model`.
#' @param n Length of the signal.
#' @param theta Optional parameter vector (not supported yet).
#' @param prep Optional output from `prepare_optim_layout(model)` (unused for now).
#' @return Variance-covariance matrix of dimension `n x n`.
#' @keywords internal
get_variance_covariance_matrix_model <- function(model, n, theta = NULL, prep = NULL) {
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n <= 0L) {
    stop("`n` must be a positive integer.", call. = FALSE)
  }

  if (!is.null(theta)) {
    stop("`theta` is not supported yet. Set model parameters and call without `theta`.",
         call. = FALSE)
  }

  if (inherits(model, "time_series_model")) {
    pars <- model$parameters
    if (is.null(pars) || any(is.null(pars))) {
      stop("Model parameters must be set (not NULL).", call. = FALSE)
    }
    return(do.call(model$get_variance_covariance_matrix_signal, c(as.list(pars), list(n = n))))
  }

  if (inherits(model, "sum_model")) {
    cov_list <- lapply(model$models, function(m) {
      pars <- m$parameters
      if (is.null(pars) || any(is.null(pars))) {
        stop("Model parameters must be set (not NULL).", call. = FALSE)
      }
      do.call(m$get_variance_covariance_matrix_signal, c(as.list(pars), list(n = n)))
    })
    return(Reduce(`+`, cov_list))
  }

  stop("`model` must be a 'time_series_model' or 'sum_model'.", call. = FALSE)
}


# get_variance_covariance_matrix_model(model = ar1(.5, 1) +rw(1), n = 10)
# get_variance_covariance_matrix_model(model = ar1() +rw(), n = 10, theta = c(.5,1,1), prep = prepare_optim_layout(ar1() +rw()))
# # phi = .8
# # sigma2 = 1
# # n=10
# fast_toeplitz_matrix_from_vector_cpp(sigma2 * (phi^(0:(n - 1))) / (1 - phi^2)) + get_sigma_mat_rw(n, 1)
