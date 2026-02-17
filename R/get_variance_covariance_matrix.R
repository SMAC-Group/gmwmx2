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
  if (is.null(theta)) {
    pars <- model$parameters
    if (is.null(pars) || any(is.null(pars))) {
      stop("Model parameters must be set (not NULL).", call. = FALSE)
    }
    do.call(model$get_variance_covariance_matrix_signal, c(as.list(pars), list(n = n)))
    # else evaluate get_variance_covariance_matrix_signal with theta, but set names to theta according to the order of parameters in model$parameters
  } else {
    # set names to theta according to the order of parameters in model$parameters
    names(theta) <- names(model$parameters)
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
    cov_mat <- Reduce(`+`, cov_list)
    return(cov_mat)
  } else {
    # throw an error if prep is null, because we need it to know how to extract the parameters for each component from theta
    if (is.null(prep)) {
      stop("`prep` must be provided when `theta` is provided, as it contains the layout information to extract parameters for each component from `theta`.", call. = FALSE)
    }

    # theta is provided and prep is provided so we construct the covariance matrix for each model in model$models with the respective parameters from theta, and sum the resulting covariance matrices
    # theta is assumed to be already in the domain of the parameters
    number_of_components <- length(prep$layout)
    cov_mat <- matrix(0, nrow = n, ncol = n)
    for (component in seq(number_of_components)) {
      # extract index of the parameters in theta vector
      index_param_in_theta_for_component <- prep$layout[[component]]$idx
      # extract parameter names for this component
      param_names_for_component <- prep$layout[[component]]$pnames
      # extract from theta and give name
      theta_component <- theta[index_param_in_theta_for_component]
      names(theta_component) <- param_names_for_component
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





