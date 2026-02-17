# -------------------------- get_autocovariance with optional theta --------------------------
#
# Goal:
#   - If theta is NULL: compute autocovariance using the DOMAIN parameters stored in the model object.
#   - If theta is provided (REAL/unconstrained parameters used by optim):
#       * transform theta -> DOMAIN using the model's transformation_function
#       * evaluate the model's autocovariance_function at these DOMAIN parameters
#   - Works for:
#       * a single "time_series_model"
#       * a composite "sum_model" = sum of several time_series_model
#
# Important:
#   - For sum_model + theta, we use a precomputed "prep" layout created by prepare_optim_layout(model).
#     This avoids any ambiguity when different components share parameter names (e.g. several sigma2).
#     It also avoids rebuilding/mutating the model at every optim iteration.

#' Compute autocovariance for a model (internal)
#'
#' If `theta` is `NULL`, uses domain parameters stored in the model.
#' If `theta` is provided, it is treated as unconstrained parameters in
#' real space and mapped to the domain via the model's transformation.
#'
#' @param object A `time_series_model` or `sum_model`.
#' @param n Length of autocovariance vector.
#' @param theta Optional real-valued parameter vector for optimization.
#' @param prep Optional output from `prepare_optim_layout` (required for sum models with `theta`).
#' @param ... Passed to methods.
#' @return Numeric vector of autocovariances of length `n`.
#' @keywords internal
get_autocovariance <- function(object, n, theta = NULL, prep = NULL, ...) {
  UseMethod("get_autocovariance")
}

# -------------------------- SINGLE MODEL CASE --------------------------
#' @keywords internal
get_autocovariance.time_series_model <- function(object, n, theta = NULL, prep = NULL, ...) {
  # Basic input checks
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n <= 0L) {
    stop("`n` must be a positive integer.", call. = FALSE)
  }

  # ---------------------------------------------------------------------
  # Case 1: theta is NULL -> use DOMAIN parameters stored in object$parameters
  # ---------------------------------------------------------------------
  # note that only a time series model has parameters, not sum model
  if (is.null(theta)) {
    pars <- object$parameters
    if (is.null(pars) || any(is.null(pars))) {
      stop("Model parameters must be set (not NULL).", call. = FALSE)
    }

    # Evaluate autocovariance in DOMAIN:
    return(do.call(object$autocovariance_function, c(as.list(pars), list(n = n))))
  }

  # ---------------------------------------------------------------------
  # Case 2: theta is provided -> theta is in REAL space (unconstrained)
  #          transform REAL -> DOMAIN, then evaluate autocovariance
  # ---------------------------------------------------------------------
  theta <- as.numeric(theta)

  # Parameter names expected by this model (e.g., c("kappa","sigma2"))
  # assign names to theta for do.call ordering, even if theta was passed unnamed (we check length below)
  pnames <- names(object$parameters)
  if (is.null(pnames) || any(pnames == "")) {
    stop("Model parameters must be named.", call. = FALSE)
  }

  # If theta is named, reorder it to match the model's parameter order.
  # If theta is unnamed, require correct length and assume correct order.
  if (!is.null(names(theta))) {
    if (!all(pnames %in% names(theta))) {
      stop("Missing theta entries: ",
           paste(setdiff(pnames, names(theta)), collapse = ", "),
           call. = FALSE)
    }
    theta <- theta[pnames]
  } else {
    if (length(theta) != length(pnames)) {
      stop("Unnamed theta must have length ", length(pnames), ".", call. = FALSE)
    }
  }

  # Transform REAL -> DOMAIN using the model-specific transformation function:
  #   transformation_function(phi_real, sigma2_real, ...) -> (phi_domain, sigma2_domain, ...)
  # assume that theta is in the correct order and set names
  names(theta) <- pnames
  dom <- do.call(object$transformation_function, as.list(theta))
  dom <- as.numeric(dom)
  # assign names to dom for do.call ordering in autocovariance_function, even if transformation_function returned unnamed
  names(dom) <- pnames

  # Evaluate autocovariance in DOMAIN with transformed parameters
  do.call(object$autocovariance_function, c(as.list(dom), list(n = n)))
}





# -------------------------- SUM MODEL CASE --------------------------
#' @keywords internal
get_autocovariance.sum_model <- function(object, n, theta = NULL, prep = NULL, return_components = FALSE, ...) {



  # Basic input checks
  n <- as.integer(n)
  if (length(n) != 1L || is.na(n) || n <= 0L) {
    stop("`n` must be a positive integer.", call. = FALSE)
  }

  # ---------------------------------------------------------------------
  # Case 1: theta is NULL -> use DOMAIN parameters stored in each component model
  #          total autocovariance = sum of component autocovariances
  # ---------------------------------------------------------------------
  if (is.null(theta)) {
    # Compute each component autocovariance using its stored DOMAIN params
    acs <- lapply(object$models, get_autocovariance, n = n, theta = NULL, prep = NULL, ...)
    # Sum them elementwise
    out <- Reduce(`+`, acs)
    if (isTRUE(return_components)) {
      return(list(sum = out, components = acs))
    }
    return(out)
  }

  # ---------------------------------------------------------------------
  # Case 2: theta is provided -> theta is in REAL space for ALL components
  #
  # For speed + correctness (no name collisions), we REQUIRE a prep layout from:
  #   prep <- prepare_optim_layout(object)
  #
  # prep$layout tells us, for each component i:
  #   - info$idx   : indices in theta belonging to this component
  #   - info$pnames: the (unprefixed) parameter names of this component
  # ---------------------------------------------------------------------
  if (is.null(prep) || prep$kind != "sum" || is.null(prep$layout)) {
    stop("For sum_model with theta, you must provide `prep <- prepare_optim_layout(model)`.",
         call. = FALSE)
  }

  theta <- as.numeric(theta)

  # Initialize summed autocovariance vector
  out <- numeric(n)
  if (isTRUE(return_components)) {
    acs <- vector("list", length(prep$layout))
  }

  # Loop over components and add their autocovariances
  number_of_components = length(prep$layout)


  # loop over components using the prep layout to extract parameters for each component from the big theta vector, transform to domain, evaluate autocovariance, and sum
  for(component in seq(number_of_components)){

    # extract index of the parameters in theta vector
    index_param_in_theta_for_component = prep$layout[[component]]$idx
    # extract parameter names for this component
    param_names_for_component = prep$layout[[component]]$pnames
    # extract from theta and give name
    theta_component = theta[index_param_in_theta_for_component]
    names(theta_component) = param_names_for_component
    # transform to domain parameters
    domain_param_component = do.call(object$models[[component]]$transformation_function, as.list(theta_component))
    # set to numeric
    domain_param_component = as.numeric(domain_param_component)
    # reset names
    names(domain_param_component) = param_names_for_component
    # Evaluate this component's autocovariance in DOMAIN
    ac_i <- do.call(object$models[[component]]$autocovariance_function, c(as.list(domain_param_component), list(n = n)))

    # Add it to the sum
    out <- out + ac_i
    if (isTRUE(return_components)) {
      acs[[component]] <- ac_i
    }
  }



  if (isTRUE(return_components)) {
    return(list(sum = out, components = acs))
  }
  out
}


# ----------------------------------------- test it
# mod = wn(sigma2=1)
# get_autocovariance(mod, n=10, theta = log(1))
# mod = wn(sigma2 = 1) + pl(kappa = -.9, sigma2 = 2)
# #
# get_autocovariance(mod, n=10)
# get_autocovariance(mod, n=10, theta = c(log(1), inv_trans_kappa_pl(-.9), log(1.9)), prep = prepare_optim_layout(mod), return_components = T)
#
