# Prefix for component i in a sum_model
#' Component name prefix for summed models
#'
#' @param i Component index.
#' @return Prefix string like `m1_`.
#' @keywords internal
.comp_prefix <- function(i) paste0("m", i, "_")

#' Prepare optimization layout for a model
#'
#' Creates a layout to map a flat parameter vector `theta` (real space)
#' into component parameters for a `sum_model`.
#'
#' @param model A `time_series_model` or `sum_model`.
#' @return A list with `kind`, `theta0`, and (for sums) a `layout`.
#' @keywords internal
#' @examples
#' mod <- wn(1) + pl(kappa = 0.5, sigma2 = 2)
prepare_optim_layout <- function(model) {

  if (inherits(model, "time_series_model")) {
    theta0 <- do.call(model$inv_transformation_function, as.list(model$parameters))
    theta0 <- as.numeric(theta0)
    names(theta0) <- names(model$parameters)

    return(list(kind = "single", theta0 = theta0, pnames = names(model$parameters)))
  }

  if (inherits(model, "sum_model")) {
    theta0 <- c()
    layout <- vector("list", length(model$models))
    idx <- 1L

    for (i in seq_along(model$models)) {
      m <- model$models[[i]]
      pnames <- names(m$parameters)
      k <- length(pnames)

      th_i <- do.call(m$inv_transformation_function, as.list(m$parameters))
      th_i <- as.numeric(th_i)
      names(th_i) <- paste0(.comp_prefix(i), pnames)

      theta0 <- c(theta0, th_i)

      layout[[i]] <- list(
        i = i,
        idx = idx:(idx + k - 1L),
        pnames = pnames
      )
      idx <- idx + k
    }

    return(list(kind = "sum", theta0 = theta0, layout = layout))
  }

  stop("model must be a 'time_series_model' or 'sum_model'.", call. = FALSE)
}

## -------------------------- fill missing parameters --------------------------
#' Fill missing model parameters using initial-parameter functions
#'
#' Ensures any NULL / missing parameters are populated from the model's
#' `get_initial_parameters_function(signal)` while preserving any user-provided
#' parameters (assumed to be in the domain).
#'
#' @param model A `time_series_model` or `sum_model`.
#' @param signal Numeric vector used to derive initial parameters.
#' @return Model with all parameters populated.
#' @keywords internal
fill_missing_parameters <- function(model, signal) {
  get_param_names <- function(m) {
    pnames <- names(formals(m$transformation_function))
    if (!is.null(pnames) && length(pnames) > 0L) return(pnames)
    pnames <- names(formals(m$inv_transformation_function))
    if (!is.null(pnames) && length(pnames) > 0L) return(pnames)
    pnames <- names(m$parameters)
    if (!is.null(pnames) && length(pnames) > 0L) return(pnames)
    NULL
  }

  fill_one <- function(m) {
    pnames <- get_param_names(m)
    if (is.null(pnames) || length(pnames) == 0L) {
      stop("Model parameters must be named.", call. = FALSE)
    }

    current <- rep(NA_real_, length(pnames))
    names(current) <- pnames

    pars <- m$parameters
    if (!is.null(pars) && length(pars) > 0L) {
      if (is.null(names(pars)) || any(names(pars) == "")) {
        stop("Model parameters must be named.", call. = FALSE)
      }
      extra <- setdiff(names(pars), pnames)
      if (length(extra) > 0L) {
        stop("Unknown parameter(s): ", paste(extra, collapse = ", "), call. = FALSE)
      }
      shared <- intersect(names(pars), pnames)
      current[shared] <- pars[shared]
    }

    missing <- is.na(current)
    if (any(missing)) {
      init_fn <- m$get_initial_parameters_function
      if (is.null(init_fn) || !is.function(init_fn)) {
        stop("Missing parameters but no `get_initial_parameters_function` available.", call. = FALSE)
      }
      init <- init_fn(signal)
      if (is.null(init) || length(init) == 0L) {
        stop("Initial parameters function returned nothing.", call. = FALSE)
      }
      if (is.null(names(init)) || any(names(init) == "")) {
        stop("Initial parameters must be named.", call. = FALSE)
      }
      needed <- pnames[missing]
      if (!all(needed %in% names(init))) {
        stop("Initial parameters missing: ",
             paste(setdiff(needed, names(init)), collapse = ", "),
             call. = FALSE)
      }
      current[needed] <- init[needed]
    }

    m$parameters <- current
    m
  }

  if (inherits(model, "time_series_model")) {
    return(fill_one(model))
  }

  if (inherits(model, "sum_model")) {
    model$models <- lapply(model$models, fill_one)
    return(model)
  }

  stop("model must be a 'time_series_model' or 'sum_model'.", call. = FALSE)
}


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
  if (is.null(theta)) {
    pars <- object$parameters
    if (is.null(pars) || any(is.null(pars))) {
      stop("Model parameters must be set (not NULL).", call. = FALSE)
    }

    # Evaluate autocovariance in DOMAIN:
    #   autocovariance_function(kappa=?, sigma2=?, n=?)
    return(do.call(object$autocovariance_function, c(as.list(pars), list(n = n))))
  }

  # ---------------------------------------------------------------------
  # Case 2: theta is provided -> theta is in REAL space (unconstrained)
  #          transform REAL -> DOMAIN, then evaluate autocovariance
  # ---------------------------------------------------------------------
  theta <- as.numeric(theta)

  # Parameter names expected by this model (e.g., c("kappa","sigma2"))
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
  dom <- do.call(object$transformation_function, as.list(theta))
  dom <- as.numeric(dom)
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
  for (info in prep$layout) {
    # Component model object
    m <- object$models[[info$i]]

    # Extract this component's REAL parameters from the big theta vector
    th_i <- theta[info$idx]
    names(th_i) <- info$pnames  # assign param names for do.call ordering

    # Transform REAL -> DOMAIN for this component
    dom_i <- do.call(m$transformation_function, as.list(th_i))
    dom_i <- as.numeric(dom_i)
    names(dom_i) <- info$pnames

    # Evaluate this component's autocovariance in DOMAIN
    ac_i <- do.call(m$autocovariance_function, c(as.list(dom_i), list(n = n)))

    # Add it to the sum
    out <- out + ac_i
    if (isTRUE(return_components)) {
      acs[[info$i]] <- ac_i
    }
  }

  if (isTRUE(return_components)) {
    return(list(sum = out, components = acs))
  }
  out
}


# ----------------------------------------- test it
# mod = wn(sigma2=1)
# get_autocovariance(mod, n=10, theta = c(0))
# mod = wn(sigma2 = 1) + pl(kappa = -.9, sigma2 = 2)
#
# get_autocovariance(mod, n=10)
# get_autocovariance(mod, n=10, theta = c(log(1), inv_trans_kappa_pl(-.9), log(2)), prep = prepare_optim_layout(mod))
#

#----------------------------------------------- try procedure to estimate parameters by matching autocovariance

#
# mod = wn(sigma2 = 5) +rw(sigma2=.1)
# x = generate(mod, n = 10000, seed=123)
# x
# prep = prepare_optim_layout(mod)
# # compute empirical wv
# wv_emp = wv::wvar(x)
# plot(wv_emp)

#' Loss function for GMWM2 optimization (internal)
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
loss_fn <- function(theta, model, n, prep, wv_obj, omega = NULL) {
  # compute autocovariance from theta

  autocov_vec <- get_autocovariance(object = model, n = n, theta = theta, prep = prep)
  # compute wv from autocovariance
  theo_wv <- autocovariance_to_wv(autocov_vec, tau = wv_obj$scales)
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


# define function to transform back parameters
# -------------------- Convert optim theta (REAL) -> DOMAIN parameters --------------------

#' Convert optimization parameters to domain parameters
#'
#' Maps a real-valued optimization vector `theta` into the constrained
#' parameter space used by each model component.
#'
#' @param model A `time_series_model` or `sum_model`.
#' @param theta Numeric vector in real space (typically from `optim`).
#' @param prep Optional output from `prepare_optim_layout` for sum models.
#' @return Named numeric vector (single model) or a named list (sum model).
#' @keywords internal
theta_to_domain <- function(model, theta, prep = NULL) {
  theta <- as.numeric(theta)

  # ---- single model ----
  if (inherits(model, "time_series_model")) {
    pnames <- names(model$parameters)
    if (length(theta) != length(pnames)) {
      stop("theta has length ", length(theta), " but expected ", length(pnames), ".", call. = FALSE)
    }

    names(theta) <- pnames
    dom <- do.call(model$transformation_function, as.list(theta))
    dom <- as.numeric(dom)
    names(dom) <- pnames
    return(dom)
  }

  # ---- sum model ----
  if (inherits(model, "sum_model")) {
    if (is.null(prep) || prep$kind != "sum" || is.null(prep$layout)) {
      stop("For sum_model, provide `prep <- prepare_optim_layout(model)`.", call. = FALSE)
    }

    out <- vector("list", length(model$models))
    for (info in prep$layout) {
      m <- model$models[[info$i]]
      th_i <- theta[info$idx]
      names(th_i) <- info$pnames

      dom_i <- do.call(m$transformation_function, as.list(th_i))
      dom_i <- as.numeric(dom_i)
      names(dom_i) <- info$pnames

      out[[info$i]] <- dom_i
    }

    # Name the list entries by model name (and index to keep unique)
    names(out) <- vapply(seq_along(model$models),
                         function(i) paste0(model$models[[i]]$model, "_", i),
                         character(1))
    return(out)
  }

  stop("model must be a 'time_series_model' or 'sum_model'.", call. = FALSE)
}




# transform back obtained parameters
# trans_param = theta_to_domain(mod, res$par, prep = prep)

# get theoretical wv implied by these parameters
#' Theoretical wavelet variance from model parameters
#'
#' Computes the theoretical WV implied by a model and parameter vector.
#'
#' @param theta Real-valued parameter vector (unconstrained).
#' @param model A `time_series_model` or `sum_model`.
#' @param n Length of autocovariance to compute.
#' @param wv_obj Optional `wv::wvar` object (uses its scales).
#' @param tau Optional vector of scales; overrides `wv_obj` if provided.
#' @param prep Optional output from `prepare_optim_layout` for sum models.
#' @return Numeric vector of theoretical wavelet variances.
#' @keywords internal
get_theoretical_wv <- function(theta, model, n, wv_obj = NULL, tau = NULL, prep = NULL, return_components = FALSE) {
  # note that theta should be in the real domain as get autocovariance expects real parameters and does the transformation internally
  # Determine the tau/scales where we want the theoretical WV
  # - Either pass tau directly
  # - Or pass a wv::wvar object and we use its scales
  if (is.null(tau)) {
    if (is.null(wv_obj)) stop("Provide either `tau` or `wv_obj` (with $scales).", call. = FALSE)
    tau <- wv_obj$scales
  }

  if (isTRUE(return_components) && inherits(model, "sum_model")) {
    acov <- get_autocovariance(object = model, n = n, theta = theta, prep = prep, return_components = TRUE)
    wv_sum <- autocovariance_to_wv(acov$sum, tau = tau)
    wv_components <- lapply(acov$components, autocovariance_to_wv, tau = tau)
    return(list(sum = wv_sum, components = wv_components))
  }

  # Compute autocovariance for these REAL parameters
  acov <- get_autocovariance(object = model, n = n, theta = theta, prep = prep)

  # Map autocovariance -> theoretical WV
  autocovariance_to_wv(acov, tau = tau)
}


# theo_wv = get_theoretical_wv(theta = res$par, model = mod, n = 10000, wv_obj = wv_emp, prep = prep)
# lines(wv_emp$scales, theo_wv, col = "red", lwd = 2)



# ----------- define gmwm2 function
# estimate parameters of a stochastic model by matching wavelet variance

#' Estimate composite stochastic model
#'
#' Implements the Generalized Method of Wavelet Moments (GMWM) estimator
#' to fit a `time_series_model` or `sum_model` to data. The estimator
#' minimizes a weighted distance between empirical wavelet variance (WV)
#' and the model-implied WV across scales.
#'
#' @param x Numeric vector, or a `generated_time_series` /
#'   `generated_composite_model_time_series` object (its `series` is used).
#' @param model A `time_series_model` or `sum_model`.
#' @param omega Optional weighting matrix. If `NULL`, a default based on the
#'   empirical WV confidence intervals is used.
#' @param method Optimization method passed to `stats::optim`.
#' @param control Optional list of control parameters for `stats::optim`.
#' @param ... Additional arguments passed to `stats::optim`.
#' @return An object of class `gmwm2_fit` with elements:
#'   `theta_hat` (real space), `theta_domain` (constrained space),
#'   `model`, `empirical_wvar`, `theoretical_wvar`, `optim`, and `n`.
#' @details
#' The GMWM estimator solves a weighted least-squares criterion of the form
#' \deqn{
#'   (\hat{\boldsymbol{\nu}} - \boldsymbol{\nu}(\boldsymbol{\theta}))^\top
#'   \boldsymbol{\Omega}
#'   (\hat{\boldsymbol{\nu}} - \boldsymbol{\nu}(\boldsymbol{\theta}))
#' }
#' where \eqn{\hat{\boldsymbol{\nu}}} denotes the \strong{empirical wavelet
#' variance (WV)} vector and \eqn{\boldsymbol{\nu}(\boldsymbol{\theta})}
#' the corresponding \strong{theoretical WV} vector implied by the model
#' parameters \eqn{\boldsymbol{\theta}}. The weighting matrix
#' \eqn{\boldsymbol{\Omega}} defaults to a diagonal matrix with entries proportional to the
#' inverse squared width of the empirical WV confidence intervals. Provide
#' `omega` to use a custom weighting (e.g., from a theoretical covariance).
#' @references
#' Guerrier, S., Skaloud, J., Stebler, Y., and Victoria-Feser, M.-P. (2013).
#' Wavelet-variance-based estimation for composite stochastic processes.
#' *Journal of the American Statistical Association*, 108(503), 1021-1030.
#' doi:10.1080/01621459.2013.799920.
#' @importFrom wv wvar
#' @importFrom stats optim
#' @examples
#' n = 50000
#' mod = wn(20) + ar1(phi = .9, sigma2 = 1)+ar1(phi = .99, sigma2 = .1)
#' y = generate(mod, n = n)
#' fit = gmwm2(y, model = wn() + ar1(phi = .9) + ar1(phi = .99))
#' plot(fit)
#' @export
gmwm2 <- function(x, model, omega = NULL, method = "L-BFGS-B", control = list(), ...) {
  # unwrap generated_* objects
  if (inherits(x, "generated_time_series") || inherits(x, "generated_composite_model_time_series")) {
    x <- x$series
  }

  if (!is.numeric(x) || !is.vector(x)) {
    stop("`x` must be a numeric vector or a generated time series object.", call. = FALSE)
  }

  if (!inherits(model, "time_series_model") && !inherits(model, "sum_model")) {
    stop("`model` must be a 'time_series_model' or 'sum_model'.", call. = FALSE)
  }

  model <- fill_missing_parameters(model, signal = x)

  n <- length(x)
  prep <- prepare_optim_layout(model)
  wv_emp <- wv::wvar(x)

  res <- optim(
    par = prep$theta0,
    fn = loss_fn,
    model = model,
    n = n,
    prep = prep,
    wv_obj = wv_emp,
    omega = omega,
    method = method,
    control = control,
    ...
  )

  theta_domain <- theta_to_domain(model, res$par, prep = prep)
  theta_init_domain <- theta_to_domain(model, prep$theta0, prep = prep)
  if (inherits(model, "sum_model")) {
    wv_out <- get_theoretical_wv(res$par, model = model, n = n, wv_obj = wv_emp, prep = prep, return_components = TRUE)
    wv_theo <- wv_out$sum
    wv_theo_components <- wv_out$components
    names(wv_theo_components) <- vapply(seq_along(model$models),
                                        function(i) model$models[[i]]$model,
                                        character(1))
  } else {
    wv_theo <- get_theoretical_wv(res$par, model = model, n = n, wv_obj = wv_emp, prep = prep)
    wv_theo_components <- NULL
  }

  out <- list(
    theta_hat = res$par,
    theta_init_domain = theta_init_domain,
    theta_domain = theta_domain,
    model = model,
    empirical_wvar = wv_emp,
    theoretical_wvar = wv_theo,
    theoretical_wvar_components = wv_theo_components,
    optim = res,
    n = n
  )
  class(out) <- "gmwm2_fit"
  out
}

#' Print method for a \code{gmwm2_fit} object
#'
#' @param x A \code{gmwm2_fit} object.
#' @param digits Significant digits for printing.
#' @param show_initial_parameters Logical; if TRUE, also show the initial parameters used for optimization.
#' @param ... Unused.
#' @examples
#' model <- wn(sigma2 = 1) + ar1(phi = 0.8, sigma2 = 0.5)
#' x <- generate(model, n = 1000, seed = 123)
#' plot(x)
#' fit <- gmwm2(x, model = wn()+ar1())
#' fit
#' @return The input object, invisibly.
#' @export
print.gmwm2_fit <- function(x, digits = 4,show_initial_parameters = FALSE, ...) {
  model <- x$model

  format_params_est <- function(pars) {
    if (is.null(pars) || length(pars) == 0L) return("(no parameters)")
    paste0(
      names(pars), " = ",
      formatC(as.numeric(pars), digits = digits, format = "g"),
      collapse = ", "
    )
  }

  cat("GMWM fit\n\n")

  if (inherits(model, "sum_model")) {
    cat("Stochastic model\n")
    cat("  Sum of", length(model$models), "processes\n\n")
    for (i in seq_along(model$models)) {
      m <- model$models[[i]]
      pnames <- names(m$parameters)
      cat(sprintf("  [%d] %s\n", i, m$model))
      cat("      Parameters : ", paste(pnames, collapse = ", "), "\n\n", sep = "")
    }
  } else {
    cat("Stochastic model\n")
    cat("  Model      :", model$model, "\n")
    cat("  Parameters :", paste(names(model$parameters), collapse = ", "), "\n\n")
  }
  if(show_initial_parameters){
    cat("Initial parameters\n")
    if (inherits(model, "sum_model")) {
      init_list <- x$theta_init_domain
      for (i in seq_along(model$models)) {
        m <- model$models[[i]]
        pars <- init_list[[i]]
        cat(sprintf("  %d) %s: %s\n", i, m$model, format_params_est(pars)))
      }
    } else {
      cat(sprintf("  1) %s: %s\n", model$model, format_params_est(x$theta_init_domain)))
    }
  }


  cat("\nEstimated parameters\n")
  if (inherits(model, "sum_model")) {
    dom_list <- x$theta_domain
    for (i in seq_along(model$models)) {
      m <- model$models[[i]]
      pars <- dom_list[[i]]
      cat(sprintf("  %d) %s: %s\n", i, m$model, format_params_est(pars)))
    }
  } else {
    cat(sprintf("  1) %s: %s\n", model$model, format_params_est(x$theta_domain)))
  }

  conv <- x$optim$convergence
  status <- if (is.numeric(conv) && length(conv) == 1L && conv == 0L) "converged" else "not converged"
  loss <- x$optim$value
  loss <- if (is.null(loss)) NA_real_ else as.numeric(loss)
  counts <- x$optim$counts
  iters <- if (is.null(counts)) NA_integer_ else sum(as.integer(counts))

  cat("\nOptimization\n")
  cat("  Convergence : ", status, " (code ", conv, ")\n", sep = "")
  cat("  Iterations  : ", iters, "\n", sep = "")
  cat("  Loss        : ", formatC(loss, digits = digits, format = "g"), "\n", sep = "")

  invisible(x)
}

#' Plot method for a \code{gmwm2_fit} object
#'
#' Plots empirical wavelet variance with the fitted theoretical curve and,
#' for sum models, component-implied theoretical curves.
#'
#' @param x A \code{gmwm2_fit} object.
#' @param show_ci Logical; if TRUE and available, show empirical CI bars.
#' @param col_emp Color for empirical WV points/line.
#' @param col_theo Color for theoretical WV line.
#' @param col_ci Color for empirical WV CI band.
#' @param lwd Line width for theoretical curve.
#' @param pch_emp Plotting character for empirical points.
#' @param pch_theo Plotting character for theoretical points.
#' @param cex_theo Size for theoretical points.
#' @param legend_pos Legend position (e.g., "topleft") or "auto".
#' @param ... Additional arguments passed to `plot()`.
#' @return The input object, invisibly.
#' @examples
#' n = 50000
#' mod = wn(20) + ar1(phi = .9, sigma2 = 1)+ar1(phi = .99, sigma2 = .1)
#' y = generate(mod, n = n)
#' fit = gmwm2(y, model = wn(20) + ar1(phi = .9, sigma2 = 1)+ar1(phi = .99, sigma2 = .1))
#' plot(fit)
#' @export
plot.gmwm2_fit <- function(x,
                           show_ci = TRUE,
                           col_emp = "black",
                           col_theo = "darkorange",
                           col_ci = "#e6f7fb",
                           lwd = 2,
                           pch_emp = 16,
                           pch_theo = 21,
                           cex_theo = 1.4,
                           legend_pos = "auto",
                           ...) {
  format_model_text <- function(model) {
    if (inherits(model, "sum_model")) {
      parts <- vapply(model$models, function(m) m$model, character(1))
      paste(parts, collapse = " + ")
    } else if (inherits(model, "time_series_model")) {
      model$model
    } else {
      "Unknown model"
    }
  }

  pick_legend_pos <- function(scales, var_emp, yl) {
    if (length(scales) == 0L || length(var_emp) == 0L) {
      return("topleft")
    }
    first_val <- var_emp[1]
    if (!is.finite(first_val) || length(yl) != 2L || !all(is.finite(yl))) {
      return("topleft")
    }
    # Use log-scale center of y-axis to decide "low vs high" positioning.
    center_log <- sqrt(yl[1] * yl[2])
    if (first_val <= center_log) "topleft" else "bottomleft"
  }

  wv_emp <- x$empirical_wvar
  if (is.null(wv_emp$scales) || is.null(wv_emp$variance)) {
    stop("Empirical wavelet variance is missing scales/variance.", call. = FALSE)
  }

  scales <- wv_emp$scales
  var_emp <- wv_emp$variance

  theo <- x$theoretical_wvar
  if (is.null(theo) || length(theo) != length(scales)) {
    prep <- prepare_optim_layout(x$model)
    theo <- get_theoretical_wv(
      theta = x$theta_hat,
      model = x$model,
      n = x$n,
      tau = scales,
      prep = prep
    )
  }
  theo_components <- NULL
  if (inherits(x$model, "sum_model")) {
    theo_components <- x$theoretical_wvar_components
    if (is.null(theo_components)) {
      prep <- prepare_optim_layout(x$model)
      wv_out <- get_theoretical_wv(
        theta = x$theta_hat,
        model = x$model,
        n = x$n,
        tau = scales,
        prep = prep,
        return_components = TRUE
      )
      theo_components <- wv_out$components
    }
  }

  if (isTRUE(show_ci) &&
      !is.null(wv_emp$ci_low) &&
      !is.null(wv_emp$ci_high) &&
      length(wv_emp$ci_low) == length(scales) &&
      length(wv_emp$ci_high) == length(scales)) {
    yl <- range(c(wv_emp$ci_low, wv_emp$ci_high, theo), finite = TRUE)
  } else {
    yl <- range(c(var_emp, theo), finite = TRUE)
  }
  if (!is.null(theo_components)) {
    yl <- range(c(yl, unlist(theo_components)), finite = TRUE)
  }

  plot(
    NA,
    ylim = yl,
    xlim = range(scales, finite = TRUE),
    main = "",
    ylab = "",
    xlab = "",
    log = "xy",
    xaxt = "n",
    yaxt = "n",
    ...
  )
  model_text <- format_model_text(x$model)
  mtext(side = 3, text = paste0("Model: ", model_text), line = 0.2)
  mtext(side = 1, text = "Scales", line = 2.3)
  mtext(side = 2, text = "Wavelet Variance", line = 3.2)

  exponents <- log2(scales)
  labels <- sapply(exponents, function(x) as.expression(bquote(2^.(x))))
  axis(1, at = scales, labels = labels)

  if (isTRUE(show_ci) &&
      !is.null(wv_emp$ci_low) &&
      !is.null(wv_emp$ci_high) &&
      length(wv_emp$ci_low) == length(scales) &&
      length(wv_emp$ci_high) == length(scales)) {
    polygon(
      x = c(scales, rev(scales)),
      y = c(wv_emp$ci_low, rev(wv_emp$ci_high)),
      col = col_ci,
      border = NA
    )
  }

  ylab <- floor(log10(yl[1])):ceiling(log10(yl[2]))
  axis(2, at = 10^ylab, labels = parse(text = sprintf("10^%.0f", ylab)), las = 1)

  abline(h = 10^ylab, col = "grey90", lt = 2)
  for (i in seq_along(exponents)) {
    abline(v = 2^exponents[i], col = "grey90", lt = 2)
  }

  lines(scales, var_emp, type = "b", col = col_emp, pch = pch_emp)
  lines(scales, theo, type = "b", col = col_theo, pch = pch_theo, cex = cex_theo, lwd = lwd)
  if (!is.null(theo_components) && length(theo_components) > 0L) {
    comp_cols <- .gmwmx2_get_plot_colors(length(theo_components))
    for (i in seq_along(theo_components)) {
      lines(scales, theo_components[[i]], type = "l", col = comp_cols[i], lwd = 1)
    }
  }

  legend_pos <- if (identical(legend_pos, "auto")) {
    pick_legend_pos(scales, var_emp, yl)
  } else {
    legend_pos
  }
  legend_labels <- c("Empirical WV", "Theoretical WV", "Empirical WV CI")
  legend_cols <- c(col_emp, col_theo, col_ci)
  legend_pch <- c(pch_emp, pch_theo, 15)
  legend_pt_cex <- c(1, cex_theo, 3)
  legend_lty <- c(1, 1, NA)

  if (!is.null(theo_components) && length(theo_components) > 0L) {
    comp_names <- names(theo_components)
    if (is.null(comp_names) || any(comp_names == "")) {
      comp_names <- vapply(seq_along(theo_components),
                           function(i) paste0("Component ", i),
                           character(1))
    }
    comp_cols <- .gmwmx2_get_plot_colors(length(theo_components))
    legend_labels <- c(legend_labels, comp_names)
    legend_cols <- c(legend_cols, comp_cols)
    legend_pch <- c(legend_pch, rep(NA_integer_, length(comp_cols)))
    legend_pt_cex <- c(legend_pt_cex, rep(1, length(comp_cols)))
    legend_lty <- c(legend_lty, rep(1, length(comp_cols)))
  }

  legend(
    legend_pos,
    legend = legend_labels,
    col = legend_cols,
    pch = legend_pch,
    pt.cex = legend_pt_cex,
    lty = legend_lty,
    horiz = FALSE,
    bty = "n",
    bg = "transparent"
  )

  box()

  invisible(x)
}



