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

