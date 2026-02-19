#' Prepare optimization layout for a model
#'
#' Build a bookkeeping structure that maps a flat optimization vector
#' (`theta`, in real/unconstrained space) to per-component parameter blocks.
#' This is used throughout the package to:
#' - Define a consistent parameter order for optimization.
#' - Provide an initial parameter vector (`theta0`) in real space.
#' - Record which subset of `theta` belongs to each component model.
#'
#' The function handles two cases:
#' - `time_series_model`: returns a single block of parameters.
#' - `sum_model`: concatenates all component blocks and records their indices.
#'
#' @param model A `time_series_model` or `sum_model`.
#' @return A list with:
#' - `kind`: `"single"` or `"sum"` depending on the model class.
#' - `theta0`: named numeric vector of initial parameters in real space.
#' - `pnames`: (single model) parameter names for the block.
#' - `layout`: (sum model) list of per-component layouts, each containing:
#'   `i` (component index), `idx` (indices into `theta0`), and `pnames`
#'   (parameter names for that component).
#' @keywords internal
prepare_optim_layout <- function(model) {
  get_param_names <- function(m) {
    # Prefer names from the forward transformation (model space -> real space).
    pnames <- names(formals(m$transformation_function))
    if (!is.null(pnames) && length(pnames) > 0L) return(pnames)
    # Fall back to inverse transformation (real space -> model space).
    pnames <- names(formals(m$inv_transformation_function))
    if (!is.null(pnames) && length(pnames) > 0L) return(pnames)
    # Last resort: names on the parameter list itself.
    pnames <- names(m$parameters)
    if (!is.null(pnames) && length(pnames) > 0L) return(pnames)
    NULL
  }

  if (inherits(model, "time_series_model")) {
    # Single model: resolve parameter names, then construct initial theta.
    pnames <- get_param_names(model)
    if (is.null(pnames) || length(pnames) == 0L) {
      stop("Model parameters must be named.", call. = FALSE)
    }

    if (is.null(model$parameters)) {
      # If no parameters are supplied, start from zeros in real space.
      theta0 <- rep(0, length(pnames))
      names(theta0) <- pnames
    } else {
      # Otherwise, map provided parameters back to real space.
      theta0 <- do.call(model$inv_transformation_function, as.list(model$parameters))
      theta0 <- as.numeric(theta0)
      names(theta0) <- names(model$parameters)
    }

    # For single models we return the initial vector and its names.
    return(list(kind = "single", theta0 = theta0, pnames = pnames))
  }

  if (inherits(model, "sum_model")) {
    # Sum model: concatenate each component block and record its indices.
    theta0 <- c()
    layout <- vector("list", length(model$models))
    idx <- 1L

    for (i in seq_along(model$models)) {
      m <- model$models[[i]]
      # Resolve parameter names for this component.
      pnames <- get_param_names(m)
      if (is.null(pnames) || length(pnames) == 0L) {
        stop("Model parameters must be named.", call. = FALSE)
      }
      k <- length(pnames)

      if (is.null(m$parameters)) {
        # If missing, initialize in real space with zeros.
        th_i <- rep(0, k)
      } else {
        # Otherwise, map provided parameters back to real space.
        th_i <- do.call(m$inv_transformation_function, as.list(m$parameters))
        th_i <- as.numeric(th_i)
      }
      # Prefix names to keep component parameters distinct.
      names(th_i) <- paste0(.comp_prefix(i), pnames)

      # Append this block to the flat vector.
      theta0 <- c(theta0, th_i)

      # Record where this block lives in the flat vector.
      layout[[i]] <- list(
        i = i,
        idx = idx:(idx + k - 1L),
        pnames = pnames
      )
      idx <- idx + k
    }

    # For sums we return the full vector and a block layout per component.
    return(list(kind = "sum", theta0 = theta0, layout = layout))
  }

  stop("model must be a 'time_series_model' or 'sum_model'.", call. = FALSE)
}

# mod = wn()+ ar1()
# prepare_optim_layout(mod)
