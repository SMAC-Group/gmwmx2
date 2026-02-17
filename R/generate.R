
# ------------------------------------------------------ function to generate data from stochastic model

#' Generate a time series from a `time_series_model` or `sum_model` object
#'
#' @param object A `time_series_model` or `sum_model`.
#' @param n Length of series to generate.
#' @param seed Optional integer seed for reproducibility.
#' @param ... Passed to method.
#' @return A `generated_time_series` (single model) or
#'   `generated_composite_model_time_series` (sum model).
#' @examples
#' # Single model
#' m1 <- ar1(phi = 0.8, sigma2 = 1)
#' y1 <- generate(m1, n = 200, seed = 123)
#' plot(y1)
#'
#' # Composite model
#' m2 <- wn(sigma2 = 1) + pl(kappa = 0.3, sigma2 = 2)
#' y2 <- generate(m2, n = 200, seed = 123)
#' plot(y2)
#' @export
generate <- function(object, n, seed = NULL, ...) UseMethod("generate")

# Generate from a single model (parameters already in domain)
#' @export
generate.time_series_model <- function(object, n =NULL, seed = NULL, ...) {
  stopifnot(inherits(object, "time_series_model"))
  # test on n
  if(is.null(n)){
    stop("`n` must be a positive integer.")
  }
  if (length(n) != 1L || n <= 0L) stop("`n` must be a positive integer.")
  n <- as.integer(n)


  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      .Random.seed
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        .Random.seed <<- old_seed
      }
    }, add = TRUE)
    set.seed(seed)
  }

  pars <- object$parameters
  x <- do.call(object$generation_function, c(as.list(pars), list(n = n)))

  res <- list(
    series = as.numeric(x),
    n = n,
    model = object$model,
    parameters = pars
  )
  class(res) <- "generated_time_series"
  res
}

# Generate from a missingness model
#' @export
generate.missingness_model <- function(object, n, seed = NULL, ...) {
  stopifnot(inherits(object, "missingness_model"))
  n <- as.integer(n)
  if (length(n) != 1L || n <= 0L) stop("`n` must be a positive integer.")

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      .Random.seed
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        .Random.seed <<- old_seed
      }
    }, add = TRUE)
    set.seed(seed)
  }

  pars <- object$parameters
  x <- do.call(object$generation_function, c(as.list(pars), list(n = n, seed = seed)))

  res <- list(
    series = as.numeric(x),
    n = n,
    model = object$model,
    parameters = pars
  )
  class(res) <- "generated_missingness"
  res
}

# Generate from a sum of models (independent components, reproducible)
#' @export
generate.sum_model <- function(object, n, seed = NULL, ...) {
  stopifnot(inherits(object, "sum_model"))
  n <- as.integer(n)
  if (length(n) != 1L || n <= 0L) stop("`n` must be a positive integer.")

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      .Random.seed
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
          rm(".Random.seed", envir = .GlobalEnv)
        }
      } else {
        .Random.seed <<- old_seed
      }
    }, add = TRUE)
    set.seed(seed)
  }

  xs <- lapply(object$models, generate, n = n, seed = NULL, ...)
  series_list <- lapply(xs, `[[`, "series")
  sum_series <- Reduce(`+`, series_list)

  res <- list(
    series = as.numeric(sum_series),
    components = series_list,
    n = n,
    model = vapply(object$models, `[[`, character(1), "model"),
    parameters = lapply(object$models, `[[`, "parameters")
  )
  class(res) <- "generated_composite_model_time_series"
  res
}
