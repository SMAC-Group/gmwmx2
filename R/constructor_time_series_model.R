
# ---------------------------------------------- define constructor for time_series_model

# container
#' Sum of stochastic models (internal)
#'
#' Creates a composite model representing a sum of independent
#' `time_series_model` components. This is an internal constructor used
#' by the `+` methods.
#'
#' @param models A list of `time_series_model` objects.
#' @return A `sum_model` object.
#' @keywords internal
#' @examples
#' mod <- pl(kappa = 0.5, sigma2 = 2) + wn(sigma2 = 1)
#' print(mod)
sum_model <- function(models) {

  # ensure that it is a list of time_series_model objects
  stopifnot(is.list(models))
  if (!all(vapply(models, inherits, logical(1), "time_series_model"))) {
    stop("All elements in `models` must have class 'time_series_model'.", call. = FALSE)
  }

  model_names <- vapply(models, function(m) m$model, character(1))
  if (anyNA(model_names) || any(model_names == "")) {
    stop("Each time_series_model must have a non-empty `$model` name.", call. = FALSE)
  }
  # check for duplicate model names (allow multiple AR(1))
  allow_multiple <- c("AR(1)")
  counts <- table(model_names[!model_names %in% allow_multiple])
  dups <- counts[counts > 1]
  # stop if more than one model of each class provided (except AR(1))
  if (length(dups) > 0) {
    stop(
      "You cannot include the same process more than once (except AR(1)).",
      call. = FALSE
    )
  }

  structure(list(models = models), class = "sum_model")
}



# helper: always returns a LIST of time_series_model objects
#' Coerce to a list of models
#'
#' Helper used by `+` methods to normalize inputs.
#'
#' @param x A `time_series_model` or `sum_model`.
#' @return List of `time_series_model` objects.
#' @keywords internal
as_model_list <- function(x) {
  if (inherits(x, "sum_model")) return(x$models)
  if (inherits(x, "time_series_model")) return(list(x))
  stop("Can only add 'time_series_model' or 'sum_model' objects.")
}

# model + (model or sum_model)
#' Add to a \code{time_series_model} object
#'
#' Combines `time_series_model` and/or `sum_model` into a `sum_model`.
#'
#' @param e1 Left operand.
#' @param e2 Right operand.
#' @return A `sum_model`.
#' @examples
#' m1 <- wn(sigma2 = 1)
#' m2 <- ar1(phi = 0.8, sigma2 = 0.5)
#' model <- m1 + m2
#' model
#' @export
`+.time_series_model` <- function(e1, e2) {
  sum_model(c(as_model_list(e1), as_model_list(e2)))
}

# sum_model + (model or sum_model)
#' Add to a \code{sum_model} object
#'
#' @param e1 Left operand.
#' @param e2 Right operand.
#' @return A `sum_model`.
#' @examples
#' m1 <- wn(sigma2 = 1)
#' m2 <- ar1(phi = 0.8, sigma2 = 0.5)
#' m3 <- pl(kappa = 0.3, sigma2 = 2)
#' model <- (m1 + m2) + m3
#' @export
`+.sum_model` <- function(e1, e2) {
  sum_model(c(as_model_list(e1), as_model_list(e2)))
}

