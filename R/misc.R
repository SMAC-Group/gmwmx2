#' Plot colors for generated time series
#'
#' Editable default palette used by the plotting methods for generated time
#' series. When more colors are needed, the palette is expanded with
#' interpolation.
#'
#' @format A character vector of hex colors.
#' @keywords internal
#' @seealso .gmwmx2_get_plot_colors
gmwmx2_plot_colors <- c(
  "#8DA0CB",  # muted lavender-blue
  "#FC8D62",  # soft coral
  "#66C2A5",  # pastel teal
  "#E78AC3",  # dusty pink
  "#A6D854",  # gentle green
  "#FFD92F"   # soft mustard yellow
)

#' Get plot colors for generated time series
#'
#' Returns a vector of plot colors of length `n`. If `n` exceeds the size of
#' `gmwmx2_plot_colors`, a larger palette is created via interpolation.
#'
#' @param n Number of colors to return.
#' @return Character vector of hex colors.
#' @keywords internal
#' @seealso gmwmx2_plot_colors
.gmwmx2_get_plot_colors <- function(n) {
  if (n <= 0L) return(character())
  if (n <= length(gmwmx2_plot_colors)) {
    return(gmwmx2_plot_colors[seq_len(n)])
  }
  grDevices::colorRampPalette(gmwmx2_plot_colors)(n)
}



# Prefix for component i in a sum_model
#' Component name prefix for summed models
#'
#' @param i Component index.
#' @return Prefix string like `m1_`.
#' @keywords internal
.comp_prefix <- function(i) paste0("m", i, "_")



#' Format model text
#'
#' Internal helper to create a human-readable model label for printing.
#'
#' @param model A \code{time_series_model} or \code{sum_model}.
#' @return A single string describing the model.
#' @keywords internal
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

.check_in_open_interval <- function(x, name, lower, upper) {
  if (!is.null(x) && (!is.finite(x) || x <= lower || x >= upper)) {
    stop(name, " must be in (", lower, ", ", upper, ").", call. = FALSE)
  }
}

.check_greater_than <- function(x, name, lower) {
  if (!is.null(x) && (!is.finite(x) || x <= lower)) {
    stop(name, " must be > ", lower, ".", call. = FALSE)
  }
}


