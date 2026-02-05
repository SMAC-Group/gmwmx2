# define transformation functions for powerlaw

# transform from real line to domain
#' Transform power-law kappa from real line to domain
#'
#' Maps unconstrained real values to the stationary domain for the
#' power-law parameter `kappa` using a `tanh` transform.
#'
#' @param x Numeric vector on the real line.
#' @return Numeric vector of transformed values in (-1, 1).
#' @keywords internal
trans_kappa_pl <- function(x) {
  # exp(x) - 1
  tanh(x)
}

# transform from domain to real line
#' Inverse transform for power-law kappa
#'
#' Maps domain values back to the real line using `atanh`.
#'
#' @param x Numeric vector in (-1, 1).
#' @return Numeric vector on the real line.
#' @keywords internal
inv_trans_kappa_pl <- function(x){
  # log(x + 1)
  atanh(x)
}


# for matern process, alpha should be greater than 1/2
# transform from real line to domain
#' Transform Matérn alpha from real line to domain
#'
#' Ensures Matérn smoothness parameter `alpha` is > 1/2 by mapping
#' real values with `exp(x) + 1/2`.
#'
#' @param x Numeric vector on the real line.
#' @return Numeric vector with values > 1/2.
#' @keywords internal
trans_alpha_matern <- function(x) {
  exp(x) + 1/2
}

# transform from domain to real line
#' Inverse transform for Matérn alpha
#'
#' Maps domain values (> 1/2) back to the real line.
#'
#' @param x Numeric vector with values > 1/2.
#' @return Numeric vector on the real line.
#' @keywords internal
inv_trans_alpha_matern <- function(x) log(x - 1/2)
