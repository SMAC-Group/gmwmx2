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



# transform from real line to domain
#' Transform power-law kappa from real line to domain
#'
#' Maps unconstrained real values to the stationary domain for the
#' power-law parameter `kappa` using a `tanh` transform.
#'
#' @param x Numeric vector on the real line.
#' @return Numeric vector of transformed values in (-1, 1).
#' @keywords internal
trans_from_real_to_minus_1_and_1 <- function(x) {
  eps = 1e-6
  (1 - eps) * tanh(x)
}

# transform from domain to real line
#' Inverse transform for power-law kappa
#'
#' Maps domain values back to the real line using `atanh`.
#'
#' @param kappa Numeric vector in (-1, 1).
#' @return Numeric vector on the real line.
#' @keywords internal
inv_trans_from_real_to_minus_1_and_1 <- function(kappa) {
  eps = 1e-6
  # safety clamp in case user gives boundary value
  kappa <- pmin(pmax(kappa, -(1 - eps)), (1 - eps))
  atanh(kappa / (1 - eps))
}


# for matern process, alpha should be greater than 1/2 and bounded above
# transform from real line to domain
#' Transform Matérn alpha from real line to domain
#'
#' Ensures Matérn smoothness parameter `alpha` is in (1/2, 10) by mapping
#' real values with a scaled logistic transform.
#'
#' @param x Numeric vector on the real line.
#' @return Numeric vector with values > 1/2.
#' @keywords internal
trans_alpha_matern <- function(x) {
  eps <- 1e-6
  lower <- 1/2 + eps
  upper <- 10 - eps
  lower + (upper - lower) * plogis(x)
}

# transform from domain to real line
#' Inverse transform for Matérn alpha
#'
#' Maps domain values in (1/2, 10) back to the real line.
#'
#' @param x Numeric vector with values > 1/2.
#' @return Numeric vector on the real line.
#' @keywords internal
inv_trans_alpha_matern <- function(x) {
  eps <- 1e-6
  lower <- 1/2 + eps
  upper <- 10 - eps
  qlogis((x - lower) / (upper - lower))
}
