

# function for matern autocovariance
#' Matern autocovariance
#'
#' Computes the Matern correlation at lags `x` for smoothness `alpha`.
#' This is used internally to build the full autocovariance vector.
#' @importFrom stats runif var
#' @param x Numeric vector of lags (non-negative).
#' @param alpha Smoothness parameter in (1/2, 10).
#' @return Numeric vector of correlations for each lag in `x`.
#' @keywords internal
Ma <- function(x, alpha){
  2/gamma(alpha-1/2)/2^(alpha-1/2)*abs(x)^(alpha-1/2)*besselK(abs(x), abs(alpha-1/2))
}

## -------------------------- parameter checks --------------------------
.check_positive <- function(x, name) {
  if (!is.null(x) && (!is.finite(x) || x <= 0)) {
    stop(name, " must be > 0.", call. = FALSE)
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



# define stationary powerlaw model
#' Stationary Power-Law process
#'
#' Constructs a `time_series_model` representing a stationary power-law
#' process with parameters `kappa` and `sigma2`.
#' In the frequency domain, a power-law process is often described by a
#' spectrum \eqn{P(f) = P_0 f^{\kappa}} (Bos et al., 2008), where \eqn{f} is the frequency, \eqn{P_0} is a constant and \eqn{\kappa} is the spectral index.
#' Note that we use the convention that the power spectral density satisfies
#' \eqn{P(f) \propto |f|^{\kappa}}, where \eqn{\kappa > -1} ensures second-order
#' stationarity. This corresponds to the alternative notation
#' \eqn{P(f) \propto |f|^{-\alpha}} with \eqn{\alpha = -\kappa}.
#' The autocovariance \eqn{\gamma(h) = \mathrm{cov}(X_t, X_{t+h})} used here (Hosking, 1981) is
#' \eqn{\gamma(0) = \sigma^{2} \frac{\Gamma(1+\kappa)}{\Gamma\left(1+\kappa/2\right)^2}},
#' and for \eqn{h > 0}{h > 0}
#' \eqn{\gamma(h) = \frac{-\kappa/2 + h - 1}{\kappa/2 + h}\,\gamma(h-1)}.
#'
#' @param kappa Power-law parameter in (-1, 1). Use `inv_trans_kappa_pl`
#'   for unconstrained optimization.
#' @param sigma2 Process variance (> 0).
#' @return A `time_series_model` object.
#' @examples
#' mod <- pl(kappa = -0.5, sigma2 = 2)
#' mod
#' @export
#' @references
#' Bos MS, Fernandes RMS, Williams SDP, Bastos L (2008). "Fast error analysis
#' of continuous GPS observations." *Journal of Geodesy*, 82, 157-166.
#'
#' Hosking JRM (1981). "Fractional differencing." *Biometrika*, 68(1), 165-176.
pl = function(kappa = NULL, sigma2 = NULL){
  .check_in_open_interval(kappa, "kappa", -1, 1)
  .check_positive(sigma2, "sigma2")
  res = list(
    "parameters" = c("kappa" = kappa, "sigma2" = sigma2),
    "model" = "Stationary PowerLaw",
    "transformation_function" = function(kappa, sigma2){
      # transform parameter from real to domain parameters
      return(c(trans_from_real_to_minus_1_and_1(kappa), exp(sigma2)))
    },
    "inv_transformation_function" = function(kappa, sigma2){
      # transform parameter from domain to real parameters
    return(c(inv_trans_from_real_to_minus_1_and_1(kappa), log(sigma2)))
    },
    "autocovariance_function" = function(kappa, sigma2, n){
      autocov = powerlaw_autocovariance(kappa, sigma2, n)
      return(autocov)
    },
    "generation_function" =  function(kappa, sigma2, n){
      autocov_vec = powerlaw_autocovariance(kappa, sigma2, n)
      x = longmemo::simGauss(autocov_vec)
      return(x)
    },
    "get_initial_parameters_function" = function(signal){
      kappa = runif(n = 1, min = -.8, max =  .8)
      sigma2 = var(signal)
      return(c("kappa" = kappa, "sigma2"=sigma2))
    },
    "get_variance_covariance_matrix_signal" = function(kappa, sigma2, n){
      autocov_vec = powerlaw_autocovariance(kappa, sigma2, n)
      cov_mat = fast_toeplitz_matrix_from_vector_cpp(autocov_vec)
      return(cov_mat)
    }

  )

  class(res) = "time_series_model"
  return(res)
}


# function to define a white noise
#' White noise process
#'
#' Constructs a `time_series_model` for white noise with variance `sigma2`.
#' The process is defined as
#' \eqn{X_t \stackrel{\text{i.i.d.}}{\sim} N(0, \sigma^2)}{X_t ~ N(0, sigma^2), i.i.d.}
#' with autocovariance
#' \eqn{\gamma(h) = \mathrm{cov}(X_t, X_{t+h}) = \sigma^2 \mathbf{1}\{h=0\}}{gamma(h) = cov(X_t, X_{t+h}) = sigma^2 * 1{h=0}}
#' @importFrom stats rnorm
#' @param sigma2 Innovation variance (> 0).
#' @return A `time_series_model` object.
#' @examples
#' mod <- wn(sigma2 = 1)
#' mod
#' @export
wn = function(sigma2 = NULL){
  .check_positive(sigma2, "sigma2")
  res = list(
    "parameters" = c("sigma2" = sigma2),
    "model" = "White Noise",
    "transformation_function" = function(sigma2){
      # transform parameter from real to domain parameters
      return(c(exp(sigma2)))
    },
    "inv_transformation_function" = function(sigma2){
      # transform parameter from domain to real parameters
      return(c(log(sigma2)))
    },
    "autocovariance_function" = function(sigma2, n){
      autocov = rep(0, n)
      autocov[1] = sigma2
      return(autocov)
    },
    "generation_function" =  function(sigma2, n){
      x <- rnorm(n = n, mean = 0, sd = sqrt(sigma2))
      return(x)
    },
    "get_initial_parameters_function" = function(signal){
      sigma2 = var(signal)
      return(c("sigma2"=sigma2))
    },
    "get_variance_covariance_matrix_signal" = function(sigma2, n){

      cov_mat = diag(n) * sigma2
      return(cov_mat)
    }

  )

  class(res) = "time_series_model"
  return(res)
}


# define matern model
#' Matern process
#'
#' Constructs a `time_series_model` for a Matern covariance process with
#' variance `sigma2`, range `lambda`, and smoothness `alpha`.
#' The autocovariance is
#' \eqn{ \gamma(h) = \mathrm{cov}(X_t, X_{t+h}) = \frac{2 \sigma^2}{\Gamma(\alpha-1 / 2) 2^{\alpha-1 / 2}}|\lambda h|^{\alpha-1 / 2} \mathcal{K}_{|\alpha-1 / 2|}(| \lambda h|)}
#' where \eqn{\mathcal{K}_\omega(x)} is the modified Bessel function of the second kind of order \eqn{\omega}.
#'
#' @param sigma2 Marginal variance (> 0).
#' @param lambda Range/scale parameter (> 0).
#' @param alpha Smoothness parameter in (1/2, 10).
#' @return A `time_series_model` object.
#' @examples
#' mod <- matern(sigma2 = 1, lambda = 0.2, alpha = 1.0)
#' mod
#' @export
#' @references
#' Lilly JM, Sykulski AM, Early JJ, Olhede SC (2017). "Fractional Brownian motion,
#' the Matérn process, and stochastic modeling of turbulent dispersion."
#' *Nonlinear Processes in Geophysics*, 24(3), 481-514.
matern = function(sigma2=NULL, lambda=NULL, alpha=NULL){
  .check_positive(sigma2, "sigma2")
  .check_positive(lambda, "lambda")
  .check_in_open_interval(alpha, "alpha", 1/2, 10)
  res = list(
    "parameters" = c("sigma2" = sigma2, "lambda" = lambda, "alpha" = alpha),
  "model" = "Matern",
    "transformation_function" = function(sigma2, lambda, alpha){
      return(c(exp(sigma2), exp(lambda), trans_alpha_matern(alpha) ) )
    },
    "inv_transformation_function" = function(sigma2, lambda, alpha){
      # transform parameter from domain to real parameters
      return(c(log(sigma2), log(lambda), inv_trans_alpha_matern(alpha) ) )
    },
    "autocovariance_function" = function(sigma2, lambda, alpha, n) {
      autocov = c(sigma2, sigma2*Ma(lambda* (1:(n-1)), alpha = alpha))
      return(autocov)
    },
    "generation_function" =  function(sigma2, lambda, alpha, n) {
      autocov_vec = c(sigma2, sigma2*Ma(lambda* (1:(n-1)), alpha = alpha))
      x = longmemo::simGauss(autocov_vec)
      return(x)
    },
  "get_initial_parameters_function" = function(signal){
    sigma2 = var(signal)
    lambda = 1
    alpha = 1
    return(c("sigma2"=sigma2, "lambda"=lambda, "alpha"=alpha))
  },
  "get_variance_covariance_matrix_signal" = function(sigma2, lambda, alpha, n){
    autocov_vec = c(sigma2, sigma2*Ma(lambda* (1:(n-1)), alpha = alpha))
    cov_mat = fast_toeplitz_matrix_from_vector_cpp(autocov_vec)
    return(cov_mat)
  }

  )

  class(res) = "time_series_model"
  return(res)
}

# define ar1
#' AR(1) process
#'
#' Constructs a `time_series_model` for a stationary AR(1) process with parameter
#' `phi` and innovation variance `sigma2`.
#' The model is
#' \eqn{X_t = \phi X_{t-1} + \varepsilon_t, \; \varepsilon_t \stackrel{\text{i.i.d.}}{\sim} N(0, \sigma^2)}{X_t = phi X_{t-1} + eps_t, eps_t ~ N(0, sigma^2)}.
#' The autocovariance is
#' \eqn{
#' \gamma(h) = \mathrm{cov}(X_t, X_{t+h})
#' = \frac{\sigma^2}{1 - \phi^2}\,\phi^{\lvert h \rvert}
#' }{
#' gamma(h) = cov(X_t, X_{t+h})
#' = sigma^2/(1 - phi^2) * phi^{|h|}
#' }.
#' @param phi AR(1) coefficient in (-1, 1).
#' @param sigma2 Innovation variance (> 0).
#' @return A `time_series_model` object.
#' @examples
#' mod <- ar1(phi = 0.8, sigma2 = 1)
#' mod
#' @export
ar1 <- function(phi = NULL, sigma2 = NULL) {
  .check_in_open_interval(phi, "phi", -1, 1)
  .check_positive(sigma2, "sigma2")
  res <- list(
    "parameters" = c("phi" = phi, "sigma2" = sigma2),
    "model" = "AR(1)",

    # kept for consistency (not used if you keep params in domain)
    "transformation_function" = function(phi, sigma2) c(trans_from_real_to_minus_1_and_1(phi), exp(sigma2)),
    "inv_transformation_function" = function(phi, sigma2) c(inv_trans_from_real_to_minus_1_and_1(phi), log(sigma2)),

    "autocovariance_function" = function(phi, sigma2, n) {
      # if (!is.finite(phi) || abs(phi) >= 1) stop("phi must be in (-1, 1) for stationarity.")
      # if (!is.finite(sigma2) || sigma2 <= 0) stop("sigma2 must be > 0 (innovation variance).")
      n <- as.integer(n)
      if (length(n) != 1L || is.na(n) || n <= 0L) stop("`n` must be a positive integer.")

      sigma2 * (phi^(0:(n - 1))) / (1 - phi^2)
    },

    "generation_function" = function(phi, sigma2, n) {
      if (!is.finite(phi) || abs(phi) >= 1) stop("phi must be in (-1, 1) for stationarity.")
      if (!is.finite(sigma2) || sigma2 <= 0) stop("sigma2 must be > 0 (innovation variance).")
      n <- as.integer(n)
      if (length(n) != 1L || is.na(n) || n <= 0L) stop("`n` must be a positive integer.")
      burnin = 200L
      x <- stats::arima.sim(
        model = list(ar = phi),
        n = n,
        sd = sqrt(sigma2),
        n.start = as.integer(burnin)
      )
      as.numeric(x)
    },
    "get_initial_parameters_function" = function(signal){
      phi = runif(n = 1, min = -.5, max = .5)
      sigma2 = var(signal)
      return(c("phi"=phi, "sigma2"=sigma2))
    },
    "get_variance_covariance_matrix_signal" = function(phi, sigma2, n){
      autocov_vec = sigma2 * (phi^(0:(n - 1))) / (1 - phi^2)
      cov_mat = fast_toeplitz_matrix_from_vector_cpp(autocov_vec)
      return(cov_mat)
    }
  )

  class(res) <- "time_series_model"
  res
}


#' Random walk process
#'
#' Constructs a `time_series_model` for a random walk with innovation
#' variance `sigma2`. The autocovariance returned is the mean of the
#' diagonal and super-diagonals of the covariance matrix.
#' The model is
#' \eqn{X_t = X_{t-1} + \varepsilon_t, \; \varepsilon_t \stackrel{\text{i.i.d.}}{\sim} N(0, \sigma^2)}{X_t = X_{t-1} + eps_t, eps_t ~ N(0, sigma^2)}.
#'
#' @param sigma2 Innovation variance (> 0).
#' @return A `time_series_model` object.
#' @examples
#' mod <- rw(sigma2 = 1)
#' mod
#' @export
rw = function(sigma2 =NULL){
  .check_positive(sigma2, "sigma2")
  res = list(
    "parameters" = c("sigma2" = sigma2),
    "model" = "Random Walk",
    "transformation_function" = function(sigma2){
      return(c(exp(sigma2)) )
    },
    "inv_transformation_function" = function(sigma2){
      # transform parameter from domain to real parameters
      return(c(log(sigma2) ) )
    },
    # note that this is not the autocovariance function of the RW process, but the average of the diagonal and super diagonals of the variance covariance matrix of a RW process
    "autocovariance_function" = function(sigma2, n) {
      autocov = get_mean_diagonal_super_diagonals_cov_mat_rw_cpp(n, sigma2)
      return(autocov)
    },
    "generation_function" =  function(sigma2, n) {
      x = cumsum(rnorm(n, mean=0, sd=sqrt(sigma2)))
      return(x)
    },
    "get_initial_parameters_function" = function(signal){
      sigma2 = var(signal)
      return(c("sigma2"=sigma2))
    },
    "get_variance_covariance_matrix_signal" = function(sigma2, n){
      cov_mat = get_sigma_mat_rw(n, sigma2)
      return(cov_mat)
    }

  )

  class(res) = "time_series_model"
  return(res)
}



# define flicker noise model
#' Flicker noise process
#'
#' Constructs a `time_series_model` for flicker noise with
#' variance `sigma2`.
#' The process has spectral density
#' \eqn{S(f) \propto \frac{1}{|f|}}{S(f) proportional to 1/|f|}. Hence, \eqn{\kappa = -1}{kappa = -1} (Bos et al., 2008).
#' The process is non-stationary and its covariance matrix is assumed to be
#' given by
#' \deqn{
#' \mathbf C = \sigma^2 \mathbf U^\top \mathbf U,
#' }
#' where \eqn{\mathbf U \in \mathbb{R}^{N \times N}} is an upper-triangular
#' Toeplitz matrix with entries
#' \deqn{
#' U_{i,j} =
#' \begin{cases}
#' h_{j-i}, & j \ge i, \\
#' 0,       & j < i,
#' \end{cases}
#' \qquad i,j = 1, \ldots, N.
#' }
#' The coefficients \eqn{\{h_i\}_{i \ge 0}} define a causal linear filter and
#' are given recursively by
#' \deqn{
#' h_0 = 1, \qquad
#' h_i = \left(i - \frac{\kappa}{2} - 1\right)\frac{h_{i-1}}{i},
#' \quad i > 0.
#' }
#' @param sigma2 Innovation variance (> 0).
#' @return A `time_series_model` object.
#' @examples
#' mod <- flicker(sigma2 = 1)
#' mod
#' @export
#' @references
#' Bos MS, Fernandes RMS, Williams SDP, Bastos L (2008). "Fast error analysis
#' of continuous GPS observations." *Journal of Geodesy*, 82, 157-166.
flicker = function(sigma2 = NULL){
  .check_positive(sigma2, "sigma2")
  res = list(
    "parameters" = c("sigma2" = sigma2),
    "model" = "Flicker",
    "transformation_function" = function(sigma2){
      # transform parameter from real to domain parameters
      return(c(exp(sigma2)))
    },
    "inv_transformation_function" = function(sigma2){
      # transform parameter from domain to real parameters
      return(c(log(sigma2)))
    },
    # note that this is not the autocovariance function of the FL process, but the average of the diagonal and super diagonals of the variance covariance matrix of a FL process
    "autocovariance_function" = function(sigma2, n){
      autocov = vec_mean_autocov_powerlaw(-1, n) * sigma2
      return(autocov)
    },
    "generation_function" =  function(sigma2, n){
      x = gen_flicker(n, sqrt(sigma2))
      return(x)
    },
    "get_initial_parameters_function" = function(signal){
      sigma2 = var(signal)
      return(c("sigma2"=sigma2))
    },
    "get_variance_covariance_matrix_signal" = function(sigma2, n){
      cov_mat = var_cov_powerlaw_cpp(sigma2, -1, n)
      return(cov_mat)
    }

  )

  class(res) = "time_series_model"
  return(res)
}






#' Markov two-state missingness model
#'
#' Constructs a `missingness_model` representing a two-state Markov process
#' for missing/observed indicators. The process takes values in \{0, 1\},
#' where 1 indicates observed and 0 indicates missing.
#'
#' @param p1 Transition probability from observed (1) to missing (0).
#' @param p2 Transition probability from missing (0) to observed (1).
#' @return A `missingness_model` object.
#' @examples
#' mod <- markov_two_states(p1 = 0.05, p2 = 0.95)
#' mod
#' z <- generate(mod, n = 200, seed = 123)
#' plot(z)
#' @importFrom stats rbinom
#' @export
markov_two_states <- function(p1 = NULL, p2 = NULL) {

  res <- list(
    "parameters" = c("p1" = p1, "p2" = p2),
    "model" = "Markov Two States",

    "generation_function" = function(p1, p2, n, seed = NULL) {
        # create vector
        c1 <- 1000
        vec_omega <- vector(mode = "numeric", length = n + c1)
        vec_omega[1] <- 1
        if (!is.null(seed)) {
          set.seed(seed)
        }
        for (i in 2:(n + c1)) {
          # generate
          if (vec_omega[i - 1] == 1) {
            vec_omega[i] <- rbinom(n = 1, size = 1, prob = (1 - p1))
          } else if (vec_omega[i - 1] == 0) {
            vec_omega[i] <- rbinom(n = 1, size = 1, prob = p2)
          }
        }
        return(tail(vec_omega, n = n))
    }
  )
  class(res) <- "missingness_model"
  return(res)
}




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
#' Add two model objects
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
#' Add to a sum_model
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



# --------------- print method for time_series_model and sum_model

# Helper: format parameters nicely
#' Format parameter display
#'
#' @param pars Named numeric vector.
#' @param digits Significant digits.
#' @return Single string for printing.
#' @keywords internal
format_params <- function(pars, digits = 4) {
  if (is.null(pars) || length(pars) == 0) return("(no parameters)")
  paste0(
    names(pars), " = ",
    formatC(pars, digits = digits, format = "g"),
    collapse = ", "
  )
}

# Print a single stochastic process
#' @export
print.time_series_model <- function(x, ...) {
  cat("Stochastic process\n")
  cat("  Model      :", x$model, "\n")
  cat("  Parameters :", format_params(x$parameters), "\n")
  invisible(x)
}

# Print a sum of stochastic processes
#' @export
print.sum_model <- function(x, ...) {
  n <- length(x$models)
  cat("Stochastic model: sum of", n, "processes\n\n")

  for (i in seq_along(x$models)) {
    m <- x$models[[i]]
    cat(sprintf(" [%d] %s\n", i, m$model))
    cat("     Parameters :", format_params(m$parameters), "\n\n")
  }

  invisible(x)
}

# Print a missingness process
#' @export
print.missingness_model <- function(x, ...) {
  cat("Missingness process\n")
  cat("  Model      :", x$model, "\n")
  cat("  Parameters :", format_params(x$parameters), "\n")
  invisible(x)
}




# ------------------------------------------------------ example of usage


# mod = pl(kappa=0.5, sigma2=2) + wn(sigma2=1)
# print(mod)



# ------------------------------------------------------ function to generate data from stochastic model

#' Generate a time series from a model
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


# --------------- plot methods for generated series

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



#' Plot a generated time series
#'
#' Produces a single line plot for a `generated_time_series` object.
#'
#' @param x A `generated_time_series`.
#' @param ... Additional arguments passed to `plot()`.
#' @return Invisibly returns `x`.
#' @examples
#' m1 <- wn(sigma2 = 1)
#' y1 <- generate(m1, n = 200, seed = 123)
#' plot(y1)
#' @export
plot.generated_time_series <- function(x, ...) {
  y <- x$series
  n <- length(y)
  line_col <- .gmwmx2_get_plot_colors(1L)
  plot(
    seq_len(n), y, type = "n",
    xlab = "Time", ylab = "Value",
    xlim = c(1, n),
    main = "",
    las = 1,
    ...
  )
  graphics::mtext(side = 3, text = paste0("Model: ", x$model), line = 0.2)
  graphics::grid(col = "grey85", lty=1)
  lines(seq_len(n), y, lty = 1, col = line_col)
  invisible(x)
}

#' Plot a generated missingness process
#'
#' Produces a step plot for a `generated_missingness` object.
#'
#' @param x A `generated_missingness`.
#' @param ... Additional arguments passed to `plot()`.
#' @return Invisibly returns `x`.
#' @examples
#' m0 <- markov_two_states(p1 = 0.05, p2 = 0.9)
#' z0 <- generate(m0, n = 200, seed = 123)
#' plot(z0)
#' @export
plot.generated_missingness <- function(x, ...) {
  y <- x$series
  n <- length(y)
  line_col <- .gmwmx2_get_plot_colors(1L)
  plot(
    seq_len(n), y, type = "n",
    xlab = "Time", ylab = "Observed (1) / Missing (0)",
    xlim = c(1, n),
    ylim = c(-0.05, 1.05),
    main = "",
    las = 1,
    ...
  )

  graphics::mtext(side = 3, text = paste0("Missingness model: ", x$model), line = 0.2)
  graphics::grid(col = "grey85", lty=1)
  lines(seq_len(n), y, type = "s", lty = 1, col = line_col)
  invisible(x)
}

#' Plot a generated composite time series
#'
#' Produces stacked line plots for each component and the sum for a
#' `generated_composite_model_time_series` object.
#'
#' @param x A `generated_composite_model_time_series`.
#' @param ... Additional arguments passed to `plot()`.
#' @return Invisibly returns `x`.
#' @examples
#' m2 <- wn(sigma2 = 1) + ar1(phi = 0.8, sigma2 = 0.5)
#' y2 <- generate(m2, n = 200, seed = 123)
#' plot(y2)
#' @export
plot.generated_composite_model_time_series <- function(x, ...) {
  comps <- x$components
  k <- length(comps)
  n <- x$n
  if (k < 1L) stop("No components to plot.", call. = FALSE)
  line_cols <- .gmwmx2_get_plot_colors(k)

  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)
  graphics::par(mfrow = c(k + 1L, 1L), mar = c(3, 4, 2, 1))

  for (i in seq_len(k)) {
    series_i <- comps[[i]]
    main_i <- x$model[[i]]
    plot(
      seq_len(n), series_i, type = "n",
      xlab = "Time", ylab = "Value",
      xlim = c(1, n),
      main = "",
      las = 1,
      ...
    )
    graphics::mtext(side = 3, text = main_i, line = 0.2)
    graphics::grid(col = "grey85", lty=1)
    lines(seq_len(n), series_i, lty = 1, col = line_cols[[i]])
  }

  sum_title <- paste0("Sum of: ", paste(x$model, collapse = " + "))
  plot(
    seq_len(n), x$series, type = "n",
    xlab = "Time", ylab = "Value",
    xlim = c(1, n),
    main = "",
    las = 1,
    ...
  )
  graphics::mtext(side = 3, text = sum_title, line = 0.2)
  graphics::grid(col = "grey85", lty=1)
  lines(seq_len(n), x$series, lty = 1, col = "grey40")

  invisible(x)
}




# ------------------------------------ exemple of usage
# model1 =ar1(phi=.91, sigma2=10) + wn(10)
# print(model1)
# x1 = generate(model1, n = 10000, seed=123)
# str(x1)
# plot(x1)
