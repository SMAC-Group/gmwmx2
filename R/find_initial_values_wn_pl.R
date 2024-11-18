#' @importFrom wv wn_to_wv
#' @importFrom stats coef fft lm
find_initial_values_wn_pl <- function(signal, wv_emp) {
  # n=5000
  # sigma2_wn = 10
  # sigma2_pl = 5
  # kappa=-.7
  # vec_autocov_pl = powerlaw_autocovariance(kappa = kappa, sigma2 = sigma2_pl, n = n)
  # plot(vec_autocov_pl)
  # signal = simGauss(autocov = vec_autocov_pl)  +rnorm(n = n, sd =sqrt(sigma2_wn) )
  # wv_emp = wvar(signal)
  # plot(wv_emp)
  # signal = eps_hat_sub
  # wv_emp = wv_emp_eps_hat_filled

  # estimate wn variance based on first scale
  sigma2_wn_start <- 2 * wv_emp$variance[1]

  # Compute power spectrum
  spectrum <- abs(fft(signal))^2
  n <- length(signal)
  freq <- (1:(n / 2)) / n

  # Log-transform frequency and spectrum for analysis
  log_freq <- log(freq)
  log_spectrum <- log(spectrum[1:(n / 2)])

  # plot(log_freq, log_spectrum)

  # Use only the first third of the spectrum estimate to obtain rought estimate of kappa
  end_index <- floor(length(log_freq) * 1 / 3) + 1 # Starting point for the last third
  log_freq_first_third <- log_freq[1:end_index]
  log_spectrum_first_third <- log_spectrum[1:end_index]

  # Estimate kappa using linear regression on the last third of the data
  fit <- lm(log_spectrum_first_third ~ log_freq_first_third)
  kappa_start <- coef(fit)[2]

  # get sigma2 start
  wv_after_removing_wn <- wv_emp$variance - wv::wn_to_wv(sigma2 = sigma2_wn_start, tau = wv_emp$scales)
  n <- tail(wv_emp$scales, 1)
  unit_wv_pl <- autocovariance_to_wv(acf = powerlaw_autocovariance(kappa = kappa_start, sigma2 = 1, n = n), tau = wv_emp$scales)
  sigma2_powerlaw_start <- tail(wv_after_removing_wn, 1) / tail(unit_wv_pl, 1)

  # check for kappa if not greater than -1
  if (kappa_start <= -1) {
    kappa_start <- -0.5
  }

  return(unname(c(abs(sigma2_wn_start), kappa_start, abs(sigma2_powerlaw_start))))
}
