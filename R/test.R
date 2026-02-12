#
# library(gmwmx2)
# library(wv)
# library(dplyr)
#
#
# n = 10*365
# X = matrix(NA, nrow = n, ncol = 4)
# # intercept
# X[, 1] = 1
# # trend
# X[, 2] = 1:n
# # annual sinusoid
# omega_1 <- (1 / 365.25) * 2 * pi
# X[, 3] <- sin((1:n) * omega_1)
# X[, 4] <- cos((1:n) * omega_1)
# beta = c(1, 0.01, 2, 1.5)
#
# # visualize the deterministic signal
# plot(x = X[, 2], y = X %*% beta, type = "l",
#      main = "Signal", xlab = "Time", ylab = "Signal")
#
#
# # Example 1: White noise + AR(1)
#
# phi_ar1 = 0.95
# sigma2_ar1 = 20
# sigma2_wn = 20
# eps = gmwmx2::generate(ar1(phi = phi_ar1, sigma2 = sigma2_ar1) + wn(sigma2_wn), n = n, seed = 123)$series
# plot(wv::wvar(eps))
# y = X %*% beta + eps
# # generate missingness
# mod_missing = markov_two_states(p1 = .05, p2 = .45)
# # generate missing
# z = generate(mod_missing, n = n, seed = 123)
# plot(z)
# mean(z$series==1)
# # mask entries in y
# y[z$series == 0] = NA
# plot(X[, 2], y, type = "l", main = "Simulated y (WN + AR1)")
#
#
# # Fit model (WN + AR1)
# fit = gmwmx2_new(X = X, y = y, model = wn() + ar1() )
# print(fit)
#
#
#
#
# B = 100
# mat_res = matrix(NA, nrow=B, ncol=19)
# for(b in seq(B)){
#   # b=1
#   eps = generate(ar1(phi=0.95, sigma2=20) + wn(20), n=n, seed = (123 + b))$series
#   y = X %*% beta + eps
#   z = generate(mod_missing, n = n, seed = 123+b)
#   # mask entries in y
#   y[z$series == 0] = NA
#   fit = gmwmx2_new(X = X, y = y, model = wn() + ar1() )
#   # mispecified model assuming white noise as the stochastic model only fitted on data we have
#   X_sub = X[z$series == 1, ]
#   y_sub = y[z$series == 1]
#   fit2 = lm(y_sub~X_sub[,2] + X_sub[,3] + X_sub[,4])
#
#   mat_res[b, ] = c(fit$beta_hat, fit$std_beta_hat,
#                    summary(fit2)$coefficients[,1],
#                    summary(fit2)$coefficients[,2],
#                    fit$theta_domain$`AR(1)_2`,
#                    fit$theta_domain$`White Noise_1`)
#   cat("Iteration ", b, " completed.\n")
# }
#
#
# # compute empirical coverage
# mat_res_df = as.data.frame(mat_res)
# colnames(mat_res_df) = c("gmwmx_beta0_hat", "gmwmx_beta1_hat",
#                          "gmwmx_beta2_hat", "gmwmx_beta3_hat",
#                          "gmwmx_std_beta0_hat", "gmwmx_std_beta1_hat",
#                          "gmwmx_std_beta2_hat", "gmwmx_std_beta3_hat",
#                          "lm_beta0_hat", "lm_beta1_hat", "lm_beta2_hat", "lm_beta3_hat",
#                          "lm_std_beta0_hat", "lm_std_beta1_hat", "lm_std_beta2_hat", "lm_std_beta3_hat",
#                          "phi_ar1","sigma_2_ar1" ,"sigma_2_wn")
#
#
# # Plot empirical distributions (beta and stochastic parameters)
#
#
# par(mfrow = c(1, 4))
# boxplot(mat_res_df[, c("gmwmx_beta0_hat")],las=1,
#         names = c("beta0"),
#         main = expression(beta[0]), ylab = "Estimate")
# abline(h = beta[1], col = "black", lwd = 2)
#
# boxplot(mat_res_df[, c("gmwmx_beta1_hat")],las=1,
#         names = c("beta1"),
#         main = expression(beta[1]), ylab = "Estimate")
# abline(h = beta[2], col = "black", lwd = 2)
#
# boxplot(mat_res_df[, c("gmwmx_beta2_hat")],las=1,
#         names = c("beta2"),
#         main = expression(beta[2]), ylab = "Estimate")
# abline(h = beta[3], col = "black", lwd = 2)
#
# boxplot(mat_res_df[, c("gmwmx_beta3_hat")],las=1,
#         names = c("beta3"),
#         main = expression(beta[3]), ylab = "Estimate")
# abline(h = beta[4], col = "black", lwd = 2)
#
#
# par(mfrow = c(1, 3))
#
# boxplot(mat_res_df$phi_ar1,las=1,
#         names = c("phi_ar1"),
#         main = expression(phi["AR1"]), ylab = "Estimate")
# abline(h = phi_ar1, col = "black", lwd = 2)
#
# boxplot(mat_res_df$sigma_2_ar1,las=1,
#         names = c("phi_ar1"),
#         main = expression(sigma["AR1"]^2), ylab = "Estimate")
# abline(h = sigma2_ar1, col = "black", lwd = 2)
#
# boxplot(mat_res_df$sigma_2_wn,las=1,
#         names = c("phi_ar1"),
#         main = expression(sigma["WN"]^2), ylab = "Estimate")
# abline(h = sigma2_wn, col = "black", lwd = 2)
# par(mfrow = c(1, 1))
#
#
# # Compute empirical coverage of confidence intervals for beta
#
# zval = qnorm(0.975)
# mat_res_df$upper_ci_gmwmx_beta0 = mat_res_df$gmwmx_beta0_hat + zval * mat_res_df$gmwmx_std_beta0_hat
# mat_res_df$lower_ci_gmwmx_beta0 = mat_res_df$gmwmx_beta0_hat - zval * mat_res_df$gmwmx_std_beta0_hat
# mat_res_df$upper_ci_gmwmx_beta1 = mat_res_df$gmwmx_beta1_hat + zval * mat_res_df$gmwmx_std_beta1_hat
# mat_res_df$lower_ci_gmwmx_beta1 = mat_res_df$gmwmx_beta1_hat - zval * mat_res_df$gmwmx_std_beta1_hat
# # empirical coverage of gmwmx beta
# dplyr::between(rep(beta[1], 500), mat_res_df$lower_ci_gmwmx_beta0, mat_res_df$upper_ci_gmwmx_beta0) %>% mean()
# dplyr::between(rep(beta[2], 500), mat_res_df$lower_ci_gmwmx_beta1, mat_res_df$upper_ci_gmwmx_beta1) %>% mean()
#
# # do the same for lm beta
# mat_res_df$upper_ci_lm_beta0 = mat_res_df$lm_beta0_hat + zval * mat_res_df$lm_std_beta0_hat
# mat_res_df$lower_ci_lm_beta0 = mat_res_df$lm_beta0_hat - zval * mat_res_df$lm_std_beta0_hat
# mat_res_df$upper_ci_lm_beta1 = mat_res_df$lm_beta1_hat + zval * mat_res_df$lm_std_beta1_hat
# mat_res_df$lower_ci_lm_beta1 = mat_res_df$lm_beta1_hat - zval * mat_res_df$lm_std_beta1_hat
# dplyr::between(rep(beta[1], 500), mat_res_df$lower_ci_lm_beta0, mat_res_df$upper_ci_lm_beta0) %>% mean()
# dplyr::between(rep(beta[2], 500), mat_res_df$lower_ci_lm_beta1, mat_res_df$upper_ci_lm_beta1) %>% mean()
#
#
#
#
# # Example 2: White noise + stationary power-law
#
# kappa_pl <- -0.8
# sigma2_wn <- 10
# sigma2_pl <- 2
# eps_pl = generate(wn(sigma2_wn) + pl(kappa = kappa_pl, sigma2 = sigma2_pl),
#                   n = n, seed = 123)
# plot(eps_pl$series, type = "l")
# plot(wv::wvar(eps_pl$series))
# y_pl = X %*% beta + eps_pl$series
# # linear model assuming white noise for stochastic model and only fitted on data we have
# # generate missing
# z = generate(mod_missing, n = n, seed = 123)
# plot(z)
# mean(z$series==1)
# # mask entries in y
# y_pl[z$series == 0] = NA
# fit_pl = gmwmx2_new(X = X, y = y_pl, model = wn() + pl())
# fit_pl
#
#
# B_pl = 100
# mat_res_pl = matrix(NA, nrow = B_pl, ncol = 19)
# for(b in seq(B_pl)){
#   eps = generate(wn(sigma2_wn) + pl(kappa = kappa_pl, sigma2 = sigma2_pl),
#                  n = n, seed = (123 + b))$series
#   y = X %*% beta + eps
#   z = generate(mod_missing, n = n, seed = 123+b)
#   # mask entries in y
#   y[z$series == 0] = NA
#   fit = gmwmx2_new(X = X, y = y, model = wn() + pl())
#
#   # linear model assuming white noise for stochastic model and only fitted on data we have
#   X_sub = X[z$series == 1, ]
#   y_sub = y[z$series == 1]
#   fit2 = lm(y_sub~X_sub[,2] + X_sub[,3] + X_sub[,4])
#
#   mat_res_pl[b, ] = c(fit$beta_hat, fit$std_beta_hat,
#                       summary(fit2)$coefficients[,1],
#                       summary(fit2)$coefficients[,2],
#                       fit$theta_domain$`Stationary PowerLaw_2`,
#                       fit$theta_domain$`White Noise_1`)
#   cat("Iteration ", b, " \n")
# }
#
# mat_res_pl_df = as.data.frame(mat_res_pl)
# colnames(mat_res_pl_df) = c("gmwmx_beta0_hat", "gmwmx_beta1_hat",
#                             "gmwmx_beta2_hat", "gmwmx_beta3_hat",
#                             "gmwmx_std_beta0_hat", "gmwmx_std_beta1_hat",
#                             "gmwmx_std_beta2_hat", "gmwmx_std_beta3_hat",
#                             "lm_beta0_hat", "lm_beta1_hat", "lm_beta2_hat", "lm_beta3_hat",
#                             "lm_std_beta0_hat", "lm_std_beta1_hat",
#                             "lm_std_beta2_hat", "lm_std_beta3_hat",
#                             "kappa_pl","sigma2_pl" ,"sigma2_wn")
#
#
# # Plot empirical distributions (beta and stochastic parameters)
#
#
# par(mfrow = c(1, 4))
# boxplot(mat_res_df[, c("gmwmx_beta0_hat")],las=1,
#         names = c("beta0"),
#         main = expression(beta[0]), ylab = "Estimate")
# abline(h = beta[1], col = "black", lwd = 2)
#
# boxplot(mat_res_df[, c("gmwmx_beta1_hat")],las=1,
#         names = c("beta1"),
#         main = expression(beta[1]), ylab = "Estimate")
# abline(h = beta[2], col = "black", lwd = 2)
#
# boxplot(mat_res_df[, c("gmwmx_beta2_hat")],las=1,
#         names = c("beta2"),
#         main = expression(beta[2]), ylab = "Estimate")
# abline(h = beta[3], col = "black", lwd = 2)
#
# boxplot(mat_res_df[, c("gmwmx_beta3_hat")],las=1,
#         names = c("beta3"),
#         main = expression(beta[3]), ylab = "Estimate")
# abline(h = beta[4], col = "black", lwd = 2)
#
# par(mfrow = c(1, 1))
#
#
#
#
# par(mfrow = c(1, 3))
#
# boxplot(mat_res_pl_df$kappa_pl,las=1,
#         main = expression(kappa["PL"]), ylab = "Estimate")
# abline(h = kappa_pl, col = "black", lwd = 2)
#
# boxplot(mat_res_pl_df$sigma2_pl,las=1,
#         main = expression(sigma["PL"]^2), ylab = "Estimate")
# abline(h = sigma2_pl, col = "black", lwd = 2)
#
# boxplot(mat_res_pl_df$sigma2_wn,las=1,
#         names = c("phi_ar1"),
#         main = expression(sigma["WN"]^2), ylab = "Estimate")
# abline(h = sigma2_wn, col = "black", lwd = 2)
# par(mfrow = c(1, 1))
#
#
#
# # Compute empirical coverage of confidence intervals for beta
#
# zval = qnorm(0.975)
# mat_res_pl_df$upper_ci_gmwmx_beta0 = mat_res_pl_df$gmwmx_beta0_hat + zval * mat_res_pl_df$gmwmx_std_beta0_hat
# mat_res_pl_df$lower_ci_gmwmx_beta0 = mat_res_pl_df$gmwmx_beta0_hat - zval * mat_res_pl_df$gmwmx_std_beta0_hat
# mat_res_pl_df$upper_ci_gmwmx_beta1 = mat_res_pl_df$gmwmx_beta1_hat + zval * mat_res_pl_df$gmwmx_std_beta1_hat
# mat_res_pl_df$lower_ci_gmwmx_beta1 = mat_res_pl_df$gmwmx_beta1_hat - zval * mat_res_pl_df$gmwmx_std_beta1_hat
# # empirical coverage of gmwmx beta
# dplyr::between(rep(beta[1], B_pl), mat_res_pl_df$lower_ci_gmwmx_beta0, mat_res_pl_df$upper_ci_gmwmx_beta0) %>% mean()
# dplyr::between(rep(beta[2], B_pl), mat_res_pl_df$lower_ci_gmwmx_beta1, mat_res_pl_df$upper_ci_gmwmx_beta1) %>% mean()
#
# # do the same for lm beta
# mat_res_pl_df$upper_ci_lm_beta0 = mat_res_pl_df$lm_beta0_hat + zval * mat_res_pl_df$lm_std_beta0_hat
# mat_res_pl_df$lower_ci_lm_beta0 = mat_res_pl_df$lm_beta0_hat - zval * mat_res_pl_df$lm_std_beta0_hat
# mat_res_pl_df$upper_ci_lm_beta1 = mat_res_pl_df$lm_beta1_hat + zval * mat_res_pl_df$lm_std_beta1_hat
# mat_res_pl_df$lower_ci_lm_beta1 = mat_res_pl_df$lm_beta1_hat - zval * mat_res_pl_df$lm_std_beta1_hat
# dplyr::between(rep(beta[1], B_pl), mat_res_pl_df$lower_ci_lm_beta0, mat_res_pl_df$upper_ci_lm_beta0) %>% mean()
# dplyr::between(rep(beta[2], B_pl), mat_res_pl_df$lower_ci_lm_beta1, mat_res_pl_df$upper_ci_lm_beta1) %>% mean()
#
#
# # Example 3: White noise + flicker
# sigma2_wn_fl <- 15
# sigma2_fl <- 10
# eps_fl = generate(wn(sigma2_wn_fl) + flicker(sigma2 = sigma2_fl),
#                   n = n, seed = 123)
# plot(eps_fl$series, type = "l")
# plot(wv::wvar(eps_fl$series))
# y_fl = X %*% beta + eps_fl$series
#
# fit_fl = gmwmx2_new(X = X, y = y_fl, model = wn() + flicker())
# fit_fl
#
#
# Monte Carlo coverage for the flicker case:
#
#
# B_fl = 500
# mat_res_fl = matrix(NA, nrow = B_fl, ncol = 18)
# for(b in seq(B_fl)){
#   eps = generate(wn(sigma2_wn_fl) + flicker(sigma2 = sigma2_fl),
#                  n = n, seed = (123 + b))$series
#   y = X %*% beta + eps
#   fit = gmwmx2_new(X = X, y = y, model = wn() + flicker())
#   fit2 = lm(y~X[,2] + X[,3] + X[,4])
#
#   mat_res_fl[b, ] = c(fit$beta_hat, fit$std_beta_hat,
#                       summary(fit2)$coefficients[,1],
#                       summary(fit2)$coefficients[,2],
#                       fit$theta_domain$`Flicker_2`,
#                       fit$theta_domain$`White Noise_1`)
#   # cat("Iteration ", b, " \n")
# }
#
# mat_res_fl_df = as.data.frame(mat_res_fl)
# colnames(mat_res_fl_df) = c("gmwmx_beta0_hat", "gmwmx_beta1_hat",
#                             "gmwmx_beta2_hat", "gmwmx_beta3_hat",
#                             "gmwmx_std_beta0_hat", "gmwmx_std_beta1_hat",
#                             "gmwmx_std_beta2_hat", "gmwmx_std_beta3_hat",
#                             "lm_beta0_hat", "lm_beta1_hat", "lm_beta2_hat", "lm_beta3_hat",
#                             "lm_std_beta0_hat", "lm_std_beta1_hat",
#                             "lm_std_beta2_hat", "lm_std_beta3_hat",
#                             "sigma2_fl" ,"sigma2_wn")
#
#
# # Plot empirical distributions (beta and stochastic parameters)
#
#
# par(mfrow = c(1, 4))
# boxplot(mat_res_df[, c("gmwmx_beta0_hat")],las=1,
#         names = c("beta0"),
#         main = expression(beta[0]), ylab = "Estimate")
# abline(h = beta[1], col = "black", lwd = 2)
#
# boxplot(mat_res_df[, c("gmwmx_beta1_hat")],las=1,
#         names = c("beta1"),
#         main = expression(beta[1]), ylab = "Estimate")
# abline(h = beta[2], col = "black", lwd = 2)
#
# boxplot(mat_res_df[, c("gmwmx_beta2_hat")],las=1,
#         names = c("beta2"),
#         main = expression(beta[2]), ylab = "Estimate")
# abline(h = beta[3], col = "black", lwd = 2)
#
# boxplot(mat_res_df[, c("gmwmx_beta3_hat")],las=1,
#         names = c("beta3"),
#         main = expression(beta[3]), ylab = "Estimate")
# abline(h = beta[4], col = "black", lwd = 2)
#
# par(mfrow = c(1, 1))
#
#
#
#
# par(mfrow = c(1, 2))
#
# boxplot(mat_res_fl_df$sigma2_fl,las=1,
#         main = expression(sigma["FL"]^2), ylab = "Estimate")
# abline(h = sigma2_fl, col = "black", lwd = 2)
#
# boxplot(mat_res_fl_df$sigma2_wn,las=1,
#         names = c("phi_ar1"),
#         main = expression(sigma["WN"]^2), ylab = "Estimate")
# abline(h = sigma2_wn_fl, col = "black", lwd = 2)
# par(mfrow = c(1, 1))
#
#
# # Interpreting the outputs
#
# The `gmwmx2_new()` object contains:
#
#   1. `beta_hat` and `std_beta_hat`: regression estimates and standard errors.
# 2. `theta_domain`: estimated stochastic parameters (e.g., AR(1), WN variance, power-law
#                                                     parameters, flicker variance).
# 3. Diagnostic information about the fitted wavelet variance and optimization output.
#
# The Monte Carlo results provide empirical coverage for the nominal 95% confidence
# intervals. This helps assess whether standard errors from `gmwmx2_new()` are well
# calibrated under different noise models.
#
# # Missing values
#
# `gmwmx2_new()` can also be used with missing observations (set entries of `y` to `NA`),
# provided the missingness model is properly specified.
