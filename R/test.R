#
# # n = 30*365
# n  = 1000
# n
# X = matrix(NA, nrow=n, ncol=4)
# # intercept
# X[,1] = 1
# # # trend
# X[,2] = 1:n
# # # add a sin signal
# omega_1 <- (1 / 365.25) * 2 * pi
# X[, 3] <- sin((1:n) * omega_1)
# X[, 4] <- cos((1:n) * omega_1)
# beta = c(1, .01, 2, 1.5)
# # plot signal
# plot(x = X[,2], y = X %*% beta, type='l', main='Signal (X %*% beta)', xlab='Time', ylab='Signal')
# eps = generate(ar1(phi=0.95, sigma2=20) + wn(20), n=n, seed = 123)
# plot(eps)
# plot(wv::wvar(eps$series))
# y = X%*% beta + eps$series
#
# fit1 = gmwmx2_new(X = X, y = y, model = wn(20) + ar1(phi=.9,sigma2 =  20) )
# fit1
# fit2 = gmwmx2_new(X = X, y = y, model = wn() + ar1() )
# fit2
#
# B = 500
# mat_res = matrix(NA, nrow=B, ncol=19)
# for(b in seq(B)){
#   eps = generate(ar1(phi=0.95, sigma2=20) + wn(20), n=n, seed = (123 + b))$series
#   # plot(wv::wvar(eps))
#   y = X %*% beta + eps
#   fit = gmwmx2_new_no_missing(X = X, y = y, model = wn() + ar1() )
#   # mispecified model assuming white noise as the stochastic model
#   fit2 = lm(y~X[,2] + X[,3] + X[,4])
#
#   mat_res[b, ] = c(fit$beta_hat, fit$std_beta_hat,
#                    summary(fit2)$coefficients[,1],
#                    summary(fit2)$coefficients[,2],
#                    fit$theta_domain$`AR(1)_2`,
#                    fit$theta_domain$`White Noise_1`)
#   cat("Iteration ", b, " \n")
# }
# #
# # # compute empirical coverage
# mat_res_df = as.data.frame(mat_res)
# colnames(mat_res_df) = c("gmwmx_beta0_hat", "gmwmx_beta1_hat",
#                          "gmwmx_beta2_hat", "gmwmx_beta3_hat",
#                          "gmwmx_std_beta0_hat", "gmwmx_std_beta1_hat",
#                          "gmwmx_std_beta2_hat", "gmwmx_std_beta3_hat",
#                          "lm_beta0_hat", "lm_beta1_hat", "lm_beta2_hat", "lm_beta3_hat",
#                          "lm_std_beta0_hat", "lm_std_beta1_hat", "lm_std_beta2_hat", "lm_std_beta3_hat",
#                          "phi_ar1","sigma_2_ar1" ,"sigma_2_wn")
# zval = qnorm(0.975)
# mat_res_df$upper_ci_gmwmx_beta0 = mat_res_df$gmwmx_beta0_hat + zval * mat_res_df$gmwmx_std_beta0_hat
# mat_res_df$lower_ci_gmwmx_beta0 = mat_res_df$gmwmx_beta0_hat - zval * mat_res_df$gmwmx_std_beta0_hat
# mat_res_df$upper_ci_gmwmx_beta1 = mat_res_df$gmwmx_beta1_hat + zval * mat_res_df$gmwmx_std_beta1_hat
# mat_res_df$lower_ci_gmwmx_beta1 = mat_res_df$gmwmx_beta1_hat - zval * mat_res_df$gmwmx_std_beta1_hat
# # empirical coverage of gmwmx beta
# dplyr::between(rep(beta[1], 500), mat_res_df$lower_ci_gmwmx_beta0, mat_res_df$upper_ci_gmwmx_beta0) %>% mean()
# dplyr::between(rep(beta[2], 500), mat_res_df$lower_ci_gmwmx_beta1, mat_res_df$upper_ci_gmwmx_beta1) %>% mean()
#
# # # do the same for lm beta
# mat_res_df$upper_ci_lm_beta0 = mat_res_df$lm_beta0_hat + zval * mat_res_df$lm_std_beta0_hat
# mat_res_df$lower_ci_lm_beta0 = mat_res_df$lm_beta0_hat - zval * mat_res_df$lm_std_beta0_hat
# mat_res_df$upper_ci_lm_beta1 = mat_res_df$lm_beta1_hat + zval * mat_res_df$lm_std_beta1_hat
# mat_res_df$lower_ci_lm_beta1 = mat_res_df$lm_beta1_hat - zval * mat_res_df$lm_std_beta1_hat
# dplyr::between(rep(beta[1], 500), mat_res_df$lower_ci_lm_beta0, mat_res_df$upper_ci_lm_beta0) %>% mean()
# dplyr::between(rep(beta[2], 500), mat_res_df$lower_ci_lm_beta1, mat_res_df$upper_ci_lm_beta1) %>% mean()
#
#
# # ---- simulation: white noise + stationary powerlaw (kappa = -0.8) ----
# kappa_pl <- -0.8
# sigma2_wn <- 15
# sigma2_pl <- 10
# n=1000
# eps_pl = generate(wn(sigma2_wn) + pl(kappa = kappa_pl, sigma2 = sigma2_pl), n = n, seed = 123)
# plot(eps_pl$series, type="l")
# plot(wv::wvar(eps_pl$series))
# y_pl = X %*% beta + eps_pl$series
#
# fit_pl = gmwmx2_new(X = X, y = y_pl, model = wn() + pl())
# fit_pl
#
# B_pl = 500
# mat_res_pl = matrix(NA, nrow = B_pl, ncol = 19)
# for(b in seq(B_pl)){
#   eps = generate(wn(sigma2_wn) + pl(kappa = kappa_pl, sigma2 = sigma2_pl), n = n, seed = (123 + b))$series
#   y = X %*% beta + eps
#   fit = gmwmx2_new_no_missing(X = X, y = y, model = wn() + pl())
#   # mispecified model assuming white noise as the stochastic model
#   fit2 = lm(y~X[,2] + X[,3] + X[,4])
#
#   mat_res_pl[b, ] = c(fit$beta_hat, fit$std_beta_hat,
#                       summary(fit2)$coefficients[,1],
#                       summary(fit2)$coefficients[,2],
#                       fit$theta_domain$`Stationary PowerLaw_2`,
#                       fit$theta_domain$`White Noise_1`)
#   cat("Iteration ", b, " \n")
# }
# #
# # # compute empirical coverage
# mat_res_pl_df = as.data.frame(mat_res_pl)
# colnames(mat_res_pl_df) = c("gmwmx_beta0_hat", "gmwmx_beta1_hat",
#                             "gmwmx_beta2_hat", "gmwmx_beta3_hat",
#                             "gmwmx_std_beta0_hat", "gmwmx_std_beta1_hat",
#                             "gmwmx_std_beta2_hat", "gmwmx_std_beta3_hat",
#                             "lm_beta0_hat", "lm_beta1_hat", "lm_beta2_hat", "lm_beta3_hat",
#                             "lm_std_beta0_hat", "lm_std_beta1_hat", "lm_std_beta2_hat", "lm_std_beta3_hat",
#                             "kappa_pl","sigma2_pl" ,"sigma2_wn")
# zval = qnorm(0.975)
# mat_res_pl_df$upper_ci_gmwmx_beta0 = mat_res_pl_df$gmwmx_beta0_hat + zval * mat_res_pl_df$gmwmx_std_beta0_hat
# mat_res_pl_df$lower_ci_gmwmx_beta0 = mat_res_pl_df$gmwmx_beta0_hat - zval * mat_res_pl_df$gmwmx_std_beta0_hat
# mat_res_pl_df$upper_ci_gmwmx_beta1 = mat_res_pl_df$gmwmx_beta1_hat + zval * mat_res_pl_df$gmwmx_std_beta1_hat
# mat_res_pl_df$lower_ci_gmwmx_beta1 = mat_res_pl_df$gmwmx_beta1_hat - zval * mat_res_pl_df$gmwmx_std_beta1_hat
# # empirical coverage of gmwmx beta
# dplyr::between(rep(beta[1], 500), mat_res_pl_df$lower_ci_gmwmx_beta0, mat_res_pl_df$upper_ci_gmwmx_beta0) %>% mean()
# dplyr::between(rep(beta[2], 500), mat_res_pl_df$lower_ci_gmwmx_beta1, mat_res_pl_df$upper_ci_gmwmx_beta1) %>% mean()
#
# # # do the same for lm beta
# mat_res_pl_df$upper_ci_lm_beta0 = mat_res_pl_df$lm_beta0_hat + zval * mat_res_pl_df$lm_std_beta0_hat
# mat_res_pl_df$lower_ci_lm_beta0 = mat_res_pl_df$lm_beta0_hat - zval * mat_res_pl_df$lm_std_beta0_hat
# mat_res_pl_df$upper_ci_lm_beta1 = mat_res_pl_df$lm_beta1_hat + zval * mat_res_pl_df$lm_std_beta1_hat
# mat_res_pl_df$lower_ci_lm_beta1 = mat_res_pl_df$lm_beta1_hat - zval * mat_res_pl_df$lm_std_beta1_hat
# dplyr::between(rep(beta[1], 500), mat_res_pl_df$lower_ci_lm_beta0, mat_res_pl_df$upper_ci_lm_beta0) %>% mean()
# dplyr::between(rep(beta[2], 500), mat_res_pl_df$lower_ci_lm_beta1, mat_res_pl_df$upper_ci_lm_beta1) %>% mean()
#
# #
# # # ---- simulation: white noise + flicker ----
# sigma2_wn_fl <- 15
# sigma2_fl <- 10
# n = 1000
# eps_fl = generate(wn(sigma2_wn_fl) + flicker(sigma2 = sigma2_fl), n = n, seed = 123)
# plot(eps_fl$series, type="l")
# plot(wv::wvar(eps_fl$series))
# y_fl = X %*% beta + eps_fl$series
#
# fit_fl = gmwmx2_new(X = X, y = y_fl, model = wn() + flicker())
# fit_fl
#
# B_fl = 500
# mat_res_fl = matrix(NA, nrow = B_fl, ncol = 18)
# for(b in seq(B_fl)){
#   eps = generate(wn(sigma2_wn_fl) + flicker(sigma2 = sigma2_fl), n = n, seed = (123 + b))$series
#   y = X %*% beta + eps
#   fit = gmwmx2_new_no_missing(X = X, y = y, model = wn() + flicker())
#   # mispecified model assuming white noise as the stochastic model
#   fit2 = lm(y~X[,2] + X[,3] + X[,4])
#
#   mat_res_fl[b, ] = c(fit$beta_hat, fit$std_beta_hat,
#                       summary(fit2)$coefficients[,1],
#                       summary(fit2)$coefficients[,2],
#                       fit$theta_domain$`Flicker_2`,
#                       fit$theta_domain$`White Noise_1`)
#   cat("Iteration ", b, " \n")
# }
#
# # compute empirical coverage
# mat_res_fl_df = as.data.frame(mat_res_fl)
# colnames(mat_res_fl_df) = c("gmwmx_beta0_hat", "gmwmx_beta1_hat",
#                             "gmwmx_beta2_hat", "gmwmx_beta3_hat",
#                             "gmwmx_std_beta0_hat", "gmwmx_std_beta1_hat",
#                             "gmwmx_std_beta2_hat", "gmwmx_std_beta3_hat",
#                             "lm_beta0_hat", "lm_beta1_hat", "lm_beta2_hat", "lm_beta3_hat",
#                             "lm_std_beta0_hat", "lm_std_beta1_hat", "lm_std_beta2_hat", "lm_std_beta3_hat",
#                             "sigma2_fl" ,"sigma2_wn")
# zval = qnorm(0.975)
# mat_res_fl_df$upper_ci_gmwmx_beta0 = mat_res_fl_df$gmwmx_beta0_hat + zval * mat_res_fl_df$gmwmx_std_beta0_hat
# mat_res_fl_df$lower_ci_gmwmx_beta0 = mat_res_fl_df$gmwmx_beta0_hat - zval * mat_res_fl_df$gmwmx_std_beta0_hat
# mat_res_fl_df$upper_ci_gmwmx_beta1 = mat_res_fl_df$gmwmx_beta1_hat + zval * mat_res_fl_df$gmwmx_std_beta1_hat
# mat_res_fl_df$lower_ci_gmwmx_beta1 = mat_res_fl_df$gmwmx_beta1_hat - zval * mat_res_fl_df$gmwmx_std_beta1_hat
# # empirical coverage of gmwmx beta
# dplyr::between(rep(beta[1], 500), mat_res_fl_df$lower_ci_gmwmx_beta0, mat_res_fl_df$upper_ci_gmwmx_beta0) %>% mean()
# dplyr::between(rep(beta[2], 500), mat_res_fl_df$lower_ci_gmwmx_beta1, mat_res_fl_df$upper_ci_gmwmx_beta1) %>% mean()
#
# # do the same for lm beta
# mat_res_fl_df$upper_ci_lm_beta0 = mat_res_fl_df$lm_beta0_hat + zval * mat_res_fl_df$lm_std_beta0_hat
# mat_res_fl_df$lower_ci_lm_beta0 = mat_res_fl_df$lm_beta0_hat - zval * mat_res_fl_df$lm_std_beta0_hat
# mat_res_fl_df$upper_ci_lm_beta1 = mat_res_fl_df$lm_beta1_hat + zval * mat_res_fl_df$lm_std_beta1_hat
# mat_res_fl_df$lower_ci_lm_beta1 = mat_res_fl_df$lm_beta1_hat - zval * mat_res_fl_df$lm_std_beta1_hat
# dplyr::between(rep(beta[1], 500), mat_res_fl_df$lower_ci_lm_beta0, mat_res_fl_df$upper_ci_lm_beta0) %>% mean()
# dplyr::between(rep(beta[2], 500), mat_res_fl_df$lower_ci_lm_beta1, mat_res_fl_df$upper_ci_lm_beta1) %>% mean()
#
