# n = 1000
# X = matrix(NA, nrow=n, ncol=4)
# # intercept
# X[,1] = 1
# # # trend
# X[,2] = 1:n
# # # add a sin signal
# omega_1 <- (1 / 365.25) * 2 * pi
# X[, 3] <- sin((1:n) * omega_1)
# X[, 4] <- cos((1:n) * omega_1)
# beta = c(1, 2, 3,4)
# eps = generate(ar1(phi=0.95, sigma2=20) + wn(20), n=n, seed = 123)$series
# plot(wv::wvar(eps))
# yy = X%*% beta + eps
# B = 500
# mat_res = matrix(NA, nrow=B, ncol=19)
# b=145
# eps = generate(ar1(phi=phi_ar1, sigma2=sigma2_ar1) + wn(sigma2_wn), n=n, seed = (123 + b))$series
# plot(wv::wvar(eps))
# y = X %*% beta + eps
#
# fit1 = gmwmx2_new(X = X, y = y, model = wn(20) + ar1(phi=.999,sigma2 =  20) )
# fit2 = gmwmx2_new(X = X, y = y, model = wn() + ar1() )

#
#
#
# #-----------
#
#
#   model <- wn(20) + ar1(phi = .95, sigma2 = 20)
#   beta_hat <- .lm.fit(y = y, x = X)$coefficients
#   eps_hat <- y - X %*% beta_hat
#   D <- diag(n) - X %*% solve(t(X) %*% X) %*% t(X)
#   quantities_D <- pre_compute_quantities_on_D_only_required_smarter_cpp(D, approx_type ="3")
#   model_filled <- fill_missing_parameters(model, signal = eps_hat)
#   prep <- prepare_optim_layout(model_filled)
#   wv_emp <- wv::wvar(eps_hat)
#
#   # evaluate loss at initial theta
#   loss_fn_gmwmx_no_missing(prep$theta0, model_filled, n, prep, wv_emp, quantities_D)
#
#
#
#   autocov_vec <- get_autocovariance(model_filled, n, prep$theta0, prep)
#   vec_mean_autocov_eps_hat <- compute_all_mean_diag_fast_w_linear_interp_only_required_cpp(
#     quantities_D$mat_D_q_term_1,
#     quantities_D$mat_D_q_term_2,
#     quantities_D$sum_on_sub_diag_of_D,
#     autocov_vec, approx_type = "3"
#   )
#   theo_wv <- autocovariance_to_wv(vec_mean_autocov_eps_hat, tau = wv_emp$scales)
#
#   any(!is.finite(autocov_vec))
#   any(!is.finite(vec_mean_autocov_eps_hat))
#   any(!is.finite(theo_wv))
#   any(!is.finite(wv_emp$ci_low - wv_emp$ci_high))
#
#
#
#
#
#
