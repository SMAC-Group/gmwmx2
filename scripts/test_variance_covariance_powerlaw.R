rm(list=ls())
library(gmwmx2)

n=5

#white noise
kappa=0
var_cov_powerlaw_cpp(kappa = kappa, sigma2 = 1, n = n)

# rw also called red noise
kappa=-2
compute_h_cpp(kappa =kappa,N = n )
var_cov_powerlaw_cpp(kappa = kappa, sigma2 = 1, n = n)

# if kappa is smaller than -1, becomes non stationary
kappa=-.6
test = var_cov_powerlaw_cpp(kappa = kappa, sigma2 = 1, n = 5000)
plot(mgcv::sdiag(test))
vec_autocov_stationary_powerlaw = powerlaw_autocovariance(kappa=kappa, sigma2 = 1, n = n)
abline(h=vec_autocov_stationary_powerlaw[1])

# stationary powerlaw exist if kappa > -1
# we gonna use the softplus function
transform_fun <- function(x) {
  exp(x) - 1
}

inv_transform_fun <- function(x) log(x + 1)
inv_transform_fun(transform_fun(100))
# matrix from


kappa=-.5
vec_autocov_stationary_powerlaw = powerlaw_autocovariance(kappa=kappa, sigma2 = 1, n = n)
mat_stationary_powerlaw = toeplitz(as.vector(vec_autocov_stationary_powerlaw))

mat_powerlaw = var_cov_powerlaw_cpp(kappa = kappa, sigma2 = 1, n = n)
