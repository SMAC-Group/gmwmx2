# invesitgate why matern is doing problems
# Ma <- function(x, alpha){
#   2/gamma(alpha-1/2)/2^(alpha-1/2)*abs(x)^(alpha-1/2)*besselK(abs(x), abs(alpha-1/2))
# }
# my_autocov = function(sigma2, lambda, alpha, n) {
#   autocov = c(sigma2, sigma2*Ma(lambda* (1:(n-1)), alpha = alpha))
#   return(autocov)
# }

