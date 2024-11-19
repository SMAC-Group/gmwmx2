#include <RcppArmadillo.h>
#include <numeric>
using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]


// Compute the stationary power law autocovariance vector
// @param kappa The spectral index kappa that needs to be > -1 so that the autocovariance of the stationary powerlaw process is defined.
// @param sigma2 The variance of the stationary powerlaw process.
// @param n An integer specifying the length of the signal
// [[Rcpp::export]]
arma::vec powerlaw_autocovariance(const double kappa , const double sigma2, const int n) {
  // construct vector of the autocovariance function as given in

  arma::vec acf(n);

  acf(0) = ::tgamma(1.0 + kappa) / pow(::tgamma(1.0 + kappa/2), 2) * sigma2;
  for (int i = 1; i < n; ++i) {
    acf(i) = (-kappa/2 + i - 1.0) * acf(i - 1) / (i +kappa/2);
  }

  return acf;
}



