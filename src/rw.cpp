#include <numeric>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat get_sigma_mat_rw(int n, double gamma2) {

  arma::mat mat(n, n, arma::fill::zeros);

  // fill the matrix with min(i,j)
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      mat(i,j) = std::min(i+1, j+1);
    }
  }
  return gamma2 * mat;
}


// [[Rcpp::export]]
arma::vec get_mean_diagonal_super_diagonals_cov_mat_rw_cpp(int n, double sigma_2_rw) {
  arma::vec vec_mean_autocovariance(n, arma::fill::zeros);

  for (int i = 0; i < n; ++i) {
    int k = i; // k = i - 1 in R (Rcpp indexing starts from 0)
    double denominator = n - k;
    double numerator = denominator * (denominator + 1) / 2.0;
    vec_mean_autocovariance[i] = (numerator / denominator) * sigma_2_rw;
  }

  return vec_mean_autocovariance;
}
