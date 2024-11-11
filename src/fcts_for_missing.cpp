#include <numeric>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::export]]
arma::vec estimate_p1_p2_mle_cpp(arma::vec omega) {
  int n = omega.n_elem;

  double n1 = sum(omega.subvec(1, n-1) % (1 - omega.subvec(0, n-2)));
  double n2 = sum((1 - omega.subvec(1, n-1)) % (1 - omega.subvec(0, n-2)));
  double n3 = sum(omega.subvec(1, n-1) % omega.subvec(0, n-2));
  double n4 = sum((1 - omega.subvec(1, n-1)) % omega.subvec(0, n-2));

  double p1_hat = n4 / (n3 + n4);
  double p2_hat = n1 / (n1 + n2);
  arma::vec result(2);
  result(0) = p1_hat;
  result(1) = p2_hat;
  return result;
}




// Function to generate vector of autocovariance of omega process, A markov process that takes values 0 or 1 indicating respectively missingness or presence of the observation.
// [[Rcpp::export]]
arma::vec create_vec_theo_autocov_omega_cpp(double p1, double p2, int n) {
  double pstar = p2 / (p1 + p2);
  int hmax = n + 1;

  // Initialize Armadillo vectors
  arma::vec v(hmax, arma::fill::zeros);
  arma::vec q(hmax, arma::fill::zeros);

  q(0) = 1.0;
  v(0) = 0.0;

  // Compute vectors v and q
  for (int h = 0; h < n; ++h) {
    v(h + 1) = p2 * q(h) + (1 - p2) * v(h);
    q(h + 1) = (1 - p1) * q(h) + p1 * v(h);
  }

  // Compute autocovariance vector
  arma::vec autocov_vec = pstar * q - pstar * pstar;

  // Update variance (autocov k = 0)
  autocov_vec(0) = (p2 / (p1 + p2)) - std::pow((p2 / (p1 + p2)), 2);

  return autocov_vec.subvec(0,n-1);
}
