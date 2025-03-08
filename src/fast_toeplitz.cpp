#include <numeric>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;


// Function to generate a Toeplitz matrix from a vector
// [[Rcpp::export]]
arma::mat fast_toeplitz_matrix_from_vector_cpp(const arma::vec& v) {
  return arma::toeplitz(v);
}
