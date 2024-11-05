#include <numeric>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::vec compute_h_cpp(double kappa, int N) {
  arma::vec vec_h(N);
  double h_prev = 1.0;
  double val_1 = kappa/2.0 + 1.0;
  for (int i = 1; i < N; ++i) {
    double h = (i - val_1) / i * h_prev;  // Eq. (25)
    vec_h(i) = h;
    h_prev = h;
  }

  vec_h(0) = 1.0;
  return vec_h;
}

// [[Rcpp::export]]
arma::vec vec_mean_autocov_powerlaw( double kappa, int N) {
  arma::vec vec_h = compute_h_cpp(kappa, N);
  int n = vec_h.n_elem;

  // Initialize variables
  arma::vec vec_average_autocov2(n, arma::fill::zeros);
  // vec_average_autocov2[0] = val_sum / n;

  // pass through diagonals
  for (int k = 0; k < n ; ++k) {
    double val_sum = 0.0;
    // Compute val_sum
    for (int i = 0; i < n - k; ++i) {
      int alpha_i = (n - i - k);
      double term = vec_h[i] * vec_h[i + k];
      val_sum += alpha_i * term;
    }
    // Assign the result to vec_average_autocov2[k+1]
    vec_average_autocov2[k] = val_sum / (n - k);
  }
  return(vec_average_autocov2);
}


// [[Rcpp::export]]
arma::mat var_cov_powerlaw_cpp(double sigma2, double kappa, int n) {
  arma::vec vec_h = compute_h_cpp(kappa, n);
  arma::mat sigma_mat(n, n, arma::fill::zeros);
  arma::vec vec_diag_k(n);

  //  for each sub diagonal up to the upper right corner
  for (int k = 0; k < n; ++k) {
    // re initialize sum
    double val = 0.0;

    // pass throught sum of product of h_{i} and h_{i+k}
    for (int i = 0; i < (n - k); ++i) {
      val += vec_h(i) * vec_h(i + k);
      vec_diag_k(i) = val;
    }

    // Assign values to the upper diagonal
    sigma_mat.diag(k) = vec_diag_k.head(n - k);
  }

  // multiply by sigma2
  sigma_mat = sigma_mat * sigma2;

  // Create symmetric matrix from upper triangular matrix
  arma::mat symmetric_matrix = arma::symmatu(sigma_mat);

  // // Create symmetric matrix
  // arma::mat symmetric_matrix = sigma_mat + sigma_mat.t();
  // // divide per two the value on the main diagonal are they are summed up two times in the previous operation
  // symmetric_matrix.diag() *= 0.5;

  return symmetric_matrix;
}


// [[Rcpp::export]]
arma::vec compute_power_of_a_base(int x, int J) {
  arma::vec power = arma::regspace(0, J);
  arma::vec out(power.n_elem);
  for (size_t i = 0; i < power.n_elem; ++i) {
    out(i) = std::pow(x, power(i));
  }
  return out;
}





// [[Rcpp::export]]
arma::vec autocovariance_to_wv(const arma::vec& acf, const arma::vec& tau) {
  // Compute max scale
  double J = std::log10(tau(tau.n_elem - 1)) / std::log10(2.0);

  // Index first element of acf
  double var_process = acf(0);
  arma::vec autocorr = acf.subvec(1, acf.n_elem - 1) / var_process;
  arma::vec ms = compute_power_of_a_base(2, J);

  // Initialize vector for theoretical wavelet variance
  arma::vec theo_wv(J);

  for (int j = 1; j <= J; ++j) {
    double m = ms(j - 1);
    double inter = m * (1 - autocorr(m - 1));

    if (m > 1) {
      for (int i = 1; i <= m - 1; ++i) {
        inter += i * (2 * autocorr(m - i - 1) - autocorr(i - 1) - autocorr(2 * m - i - 1));
      }
    }

    theo_wv(j - 1) = inter / (m * m) * var_process / 2.0;
  }

  return theo_wv;
}




/*** R


*/
