#include <numeric>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;




// [[Rcpp::export]]
arma::vec gen_flicker(int N, double sigma) {
  double kappa=-1;
  
  // // Set seed for reproducibility
  // arma::arma_rng::set_seed(0);
  
  // Flicker noise parameters
  double h_prev = 1.0;
  arma::vec h = arma::zeros<arma::vec>(2 * N);
  h(0) = 1.0;  // Eq. (25)
  
  // Generate h values
  for (int i = 1; i < N; ++i) {
    h(i) = (i - kappa / 2 - 1) / static_cast<double>(i) * h_prev;  // Eq. (25)
    h_prev = h(i);
  }
  
  // Generate v values
  arma::vec v = arma::zeros<arma::vec>(2 * N);
  v.subvec(0, N - 1) = arma::randn<arma::vec>(N) * sigma;
  
  // Compute w values
  arma::cx_vec fft_v = arma::fft(v);
  arma::cx_vec fft_h = arma::fft(h);
  arma::cx_vec w1 = arma::ifft(fft_v % fft_h);
  
  // Extract real part of w
  arma::vec w2 = arma::real(arma::vectorise(w1));
  
  // Extract the first N elements of w
  arma::vec y = w2.subvec(0, N - 1);
  
  return y;
}






// function to compute cov mat of flicker
// [[Rcpp::export]]
arma::mat create_var_cov_flicker(double sigma_pl,  int N) {
  // use var_cov_powerlaw_cpp, more general and faster !
  double kappa = -1;
  arma::mat U(N, N, arma::fill::eye); 
  
  arma::vec h(N);
  double h_prev = 1.0;
  
  for (int i = 1; i < N; ++i) {
    h(i) = (i - kappa / 2 - 1) / static_cast<double>(i) * h_prev; // Eq. (25)
    for (int j = 0; j < N - i; ++j) {
      U(j, j + i) = h(i);
    }
    h_prev = h(i);
  }
  
  U *= sigma_pl; // scale noise
  // long huge matrix multiplicaiton
  return trans(U) * U; // Eq. (26)
}


// [[Rcpp::export]]
arma::mat create_U_flicker(int N) {
  double kappa = -1;
  arma::mat U(N, N, arma::fill::eye); 
  
  arma::vec h(N);
  double h_prev = 1.0;
  
  for (int i = 1; i < N; ++i) {
    h(i) = (i - kappa / 2 - 1) / static_cast<double>(i) * h_prev; // Eq. (25)
    for (int j = 0; j < N - i; ++j) {
      U(j, j + i) = h(i);
    }
    h_prev = h(i);
  }
  
  // U *= sigma_pl; // scale noise
  
  // return trans(U) * U; // Eq. (26)
  return(U);
}


// Function to calculate mean of all diagonals of a matrix
// [[Rcpp::export]]
arma::vec mean_all_diagonals(arma::mat matrix) {
  
  // Ensure that the input is a square matrix
  if (matrix.n_rows != matrix.n_cols) {
    throw std::invalid_argument("Input must be a square matrix");
  }
  
  // Get the size of the matrix
  int size = matrix.n_rows;
  
  // Vector to store mean values for each diagonal
  arma::vec diagonal_means(size, arma::fill::zeros);
  
  // Iterate over each diagonal
  for (int i = 0; i < size; ++i) {
    arma::vec diagonal_1 = diagvec(matrix, i);
    arma::vec diagonal_2 = diagvec(matrix, -i);
    arma::vec diagonal = join_cols(diagonal_1, diagonal_2);
    double diagonal_mean = arma::mean(diagonal);
    diagonal_means(i) = diagonal_mean;
  }
  
  return diagonal_means;
}



// Function to calculate mean of all diagonals of a matrix
// [[Rcpp::export]]
arma::vec mean_all_upper_diagonals(arma::mat matrix) {
  
  // Ensure that the input is a square matrix
  if (matrix.n_rows != matrix.n_cols) {
    throw std::invalid_argument("Input must be a square matrix");
  }
  
  // Get the size of the matrix
  int size = matrix.n_rows;
  
  // Vector to store mean values for each diagonal
  arma::vec diagonal_means(size, arma::fill::zeros);
  
  // Iterate over each diagonal
  for (int i = 0; i < size; ++i) {
    arma::vec diagonal_1 = diagvec(matrix, i);
    double diagonal_mean = arma::mean(diagonal_1);
    diagonal_means(i) = diagonal_mean;
  }
  
  return diagonal_means;
}



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
arma::vec return_cumsum_product_of_shifted_vector_cpp(const arma::vec& vec_x, int shift) {
  int n = vec_x.n_elem;
  arma::vec vec_y = vec_x;
  arma::vec vec_x_tail = vec_x.tail(n - shift);
  arma::vec vec_y_head = vec_y.head(n - shift);
  arma::vec result = arma::cumsum(vec_x_tail % vec_y_head);
  return result;
}



// Function to compute vec_mean_autocov_per_diag
// [[Rcpp::export]]
arma::vec vec_mean_autocov_per_diag_cpp( double kappa, int N) {
  arma::vec vec_h = compute_h_cpp(kappa, N);
  int n = vec_h.n_elem;
  
  arma::vec vec_mean_autocov_half(n);
  vec_mean_autocov_half(0) = arma::mean(arma::cumsum(vec_h % vec_h));
  
  for (int i = 1; i < n; ++i) {
    vec_mean_autocov_half(i) = arma::mean(return_cumsum_product_of_shifted_vector_cpp(vec_h, i));
  }
  return vec_mean_autocov_half;
}


// [[Rcpp::export]]
arma::vec vec_mean_autocov_per_diag_2_cpp( double kappa, int N) {
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

// [[Rcpp::depends(RcppArmadillo)]]
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
arma::vec wn_autocovariance(double sigma2_wn, int n) {
  arma::vec acf(n, fill::zeros);  // Create a vector of length n filled with zeros
  acf(0) = sigma2_wn;            // Set the first element to sigma2_wn
  return acf;
}
// Function to generate a Toeplitz matrix from a vector
// [[Rcpp::export]]
arma::mat fast_toeplitz_matrix_from_vector_cpp(const arma::vec& v) {
  return arma::toeplitz(v);
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




// [[Rcpp::export]]
arma::vec get_mean_per_diag_I_minus_H_Sigma_cpp_4(const arma::mat& mat_I_minus_H, const arma::mat& mat_Sigma) {
  // Get n
  int n = mat_I_minus_H.n_rows;
  
  // get max J
  int max_j = std::floor(log2(n))-1;
  
  // get max mean autocov to get
  int max_index_mean_autocov_to_compute = 2 * std::pow(2, max_j - 1) - 1;
  
  // Create vector
  arma::vec vec_mean_autocov_per_diag((max_index_mean_autocov_to_compute + 1), arma::fill::zeros);
  
  for (int i = 0; i <= max_index_mean_autocov_to_compute; ++i) {
    // Compute the sum and mean for the current diagonal
    vec_mean_autocov_per_diag(i) = arma::accu( mat_I_minus_H.submat(0, 0, n - i - 1, n - 1) % mat_Sigma.submat(i, 0, n - 1, n - 1)) / (n-i);
  }
  return vec_mean_autocov_per_diag;
}


// [[Rcpp::export]]
arma::vec objective_function_wn_flicker_eps_hat_cpp(const arma::vec theta,const arma::vec scales, const arma::vec wv_var, const arma::vec wv_ci_up, const arma::vec wv_ci_low, const int n, const arma::mat& Q_matrix) {
  double sigma2_wn = exp(theta(0));
  double sigma2_fl = exp(theta(1));
  
  arma::vec wn_autocov_vec = wn_autocovariance(sigma2_wn, n);
  arma::mat sigma_mat_wn = fast_toeplitz_matrix_from_vector_cpp(wn_autocov_vec); //sigma White Noise
  
  arma::mat sigma_mat_fl = var_cov_powerlaw_cpp(sigma2_fl,-1,n); //Sigma Flicker Noise
  
  arma::mat sigma_mat = sigma_mat_wn + sigma_mat_fl; //Sum both process
  
  arma::vec mean_autocov_eps_hat_vec = get_mean_per_diag_I_minus_H_Sigma_cpp_4(Q_matrix, sigma_mat);
  
  arma::vec theo_wv = autocovariance_to_wv(mean_autocov_eps_hat_vec, scales);
  arma::vec variance = wv_var;
  arma::vec ci_diff = wv_ci_low - wv_ci_up;
  arma::mat omega = diagmat(1 / square(ci_diff));
  arma::vec difference = variance - theo_wv;
  arma::vec objective = difference.t() * omega * difference;
  return objective ;
}


// [[Rcpp::export]]
arma::vec objective_function_wn_flicker_eps_hat_cpp_2(const arma::vec theta,const arma::vec scales, const arma::vec wv_var, const arma::vec wv_ci_up, const arma::vec wv_ci_low, const int n, const arma::mat& Q_matrix) {
  double sigma2_wn = exp(theta(0));
  double sigma2_fl = exp(theta(1));
  
  arma::vec wn_autocov_vec = wn_autocovariance(sigma2_wn, n);
  arma::mat sigma_mat_wn = fast_toeplitz_matrix_from_vector_cpp(wn_autocov_vec);
  
  arma::mat sigma_mat_fl = var_cov_powerlaw_cpp(sigma2_fl,-1,n);
  
  arma::mat sigma_mat = sigma_mat_wn + sigma_mat_fl;
  
  arma::mat sigma_mat_eps_hat = Q_matrix * sigma_mat * Q_matrix.t();
  arma::vec mean_autocov_eps_hat_vec = mean_all_diagonals(sigma_mat_eps_hat);
  
  arma::vec theo_wv = autocovariance_to_wv(mean_autocov_eps_hat_vec, scales);
  arma::vec variance = wv_var;
  arma::vec ci_diff = wv_ci_low - wv_ci_up;
  arma::mat omega = diagmat(1 / square(ci_diff));
  arma::vec difference = variance - theo_wv;
  arma::vec objective = difference.t() * omega * difference;
  return objective ;
}











/*** R


*/
