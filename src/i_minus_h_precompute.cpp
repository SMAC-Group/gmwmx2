#include <numeric>
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;


// Function to calculate mean of all diagonals of a matrix
// [[Rcpp::export]]
arma::vec sum_all_upper_diagonals(const arma::mat& matrix) {

  // Ensure that the input is a square matrix
  if (matrix.n_rows != matrix.n_cols) {
    throw std::invalid_argument("Input must be a square matrix");
  }

  // Get the size of the matrix
  int size = matrix.n_rows;

  // Vector to store mean values for each diagonal
  arma::vec diagonal_sums(size, arma::fill::zeros);

  // Iterate over each diagonal
  for (int i = 0; i < size; ++i) {
    arma::vec diagonal_i = diagvec(matrix, i);
    double diagonal_sum = arma::accu(diagonal_i);
    diagonal_sums(i) = diagonal_sum;
  }

  return diagonal_sums;
}




// [[Rcpp::export]]
arma::vec linear_interp_cpp(const arma::uvec& x, const arma::vec& y, const arma::vec& xout) {
  int n = xout.n_elem;
  arma::vec yout(n);

  for (int i = 0; i < n; ++i) {
    double xi = xout[i];
    if (xi <= x[0]) {
      yout[i] = y[0];
    } else if (xi >= x[x.n_elem - 1]) {
      yout[i] = y[y.n_elem - 1];
    } else {
      for (unsigned j = 1; j < x.n_elem; ++j) {
        if (xi < x[j]) {
          double t = (xi - x[j-1]) / (x[j] - x[j-1]);
          yout[i] = y[j-1] + t * (y[j] - y[j-1]);
          break;
        }
      }
    }
  }

  return yout;
}


// [[Rcpp::export]]
Rcpp::List pre_compute_quantities_on_D_only_required_smarter_cpp(const arma::mat& D, std::string approx_type="3") {
  // Extract dimensions
  int n = D.n_rows;

  // Compute number of points
  int num_points = 0;

  if (approx_type == "1") {
    num_points = std::floor(std::log2(n));
  } else if (approx_type == "2") {
    num_points = std::floor(std::log2(n)) * 2;
  } else if (approx_type == "3") {
    num_points = std::floor(std::log2(n)) * 3;
  } else if (approx_type == "all") {
    num_points = n;
  } else {
    Rcpp::stop("Invalid approx_type argument. Must be '1', '2', '3', or 'all'.");
  }

  // Create a sequence from 1 to n with length `num_points`
  arma::vec seq = arma::round(arma::linspace(1, n, num_points)) - 1;
  arma::uvec indices_l_values = arma::conv_to<arma::uvec>::from(seq);

  // Initialize matrices
  arma::mat mat_D_q_term_1_new = arma::zeros(num_points, n);
  arma::mat mat_D_q_term_2_new = arma::zeros(num_points, n);

  // Fill matrices
  for (int j = 0; j < n; ++j) {
    arma::vec diag_j = D.diag(j);
    arma::vec cumsum_diag_j = arma::cumsum(diag_j); // Compute cumulative sum of diag_j

    for (size_t l_index = 0; l_index < indices_l_values.n_elem; ++l_index) {
      int l_value = indices_l_values[l_index];

      if (l_value >= j) {
        mat_D_q_term_2_new(l_index, l_value - j) = cumsum_diag_j(n - l_value - 1);
      }

      if (l_value < static_cast<int>(diag_j.n_elem) && (j > 0)) {
        mat_D_q_term_1_new(l_index, j-1) = cumsum_diag_j(diag_j.n_elem - l_value - 1);
      }
    }
  }

  // Compute the vector required for term 3
  arma::vec sum_on_sub_diag_of_D = sum_all_upper_diagonals(D);

  // Return the two matrices and the vector
  return Rcpp::List::create(
    Rcpp::Named("mat_D_q_term_1") = mat_D_q_term_1_new,
    Rcpp::Named("mat_D_q_term_2") = mat_D_q_term_2_new,
    Rcpp::Named("sum_on_sub_diag_of_D") = sum_on_sub_diag_of_D
  );
}




// [[Rcpp::export]]
arma::vec compute_all_mean_diag_fast_w_linear_interp_only_required_cpp(const arma::mat& mat_D_q_term_1,
                                                         const arma::mat& mat_D_q_term_2,
                                                         const arma::vec& sum_on_sub_diag_of_D,
                                                         const arma::vec& vec_autocov,
                                                         std::string approx_type="3") {

  int n = vec_autocov.n_elem;
  int num_points = 0;

  if(approx_type == "1"){
    num_points = std::floor(std::log2(n));
  }else if(approx_type == "2"){
    num_points = std::floor(std::log2(n)) * 2;
  }else if(approx_type == "3"){
    num_points = std::floor(std::log2(n)) * 3;
  }else if(approx_type == "all"){
    num_points = n;
  } else {
    Rcpp::stop("Invalid approx_type argument. Must be '1', '2', '3', or 'all'.");
  }

  // Create a sequence from 1 to n with length `num_points`
  arma::vec seq = arma::round(arma::linspace(1, n, num_points))-1;
  arma::uvec indices_l_values = arma::conv_to<arma::uvec>::from(seq);

  arma::vec vec_values_at_indices(indices_l_values.n_elem, arma::fill::zeros);
  //
  for (arma::uword l_index = 0; l_index < indices_l_values.n_elem; ++l_index) {
    int l = indices_l_values[l_index];
    int m = n - l;
    if (m == 1) {
      double term_2 = 0.0;
      for (int j = 0; j <= l; ++j) {
        term_2 += mat_D_q_term_2(l_index, j) * vec_autocov[std::abs(-j)];
      }
      vec_values_at_indices[l_index] = term_2 / (n - l);
    } else {
      double term_1 = 0.0;
      for (int j = 1; j < m; ++j) {
        term_1 += vec_autocov[std::abs(-l - j)] * mat_D_q_term_1(l_index, j-1);
      }
      double term_2 = 0.0;
      for (int j = 0; j <= l; ++j) {
        term_2 += mat_D_q_term_2(l_index, j) * vec_autocov[std::abs(-j)];
      }
      double term_3 = 0.0;
      for (int j = 1; j < m; ++j) {
        term_3 += vec_autocov[std::abs(j)] * sum_on_sub_diag_of_D[j + l];
      }
      vec_values_at_indices[l_index] = (term_1 + term_2 + term_3) / (n - l);
    }
  }
  //
  arma::vec errors = vec_values_at_indices - vec_autocov.elem(indices_l_values);

  //  Perform linear interpolation of the error over the entire range
  arma::vec xout =arma::regspace<arma::vec>(1, 1, n);
  arma::vec interpolated_errors =  linear_interp_cpp(indices_l_values+1, errors, xout );
  // correct vector
  arma::vec vec_autocov_approx = vec_autocov + interpolated_errors;
  // output
  return vec_autocov_approx;
}




