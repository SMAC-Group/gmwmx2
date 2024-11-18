#include <numeric>
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;



//  function that returns the sum of power of 2 from from to to
// [[Rcpp::export]]
int sum_of_powers_of_2(int from, int to) {
  int sum = 0;
  for (int i = from; i <= to; ++i) {
    sum += pow(2, i);
  }
  return sum;
}

// return power of 2
int pow2(int x) {
  return pow(x, 2);
}



// function that compute the autocovariance of the Wavelet coefficient at scale 1 from the autocovariance vector of the signal, up to a given lag
// [[Rcpp::export]]
double get_cov_W_scale_1_from_autocov_cpp_with_treshold(int h, const arma::vec autocov_vec, int lag_treshold=-1) {

  // h is the lag at which we want to compute the covariance of wavelet coefficient
  // autocov_vec is a vector of the autocivariance of x from lag 0 to lag
  // lag_treshold is an optimal parameter that allows to ignore the

  int j = 1;
  int L_j = 2;
  double out = 0;
  if(lag_treshold != -1 && lag_treshold <= 0){
    Rcpp::stop("lag_treshold should be a positive value greater than 0.");
 // only compute out if lag is greater than h, else return out = 0
  }else if(lag_treshold == -1 || lag_treshold >= h+1){
    double tmp = 2 * autocov_vec(std::abs(h)) - autocov_vec(std::abs(h + 1)) - autocov_vec(std::abs(h - 1));
    out = 1 / std::pow(2, 2 * j) * (L_j / 2) * tmp;
  }
  return(out);
}



// return a matrix with the autocovariance of the wavelet coefficient from 0 to n by column for each scale j, from the autocovariance vector of the signal
// [[Rcpp::export]]
arma::mat compute_all_cov_wv_recursive_2_cpp_with_mat(int n, const arma::vec autocov_vec, int lag_treshold = -1) {
  int max_j = floor(log2(n)) - 1;
  double sum_power_2_max_j_minus_1 = sum_of_powers_of_2(1, max_j - 1);

  if (autocov_vec.n_elem < n + sum_power_2_max_j_minus_1 + 2) {
    stop("autocov vec is not large enough, it should be at least of length " + std::to_string(n + sum_power_2_max_j_minus_1 + 2));
  }
  // max h to compute for scale 1
  int max_h_to_compute_scale_1 = n + sum_power_2_max_j_minus_1;
  // create arma mat, save by column the autocovariance of the wavelet coefficient from 0 to n for each scale j
  arma::mat mat_cov_W(n+1, max_j, arma::fill::zeros);
  // int j = 1;
  // int L_j = pow(2, j);
  // int M_j = n - L_j + 1;
  // int M_j1 = n - pow(2, (j + 1)) + 1;
  arma::vec cov_W_vec_lag(max_h_to_compute_scale_1 + 1, arma::fill::zeros);
  for (int t = 0; t <= max_h_to_compute_scale_1; ++t) {
    cov_W_vec_lag(t) = get_cov_W_scale_1_from_autocov_cpp_with_treshold(t,autocov_vec, lag_treshold);
  }
  // assign to first column of mat_cov_W
  mat_cov_W.col(0) = cov_W_vec_lag.subvec(0, n);
  // compute for all other scales
  for (int j = 2; j <= max_j; ++j) {
    // int L_j = pow(2, j);
    // int M_j = n - L_j + 1;
    // int M_j1 = n - pow(2, (j + 1)) + 1;
    int max_h_to_compute = max_h_to_compute_scale_1 - sum_of_powers_of_2(1, (j-1));
    int L_j_previous_scale =  pow(2, (j-1));
    // save cov_W at previous scale
    arma::vec tmp = cov_W_vec_lag;
      for (int t = 0; t <= max_h_to_compute; ++t) {
        cov_W_vec_lag[t] =  1.5*tmp(t) + tmp(abs(t+L_j_previous_scale/2)) + tmp(abs(t-L_j_previous_scale/2)) + 0.25 * tmp(abs(t+L_j_previous_scale)) + 0.25 * tmp(abs(t-L_j_previous_scale));
      }
    // save in matrix
    mat_cov_W.col(j-1) = cov_W_vec_lag.subvec(0, n);
    }
  return mat_cov_W;
}

// function that computes the autocovariance of the wavelet coefficients at scale j = 1 only.
// Return a vector of length n + sum_power_2_max_j_minus_1
// [[Rcpp::export]]
arma::vec compute_autocov_W_j_equal_1_from_autocov_X(const arma::vec autocov_vec, int n, int lag_treshold = -1) {
  int max_j = floor(log2(n)) - 1;
  double sum_power_2_max_j_minus_1 = sum_of_powers_of_2(1, max_j - 1);

  if (autocov_vec.n_elem < n + sum_power_2_max_j_minus_1 + 2) {
    stop("autocov vec is not large enough, it should be at least of length " + std::to_string(n + sum_power_2_max_j_minus_1 + 2));
  }
  // max h to compute for scale 1
  int max_h_to_compute_scale_1 = n + sum_power_2_max_j_minus_1;
  // create arma mat, save by column the autocovariance of the wavelet coefficient from 0 to n for each scale j
  arma::mat mat_cov_W(n+1, max_j, arma::fill::zeros);
  // int j = 1;
  // int L_j = pow(2, j);
  // int M_j = n - L_j + 1;
  // int M_j1 = n - pow(2, (j + 1)) + 1;
  arma::vec cov_W_vec_lag(max_h_to_compute_scale_1 + 1, arma::fill::zeros);
  for (int t = 0; t <= max_h_to_compute_scale_1; ++t) {
    cov_W_vec_lag(t) = get_cov_W_scale_1_from_autocov_cpp_with_treshold(t,autocov_vec, lag_treshold);
  }
  return(cov_W_vec_lag);

}


// [[Rcpp::export]]
arma::mat compute_all_cov_W_recursive_from_j_2(int n, arma::vec autocov_W_j_equal_1) {
  int max_j = floor(log2(n)) - 1;
  double sum_power_2_max_j_minus_1 = sum_of_powers_of_2(1, max_j - 1);
  // max h to compute for scale 1
  int max_h_to_compute_scale_1 = n + sum_power_2_max_j_minus_1;
  if (static_cast<int>(autocov_W_j_equal_1.n_elem) < max_h_to_compute_scale_1) {
    stop("autocov vector of the wavelet coefficient is not large enough, it should be at least of length " + std::to_string(max_h_to_compute_scale_1));
  }
  // create arma mat, save by column the autocovariance of the wavelet coefficient from 0 to n for each scale j
  arma::mat mat_cov_W(n+1, max_j, arma::fill::zeros);
  // assign to first column of mat_cov_W
  mat_cov_W.col(0) = autocov_W_j_equal_1.subvec(0, n);
  //  save autocov of W at scale 1 in temporary vector that we upadte at each scale
  arma::vec cov_W_vec_lag = autocov_W_j_equal_1;
  // compute for all other scales
  for (int j = 2; j <= max_j; ++j) {
    // int L_j = pow(2, j);
    // int M_j = n - L_j + 1;
    // int M_j1 = n - pow(2, (j + 1)) + 1;
    int max_h_to_compute = max_h_to_compute_scale_1 - sum_of_powers_of_2(1, (j-1));
    int L_j_previous_scale =  pow(2, (j-1));
    // save cov_W at previous scale
    arma::vec tmp = cov_W_vec_lag;
    for (int t = 0; t <= max_h_to_compute; ++t) {
      cov_W_vec_lag[t] =  1.5*tmp(t) + tmp(abs(t+L_j_previous_scale/2)) + tmp(abs(t-L_j_previous_scale/2)) + 0.25 * tmp(abs(t+L_j_previous_scale)) + 0.25 * tmp(abs(t-L_j_previous_scale));
    }
    // save in matrix
    mat_cov_W.col(j-1) = cov_W_vec_lag.subvec(0, n);
  }

  return mat_cov_W;
}



// [[Rcpp::export]]
double get_cov_wvar_cpp(int j, int k, int n, arma::mat mat_autocov_W){
  int J = std::floor(std::log2(n)) - 1;
  if (j + k > J) {
    stop("This function requires j+k <= J, i.e., not exceed the maximum scale of wavelet!");
  }

  // define constants
  double pow2j = std::pow(2, j);
  double pow2k = std::pow(2, k);
  double pow2jk = std::pow(2, j + k);
  double M_j = n - pow2j + 1;
  double M_jk = n - pow2jk + 1;
  double pow2jk_minus_1 = std::pow(2, j + k - 1);
  double pow2j_minus_1 = std::pow(2, j - 1);

  // int index_max = std::max(M_jk - 1 + (pow2k - 2) * pow2j_minus_1 + pow2jk_minus_1,
  //                          M_jk - 1 + pow2jk_minus_1 - pow2j_minus_1);
  // int index_min = std::min(-M_j + 1, -M_j + 1 + pow2jk_minus_1 - pow2j_minus_1);

  // Extract the vector of autocov of the wavelet coefficients from the Rcpp list
  arma::vec cov_W_vec = mat_autocov_W.col(j-1);

  double out = 0;
  for (int t = -M_j + 1; t <= M_jk - 1; ++t) {
    double C_lm = 0;
    for (int p = 0; p <= pow2k - 2; ++p) {
      C_lm += (p + 1) / pow2k * cov_W_vec(std::abs(t + p * pow2j_minus_1))
      + (pow2k - 1 - p) / pow2k * cov_W_vec(std::abs(t + p * pow2j_minus_1 + pow2jk_minus_1));
    }
    C_lm += cov_W_vec(std::abs(t + pow2jk_minus_1 - pow2j_minus_1));

    double tmp = M_jk * (t <= 0 && t >= (pow2j - pow2jk))
      + (M_jk - std::abs(t)) * (t > 0) + (M_jk - std::abs(t - pow2j + pow2jk)) * (t < pow2j - pow2jk);
    out += tmp * C_lm * C_lm;
  }
  return 2.0 / (M_j * M_jk) * out;
}


// Function that computes the variance of the wavelet variance at scale j from the autocovariance of the wavelet coefficients at scale j
// [[Rcpp::export]]
double get_var_wvar_j_from_autcov_W_j_cpp(int j, int n, arma::vec autocov_W) {
  int L_j = pow(2, j);
  int M_j = n - L_j + 1;
  double out = 0;

  for (int i = -(M_j - 1); i < M_j; ++i) {
    out += (M_j - std::abs(i)) * pow(autocov_W[std::abs(i)], 2);
  }

  return 2.0 / pow(M_j, 2) * out;
}

//  function to compute the theoretical variance covariance matrix of the wavelet variance based on either the autocovariance or the process X or the autcovariance of the wavelet coefficent at scale 1
// [[Rcpp::export]]
arma::mat get_theo_cov_matrix_wvar_cpp(int n, Rcpp::Nullable<arma::vec> autocov_vec_X = R_NilValue, Rcpp::Nullable<arma::vec> autocov_vec_W = R_NilValue, int num_off_diagonal = -1, int lag_treshold = -1){

  // compute max J
  int max_j = std::floor(std::log2(n)) - 1;

  // check number of diagonal parameters, if equal to -1, set to max j
  if (num_off_diagonal == -1) {
    num_off_diagonal = max_j; // Set it to max_j by default
  }
  //  make sure that either autocov_vec_X or autocov_vec_W is provided
  if (autocov_vec_X.isNull() && autocov_vec_W.isNull()) {
    Rcpp::stop("Provide either the vector of autocovariance of X or the vector of autocovariance of the Wavelet coefficient at scale j=1.");
  }
  // Throw an error if both autocov_vec_X and autocov_vec_W are provided
  if (autocov_vec_X.isNotNull() && autocov_vec_W.isNotNull()) {
    Rcpp::stop("Error, both vector are provided. Provide either the vector of autocovariance of X or the vector of autocovariance of the Wavelet coefficient at scale j=1.");
  }
  //  create arma mat of n+ 1 by max j to save the vector of autocovariance of W per scale j
  arma::mat mat_autocov_W(n+1, max_j, arma::fill::zeros);

  // check length of autocovariance vector of X and use it directly if not null to compute the autocvariance of the Wavelet coefficents
  if(autocov_vec_X.isNotNull()){
    arma::vec autocov_vec_X2 = Rcpp::as<arma::vec>(autocov_vec_X);
    if (static_cast<int>(autocov_vec_X2.n_elem) < n + sum_of_powers_of_2(1, max_j - 1) + 1) {
      stop("autocov vec is not large enough, it should be at least of length " + std::to_string(n + sum_of_powers_of_2(1, max_j - 1) + 1));
    }
    mat_autocov_W = compute_all_cov_wv_recursive_2_cpp_with_mat(n, autocov_vec_X2, lag_treshold);
  }
  //  Use autocovariance of the Wavelet coeffient if not null
  if(autocov_vec_W.isNotNull()){
    arma::vec autocov_vec_W2 = Rcpp::as<arma::vec>(autocov_vec_W);
    // check length of autocovariance vector of W
    if (static_cast<int>(autocov_vec_W2.n_elem) < n + sum_of_powers_of_2(1, max_j - 1) + 1) {
      stop("autocov vec of the Wavelet coefficient is not large enough, it should be at least of length " + std::to_string(n + sum_of_powers_of_2(1, max_j - 1)+ 1));
    }
    //  compute autocovariance of W
    mat_autocov_W = compute_all_cov_W_recursive_from_j_2(n,  autocov_vec_W2);
  }
  //  create armadillo vector of length max j
  arma::vec var_wv(max_j);

  // compute variance of wavelet variance for all scales based on the autocov of the wavelet coefficient at scale j
  for(int i = 1; i <= max_j; i++){
    var_wv(i-1)  = get_var_wvar_j_from_autcov_W_j_cpp(i, n, mat_autocov_W.col(i-1));
  }

  // create variance covariance matrix
  arma::mat cov_mat_wvar(max_j, max_j);

  //  assign variance of wavelet variance as diagonal of the matrix
  cov_mat_wvar.diag() = var_wv;

  // Calculate off-diagonal elements of the covariance matrix
  for (int j = 0; j < max_j - 1; ++j) {
    for (int k = j + 1; k < std::min(j + num_off_diagonal + 1, max_j); ++k) {
      cov_mat_wvar(j, k) = cov_mat_wvar(k, j) = get_cov_wvar_cpp(j + 1, k - j, n, mat_autocov_W );
    }
  }

  return cov_mat_wvar;
}

