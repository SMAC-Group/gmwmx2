#include <numeric>
//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;


// [[Rcpp::export]]
arma::mat compute_cov_W_scale_1_cpp(arma::mat Sigma_X){
  int n = Sigma_X.n_rows;
  int M = n-1;
  arma::mat Sigma_W_j_1(M, M, arma::fill::zeros);
  for(int k = 0; k < M; ++k){
    for(int l = 0 ; l < M; ++l){
      Sigma_W_j_1(k,l) =  0.25*Sigma_X(k,l) - 0.25*Sigma_X(k,l+1) - 0.25*Sigma_X(k+1,l) + 0.25*Sigma_X(k+1,l+1);
      // Sigma_W_j_1(k,l) =  1;

    }
  }
  return(Sigma_W_j_1);
}

// [[Rcpp::export]]
Rcpp::List compute_cov_W_all_scales_cpp(arma::mat Sigma_X){
  // get n
  int n = Sigma_X.n_rows;
  // get max J
  int max_j = std::floor(std::log2(n)) - 1;
  // create List
  Rcpp::List lst_cov_W(max_j);
  // Compute cov Wavelet coefficients at scale 1
  int M = n-1;
  arma::mat Sigma_W_j_1(M, M, arma::fill::zeros);
  for(int k = 0; k < M; ++k){
    for(int l = 0 ; l < M; ++l){
      Sigma_W_j_1(k,l) =  0.25*Sigma_X(k,l) - 0.25*Sigma_X(k,l+1) - 0.25*Sigma_X(k+1,l) + 0.25*Sigma_X(k+1,l+1);
    }
  }
  //  assign to first entry of list
  lst_cov_W(0) = Sigma_W_j_1;
  //  compute the other matrices of covariance of wavelet coefficients for all other scales
  for(int j = 1; j <= (max_j - 1); j++){
    //  extract previous matrix
    // arma::mat Sigma_W_j_minus_1 = Rcpp::as<arma::mat>(lst_cov_W(j-1));
    arma::mat Sigma_W_j_minus_1 = lst_cov_W(j-1);
    int L = std::pow(2, j);
    int M = n - std::pow(2, (j+1))+1;
    //  create matrix to save for covariance of wavelet coefficients at scale j
    arma::mat Sigma_W_j(M, M, arma::fill::zeros);
    for(int k = 0; k < M; ++k){
      for(int l = 0 ; l < M; ++l){
        double tmp1 = 0.25 * Sigma_W_j_minus_1(k,l) + 0.5*Sigma_W_j_minus_1(k,l+L/2) + 0.25 * Sigma_W_j_minus_1(k,l+L);
        double tmp2 = 0.5*Sigma_W_j_minus_1(k+L/2,l) + Sigma_W_j_minus_1(k+L/2,l+L/2) + 0.5*Sigma_W_j_minus_1(k+L/2,l+L);
        double tmp3 = 0.25*Sigma_W_j_minus_1(k+L,l) + 0.5*Sigma_W_j_minus_1(k+L,l+L/2) + 0.25*Sigma_W_j_minus_1(k+L,l+L);
        Sigma_W_j(k,l) = tmp1 + tmp2 + tmp3;
      }
    }
    //  assign to entry j
    lst_cov_W(j) = Sigma_W_j;
  }

  return(lst_cov_W);
}



// Function to compute all pairs of scales

// [[Rcpp::export]]
arma::mat compute_index_pairs_scales_cpp(int max_j){
  //  get number of combination
  int num_combinations= (max_j * (max_j + 1)) / 2;
  // create matrix
  arma::mat mat_index(num_combinations, 2, arma::fill::zeros);
  int counter = 0;
  for (int i = 1; i <= max_j; i++) {
    for (int j = i; j <= max_j; j++) {
      mat_index(counter, 0) =i;
      mat_index(counter, 1) =j;
      // update counter
      counter += 1;
    }
  }
 return(mat_index) ;
}



// function to compute all covariance matrix between wavelet coefficients at different scales
// [[Rcpp::export]]
Rcpp::List compute_cov_W_between_scales_all_scales_cpp(arma::mat mat_index, int n, Rcpp::List lst_cov_W){
  int nrow_mat_index = mat_index.n_rows;
  // create List
  Rcpp::List lst_cov_W_j_k(nrow_mat_index);
  for (int count = 0; count < nrow_mat_index; count++) {
    int j = mat_index(count, 0);
    int k = mat_index(count, 1) - mat_index(count, 0);
    int Mj = n - std::pow(2, j) +1;
    int Mjk = n - std::pow(2, j+k)+1;
    int Cp_size = std::pow(2, (k+1))-1;
    // create vector Cp
    arma::vec Cp(Cp_size,  arma::fill::zeros);
    // fill in Cp
    for (int p = 0; p <= Cp_size - 1; p++) {
      if (p >= 0 && p <= std::pow(2, k) - 1) {
        Cp(p) = (p + 1) / std::pow(2, k);
      } else if (p >= std::pow(2, k) && p <= std::pow(2, k + 1) - 2) {
        Cp(p)  = (std::pow(2, k + 1) - 1 - p) / std::pow(2, k);
      }
    }
    // create matrix
    arma::mat cov_W_j_k(Mjk, Mj, arma::fill::zeros);
    //  extract Sigma_W_j mat
    arma::mat Sigma_Wj = lst_cov_W(j-1);
    //  fill matrix
    for (int l = 0; l <= (Mjk-1) ; l++) {
      for(int m = 0; m <= (Mj-1); m++){
          for (int p = 0; p <= (std::pow(2, (k+1)) -2); p++) {
            cov_W_j_k(l,m) = cov_W_j_k(l,m) + Cp(p)* Sigma_Wj(l+p*std::pow(2, (j-1)), m);
          }

      }
    }

    // append covariance matrix to lst_cov_W_j_k
    lst_cov_W_j_k(count) = cov_W_j_k;
  }
  return(lst_cov_W_j_k);
}

// function f_jk
// [[Rcpp::export]]
double f_jk_cpp(int h, arma::mat cov_Wj_Wjk, int Mjk, int Mj){
  double out = 0.0;

  if (h >= Mjk - Mj && h <= 0) {
    for (int l = 0; l < Mjk; l++) {
      out += (1.0 / Mjk) * std::pow(cov_Wj_Wjk(l, l - h), 2);
    }
    return std::sqrt(out);
  }

  if (h >= 1 && h <= Mjk - 1) {
    for (int l = h; l < Mjk; l++) {
      out += (1.0 / (Mjk - h)) * std::pow(cov_Wj_Wjk(l, l - h), 2);
    }
    return std::sqrt(out);
  }

  if (h >= 1 - Mj && h <= Mjk - Mj - 1) {
    int limit = Mjk - std::abs(h - Mjk + Mj);
    for (int l = 0; l < limit; l++) {
      out += (1.0 / limit) * std::pow(cov_Wj_Wjk(l, l - h), 2);
    }
    return std::sqrt(out);
  }else{
    return(0);
  }
}




// Create f_j_k approx
// [[Rcpp::export]]
double f_jk_approx_cpp(int h, arma::mat cov_Wj_Wjk, int Mjk, int Mj, std::string approx_type="3"){
  double out = 0.0;

  if (h >= Mjk - Mj && h <= 0) {

    // create full sequence of l
    arma::vec vec_full_l = arma::linspace(0, Mjk-1);
    int num_total_points = vec_full_l.n_elem;
    int num_points_on_which_to_compute = std::floor(std::log2(num_total_points)) ;

    // Create a sequence from 0 to Mjk-1 with length `num_points_on_which_to_compute`
    arma::vec myseq = arma::round(arma::linspace(0, Mjk-1, num_points_on_which_to_compute));
    arma::uvec indices_l_values = arma::conv_to<arma::uvec>::from(myseq);

    int num_elements = static_cast<int>(indices_l_values.n_elem);

    for (int l_index = 0; l_index < num_elements; l_index++) {
      int l = indices_l_values(l_index);
      out += (1.0 / num_points_on_which_to_compute) * std::pow(cov_Wj_Wjk(l, l - h), 2);
    }
    return std::sqrt(out);
  }

  if (h >= 1 && h <= Mjk - 1) {

    // create full sequence of l
    arma::vec vec_full_l = arma::linspace(h, Mjk-1);
    int num_total_points = vec_full_l.n_elem;
    int num_points_on_which_to_compute = std::floor(std::log2(num_total_points)) ;

    // Create a sequence from 0 to Mjk-1 with length `num_points_on_which_to_compute`
    arma::vec myseq = arma::round(arma::linspace(h, Mjk-1, num_points_on_which_to_compute));
    arma::uvec indices_l_values = arma::conv_to<arma::uvec>::from(myseq);
    int num_elements = static_cast<int>(indices_l_values.n_elem);


    for (int l_index = 0; l_index < num_elements; l_index++) {
      int l =  indices_l_values(l_index);
      out += (1.0 / num_points_on_which_to_compute) * std::pow(cov_Wj_Wjk(l, l - h), 2);
    }
    return std::sqrt(out);
  }

  if (h >= 1 - Mj && h <= Mjk - Mj - 1) {
    int limit = Mjk - std::abs(h - Mjk + Mj);

    // create full sequence of l
    arma::vec vec_full_l = arma::linspace(0, limit-1);
    int num_total_points = vec_full_l.n_elem;
    int num_points_on_which_to_compute = std::floor(std::log2(num_total_points)) ;

    // Create a sequence from 0 to Mjk-1 with length `num_points_on_which_to_compute`
    arma::vec myseq = arma::round(arma::linspace(0, limit-1, num_points_on_which_to_compute));
    arma::uvec indices_l_values = arma::conv_to<arma::uvec>::from(myseq);
    int num_elements = static_cast<int>(indices_l_values.n_elem);
    for (int l_index = 0; l_index < num_elements; l_index++) {
      int l = indices_l_values(l_index);
      out += (1.0 / num_points_on_which_to_compute) * std::pow(cov_Wj_Wjk(l, l - h), 2);
    }
    return std::sqrt(out);
  }else{
    return(0);
  }
}



// function compute_cov_wv_cpp_1
// [[Rcpp::export]]
arma::mat compute_cov_wv_cpp_1(Rcpp::List cov_W_j, int J, int n, arma::mat mat_index){
  //  create matrix
  arma::mat cov_wv_res(J, J, arma::fill::zeros);
  for (int j = 1; j <= J; j++) {
    for (int k1 = j; k1 <= J; k1++) {
      int k = k1 - j;
      int Mj = n - std::pow(2, j) + 1;
      int Mjk = n - std::pow(2, (j + k)) + 1;
      // Find the first index where index(:,1) == j and index(:,2) == k1
      int my_index = -1;  // Initialize to an invalid index
      for (arma::uword i = 0; i < mat_index.n_rows; ++i) {
        if (mat_index(i, 0) == j && mat_index(i, 1) == k1) {
          my_index = i;
          break;  // Exit the loop after finding the first match
        }
      }
      // // Check if a valid index was found
      // if (my_index == -1) {
      //   Rcpp::Rcout << "Warning: No valid index found for j = " << j << " and k1 = " << k1 << std::endl;
      //   continue;  // Skip if no matching index was found
      // }


      // Print the current index and values of j, k1, and my_index
      // Rcpp::Rcout << "j = " << j << ", k1 = " << k1 << ", my_index = " << my_index << std::endl;


      // Create the sequence my_h from 1-Mj to Mjk-1
      std::vector<int> my_h;
      for (int h = 1 - Mj; h <= Mjk - 1; ++h) {
        my_h.push_back(h);
      }

      // Print the contents of my_h
      // Rcpp::Rcout << "my_h = ";
      // for (std::size_t i = 0; i < my_h.size(); ++i) {
      //   Rcpp::Rcout << my_h[i] << " ";
      // }
      // Rcpp::Rcout << std::endl;



      // initialize out
      double out = 0.0;
      for (arma::uword i = 0; i < my_h.size(); ++i) {
        int h = my_h[i];
        //  extract matrix
        arma::mat cov_Wj_Wjk = Rcpp::as<arma::mat>(cov_W_j[my_index]);

        // // Debug: Print the matrix dimensions
        // Rcpp::Rcout << "cov_Wj_Wjk (" << my_index << ") dimensions: "
        //             << cov_Wj_Wjk.n_rows << " x " << cov_Wj_Wjk.n_cols << std::endl;

        // // Print the first row of cov_Wj_Wjk
        // Rcpp::Rcout << "First row of cov_Wj_Wjk: ";
        // for (arma::uword col = 0; col < cov_Wj_Wjk.n_cols; ++col) {
        //   Rcpp::Rcout << cov_Wj_Wjk(0, col) << " "; // Print each element in the first row
        // }
        // Rcpp::Rcout << std::endl;


        double my_f_jk = f_jk_cpp(h, cov_Wj_Wjk, Mjk, Mj);

        int ch = Mjk * (h >= Mjk - Mj && h <= 0) +
          (Mjk - h) * (h >= 1 && h <= Mjk - 1) +
          (Mjk - std::abs(h - Mjk + Mj)) * (h >= 1 - Mj && h <= Mjk - Mj - 1);

        out += ch * std::pow(my_f_jk, 2);
      }

      double value = 2.0 / (Mj * Mjk) * out;
      cov_wv_res(j-1, k1-1) = value;
      cov_wv_res(k1-1, j-1) = value;  // Fill in the symmetric part
    }
  }
  return cov_wv_res;
}


// write function to compute covariance of wavelet variance directly
// [[Rcpp::export]]
arma::mat compute_cov_wv_cpp(arma::mat Sigma_X){

  // --------------------------------- get variance covariance matrix of the wavelet coefficients for all scales
  // get n
  int n = Sigma_X.n_rows;
  // get max J
  int max_j = std::floor(std::log2(n)) - 1;
  // create List containing the variance covariance matrix of the wavelet coefficients for all scales
  Rcpp::List lst_cov_W(max_j);
  // Compute cov Wavelet coefficients at scale 1
  int M = n-1;
  arma::mat Sigma_W_j_1(M, M, arma::fill::zeros);
  for(int k = 0; k < M; ++k){
    for(int l = 0 ; l < M; ++l){
      Sigma_W_j_1(k,l) =  0.25*Sigma_X(k,l) - 0.25*Sigma_X(k,l+1) - 0.25*Sigma_X(k+1,l) + 0.25*Sigma_X(k+1,l+1);
    }
  }
  //  assign to first entry of list
  lst_cov_W(0) = Sigma_W_j_1;
  //  compute the other matrices of covariance of wavelet coefficients for all other scales
  for(int j = 1; j <= (max_j - 1); j++){
    //  extract previous matrix
    // arma::mat Sigma_W_j_minus_1 = Rcpp::as<arma::mat>(lst_cov_W(j-1));
    arma::mat Sigma_W_j_minus_1 = lst_cov_W(j-1);
    int L = std::pow(2, j);
    int M = n - std::pow(2, (j+1))+1;
    //  create matrix to save for covariance of wavelet coefficients at scale j
    arma::mat Sigma_W_j(M, M, arma::fill::zeros);
    for(int k = 0; k < M; ++k){
      for(int l = 0 ; l < M; ++l){
        double tmp1 = 0.25 * Sigma_W_j_minus_1(k,l) + 0.5*Sigma_W_j_minus_1(k,l+L/2) + 0.25 * Sigma_W_j_minus_1(k,l+L);
        double tmp2 = 0.5*Sigma_W_j_minus_1(k+L/2,l) + Sigma_W_j_minus_1(k+L/2,l+L/2) + 0.5*Sigma_W_j_minus_1(k+L/2,l+L);
        double tmp3 = 0.25*Sigma_W_j_minus_1(k+L,l) + 0.5*Sigma_W_j_minus_1(k+L,l+L/2) + 0.25*Sigma_W_j_minus_1(k+L,l+L);
        Sigma_W_j(k,l) = tmp1 + tmp2 + tmp3;
      }
    }
    //  assign to entry j
    lst_cov_W(j) = Sigma_W_j;
  }

  // -------------------- compute index mat that will save all combinations between scales

  //  get number of combination
  int num_combinations= (max_j * (max_j + 1)) / 2;
  // create matrix
  arma::mat mat_index(num_combinations, 2, arma::fill::zeros);
  int counter = 0;
  for (int i = 1; i <= max_j; i++) {
    for (int j = i; j <= max_j; j++) {
      mat_index(counter, 0) =i;
      mat_index(counter, 1) =j;
      // update counter
      counter += 1;
    }
  }

  //  -------------------------- compute the covariance between wavelet coefficients of different scales
  int nrow_mat_index = mat_index.n_rows;
  // create List
  Rcpp::List lst_cov_W_j_k(nrow_mat_index);
  for (int count = 0; count < nrow_mat_index; count++) {
    int j = mat_index(count, 0);
    int k = mat_index(count, 1) - mat_index(count, 0);
    // do not recompute covariance matrix of wavelet coefficients at a certain scale if its cov between W_j1 and W_j2 with j1 = j2
    // if(k == 0){
    //   lst_cov_W_j_k(count) = lst_cov_W(j-1);
    //   continue;
    // }
    int Mj = n - std::pow(2, j) +1;
    int Mjk = n - std::pow(2, j+k)+1;
    int Cp_size = std::pow(2, (k+1))-1;
    // create vector Cp
    arma::vec Cp(Cp_size,  arma::fill::zeros);
    // fill in Cp
    for (int p = 0; p <= Cp_size - 1; p++) {
      if (p >= 0 && p <= std::pow(2, k) - 1) {
        Cp(p) = (p + 1) / std::pow(2, k);
      } else if (p >= std::pow(2, k) && p <= std::pow(2, k + 1) - 2) {
        Cp(p)  = (std::pow(2, k + 1) - 1 - p) / std::pow(2, k);
      }
    }
    // create matrix
    arma::mat cov_W_j_k(Mjk, Mj, arma::fill::zeros);
    //  extract Sigma_W_j mat
    arma::mat Sigma_Wj = lst_cov_W(j-1);
    //  fill matrix
    for (int l = 0; l <= (Mjk-1) ; l++) {
      for(int m = 0; m <= (Mj-1); m++){
        for (int p = 0; p <= (std::pow(2, (k+1)) -2); p++) {
          cov_W_j_k(l,m) = cov_W_j_k(l,m) + Cp(p)* Sigma_Wj(l+p*std::pow(2, (j-1)), m);
        }
      }
    }

    // append covariance matrix to lst_cov_W_j_k
    lst_cov_W_j_k(count) = cov_W_j_k;
  }

  // -------------------- compute covariance of wavelet variance
  //  create matrix
  arma::mat cov_wv_res(max_j, max_j, arma::fill::zeros);
  for (int j = 1; j <= max_j; j++) {
    for (int k1 = j; k1 <= max_j; k1++) {
      int k = k1 - j;
      int Mj = n - std::pow(2, j) + 1;
      int Mjk = n - std::pow(2, (j + k)) + 1;
      // Find the first index where index(:,1) == j and index(:,2) == k1
      int my_index = -1;  // Initialize to an invalid index
      for (arma::uword i = 0; i < mat_index.n_rows; ++i) {
        if (mat_index(i, 0) == j && mat_index(i, 1) == k1) {
          my_index = i;
          break;  // Exit the loop after finding the first match
        }
      }

      // Create the sequence my_h from 1-Mj to Mjk-1
      std::vector<int> my_h;
      for (int h = 1 - Mj; h <= Mjk - 1; ++h) {
        my_h.push_back(h);
      }


      // initialize out
      double out = 0.0;
      for (arma::uword i = 0; i < my_h.size(); ++i) {
        int h = my_h[i];
        //  extract matrix
        arma::mat cov_Wj_Wjk = Rcpp::as<arma::mat>(lst_cov_W_j_k[my_index]);

        double my_f_jk = f_jk_cpp(h, cov_Wj_Wjk, Mjk, Mj);

        int ch = Mjk * (h >= Mjk - Mj && h <= 0) +
          (Mjk - h) * (h >= 1 && h <= Mjk - 1) +
          (Mjk - std::abs(h - Mjk + Mj)) * (h >= 1 - Mj && h <= Mjk - Mj - 1);

        out += ch * std::pow(my_f_jk, 2);
      }

      double value = 2.0 / (Mj * Mjk) * out;
      cov_wv_res(j-1, k1-1) = value;
      cov_wv_res(k1-1, j-1) = value;  // Fill in the symmetric part
    }
  }
  return cov_wv_res;

}






// [[Rcpp::export]]
arma::uvec compute_l_index_to_compute_cpp(int Mjk, std::string approx_type) {
  int nbr_l_index_to_compute = 0;
  // Compute the number of indices to compute
  if(approx_type == "1"){
    nbr_l_index_to_compute = std::round(std::log10(Mjk)) + 1;
  }else if(approx_type == "2"){
    nbr_l_index_to_compute = std::round(std::log2(Mjk)) + 1;
  }else if(approx_type == "all"){
    return arma::regspace<arma::uvec>(1, Mjk);
  }

  // Determine the two boundary terms
  int term1 = nbr_l_index_to_compute;
  int term2 = Mjk - nbr_l_index_to_compute;

  // Generate the sequence between term1 + 1 and term2 - 1
  arma::uvec myseq = arma::conv_to<arma::uvec>::from(arma::linspace(term1 + 1, term2 - 1, nbr_l_index_to_compute));

  // Combine the different parts into one final vector
  arma::uvec l_index_to_compute(term1 + myseq.n_elem + (Mjk - term2 + 1));

  // First part: 1 to term1
  l_index_to_compute.head(term1) = arma::regspace<arma::uvec>(1, term1);

  // Middle part: myseq
  l_index_to_compute.subvec(term1, term1 + myseq.n_elem - 1) = myseq;

  // Last part: term2 to Mjk
  l_index_to_compute.tail(Mjk - term2 + 1) = arma::regspace<arma::uvec>(term2, Mjk);

  return l_index_to_compute;
}

// Function to find indices present in all_l_values
arma::uvec find_indices_present(const arma::uvec& all_l_values, const arma::uvec& lst_index_row) {
  // Initialize a vector to store indices of present values
  std::vector<arma::uword> indices_present;

  // Check each value in all_l_values
  for (arma::uword i = 0; i < all_l_values.n_elem; ++i) {
    // Use std::find to check if the current value exists in lst_index_row
    if (std::find(lst_index_row.begin(), lst_index_row.end(), all_l_values[i]) != lst_index_row.end()) {
      indices_present.push_back(i);  // Store index if present
    }
  }

  // Convert to Armadillo uvec for returning
  return arma::conv_to<arma::uvec>::from(indices_present);
}

// [[Rcpp::export]]
arma::uvec compute_indices(const arma::uvec& all_l_values, const arma::uvec& lst_index_row) {
  return find_indices_present(all_l_values, lst_index_row);
}



// write function to compute covariance of wavelet variance approximation
// [[Rcpp::export]]
arma::mat compute_cov_wv_cpp_approx_faster(arma::mat Sigma_X, std::string approx_type = "1"){

  // --------------------------------- get variance covariance matrix of the wavelet coefficients for all scales

  // get n
  int n = Sigma_X.n_rows;
  // get max J
  int max_j = std::floor(std::log2(n)) - 1;
  // create List containing the variance covariance matrix, for each j its a matrix Mj by Mj
  Rcpp::List lst_cov_W(max_j);
  // Compute number of wavelet coefficients at scale 1
  int M = n-1;
  //  create matrix of variance covariance of of wavelet coefficients at scale 1
  arma::mat Sigma_W_j_1(M, M, arma::fill::zeros);
  for(int k = 0; k < M; ++k){
    for(int l = 0 ; l < M; ++l){
      double val1 = Sigma_X(k,l);
      double val2 = Sigma_X(k,l+1);
      double val3 = Sigma_X(k+1,l);
      double val4 = Sigma_X(k+1,l+1);
      Sigma_W_j_1(k,l) =  0.25*(val1 -val2 - val3 + val4);
    }
  }
  //  assign this matrix to first entry of list
  lst_cov_W(0) = Sigma_W_j_1;

  //  compute the other matrices of covariance of wavelet coefficients for all other scales to 2 , ... J
  for(int j = 1; j <= (max_j - 1); j++){
    //  extract previous matrix and create a reference and define as const, do not perform a copy
    const arma::mat& Sigma_W_j_minus_1 = lst_cov_W[j - 1];
    // int L = std::pow(2, j);
    int L = 1 << j; // here we use  bitwise shifts rather than std::pow for faster computation
    // int M = n - std::pow(2, (j+1))+1;
    int M = n - (L << 1) + 1;
    //  create matrix to save for covariance of wavelet coefficients at scale j
    arma::mat Sigma_W_j(M, M, arma::fill::zeros);
    for(int k = 0; k < M; ++k){
      for(int l = 0 ; l < M; ++l){
        double val1 = Sigma_W_j_minus_1(k,l);
        double val2 = Sigma_W_j_minus_1(k,l+L/2);
        double val3 = Sigma_W_j_minus_1(k,l+L);
        double tmp1 = 0.25* (val1 + val3)  + 0.5* val2  ;
        double tmp2 = 0.5*Sigma_W_j_minus_1(k+L/2,l) + Sigma_W_j_minus_1(k+L/2,l+L/2) + 0.5*Sigma_W_j_minus_1(k+L/2,l+L);
        double tmp3 = 0.25*Sigma_W_j_minus_1(k+L,l) + 0.5*Sigma_W_j_minus_1(k+L,l+L/2) + 0.25*Sigma_W_j_minus_1(k+L,l+L);
        Sigma_W_j(k,l) = tmp1 + tmp2 + tmp3;
      }
    }
    //  assign to entry j
    lst_cov_W(j) = Sigma_W_j;
  }

  // -------------------- compute index mat that will save all combinations between scales

  //  get number of combination
  int num_combinations= (max_j * (max_j + 1)) / 2;
  // create matrix
  arma::mat mat_index(num_combinations, 2, arma::fill::zeros);
  int counter = 0;
  for (int i = 1; i <= max_j; i++) {
    for (int j = i; j <= max_j; j++) {
      mat_index(counter, 0) =i;
      mat_index(counter, 1) =j;
      // update counter
      counter += 1;
    }
  }

  //  -------------------------- compute the covariance between wavelet coefficients of different scales

  // get total number of combination
  int nrow_mat_index = mat_index.n_rows;
  //  create list for index
  Rcpp::List lst_index_jk(nrow_mat_index);
  // create List
  Rcpp::List lst_cov_W_j_k(nrow_mat_index);
  for (int count = 0; count < nrow_mat_index; count++) {
    int j = mat_index(count, 0);
    int k = mat_index(count, 1) - mat_index(count, 0);

    int Mj = n - (1 << j) +1;
    int Mjk = n - (1 << (j+k)) + 1;
    // int Cp_size = std::pow(2, (k+1))-1; std::pow(2, j) == 1 << j
    int Cp_size = (1 << (k+1)) - 1;
    // create vector Cp
    arma::vec Cp(Cp_size,  arma::fill::zeros);
    // fill in Cp
    double pow_2_k = 1<<k;
    double pow_2_k_1 = 1 << (k + 1);
    for (int p = 0; p <= Cp_size - 1; p++) {
      if (p >= 0 && p <= pow_2_k - 1) {
        Cp(p) = (p + 1) / pow_2_k;
      } else if (p >= pow_2_k && p <= pow_2_k_1 - 2) {
        Cp(p)  = (pow_2_k_1 - 1 - p) / pow_2_k;
      }
    }
    // starting approx
    arma::uvec l_index_to_compute = compute_l_index_to_compute_cpp(Mjk, approx_type);
    lst_index_jk[count] = l_index_to_compute ;

    // create matrix cov_W_j, j+k
    arma::mat cov_W_j_k(l_index_to_compute.size(), Mj, arma::fill::zeros);

    //  Create reference to variance covariance of wavelet coef at scale j-1
    const arma::mat& Sigma_Wj = lst_cov_W[j-1];
    //
    double val_2_pow_j_1 = std::pow(2, (j-1));
    //  fill matrix
    for (int l_index = 0; l_index < static_cast<int>(l_index_to_compute.size()) ; l_index++) {
      int l = l_index_to_compute[l_index]-1;
      for(int m = 0; m <= (Mj-1); m++){
        for (int p = 0; p <= (std::pow(2, (k+1)) -2); p++) {
          cov_W_j_k(l_index,m) += Cp(p)* Sigma_Wj(l+p*val_2_pow_j_1, m);
        }
      }
    }

    // append covariance matrix to lst_cov_W_j_k
    lst_cov_W_j_k(count) = cov_W_j_k;
  }

  // -------------------- compute covariance of wavelet variance

  //
  //  create matrix
  arma::mat cov_wv_res(max_j, max_j, arma::fill::zeros);
  for (int j = 1; j <= max_j; j++) {
    for (int k1 = j; k1 <= max_j; k1++) {
      int k = k1 - j;
      int Mj = n - (1 << j) +1;
      int Mjk = n - (1 << (j+k)) + 1;
      // int Mj = n - std::pow(2, j) + 1;
      // int Mjk = n - std::pow(2, (j + k)) + 1;
      // Find the first index where index(:,1) == j and index(:,2) == k1
      int my_index = -1;  // Initialize to an invalid index
      for (arma::uword i = 0; i < mat_index.n_rows; ++i) {
        if (mat_index(i, 0) == j && mat_index(i, 1) == k1) {
          my_index = i;
          break;  // Exit the loop after finding the first match
        }
      }

      // extract variance covariance of Wj, W_j+k at this index
      //  Create reference to variance covariance of wavelet coef at scale j-1
      const arma::mat& cov_W_j_k_my_index = lst_cov_W_j_k[my_index];

      // Create the sequence my_h from 1-Mj to Mjk-1
      std::vector<int> my_h(Mjk + Mj - 1);
      std::iota(my_h.begin(), my_h.end(), 1 - Mj);
      // std::vector<double> my_h;
      // my_h.reserve(Mjk + Mj - 1); // pre allocate memory
      // for (int h = 1 - Mj; h <= Mjk - 1; ++h) {
      //   my_h.push_back(h);
      // }

      //



      // initialize out
      double out = 0.0;
      for (arma::uword i = 0; i < my_h.size(); ++i) {
        int h = my_h[i];

        //  get all l value that we would compute normally
        arma::uvec all_l_values;

        if (h >= Mjk - Mj && h <= 0) {
          all_l_values = arma::regspace<arma::uvec>(1, Mjk);
        } else if (h >= 1 && h <= Mjk - 1) {
          all_l_values = arma::regspace<arma::uvec>(h + 1, Mjk);
        } else if (h >= 1 - Mj && h <= Mjk - Mj - 1) {
          all_l_values = arma::regspace<arma::uvec>(1, Mjk - std::abs(h - Mjk + Mj));
        }

        // find which value are found:
        // extract which are the row that are stored for this variance covariance (Wj, Wj+k)
        arma::uvec index_l_at_this_index = lst_index_jk[my_index];
        // find which of all the suposed row indices to sum are stored and their position in the vector all_l_values
        arma::uvec id_l_values_present = compute_indices(all_l_values, index_l_at_this_index );
        //  extract
        arma::uvec l_values_present = all_l_values.elem(id_l_values_present);
        // arma::uvec id_l_values_present = arma::find(arma::any(all_l_values == index_l_at_this_index.t() , 1));
        double inv_l_size = (1.0 / l_values_present.n_elem);
        double out_fjk = 0.0;
        for (int l_index = 0; l_index < static_cast<int>(l_values_present.size()); ++l_index) {
          int l = l_values_present[l_index];
          arma::uvec l_index_in_mat = arma::find(index_l_at_this_index == l);
          out_fjk += inv_l_size * std::pow(cov_W_j_k_my_index(l_index_in_mat(0), l-h-1), 2);
        }
        out_fjk = std::sqrt(out_fjk);
        // compute ch
        int ch;
        if (h >= Mjk - Mj && h <= 0) {
          ch = Mjk;
        } else if (h >= 1 && h <= Mjk - 1) {
          ch = Mjk - h;
        } else {
          ch = Mjk - std::abs(h - Mjk + Mj);
        }
        out += ch * std::pow(out_fjk, 2);

        // int ch = Mjk * (h >= Mjk - Mj && h <= 0) +
        //   (Mjk - h) * (h >= 1 && h <= Mjk - 1) +
        //   (Mjk - std::abs(h - Mjk + Mj)) * (h >= 1 - Mj && h <= Mjk - Mj - 1);

      }

      double value = 2.0 / (Mj * Mjk) * out;
      cov_wv_res(j-1, k1-1) = value;
      cov_wv_res(k1-1, j-1) = value;  // Fill in the symmetric part
    }
  }

  // return cov_wv_res;
  return(cov_wv_res);

}


