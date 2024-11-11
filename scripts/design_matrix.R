create_X_matrix = function (all_mjd_index,
                            jumps=NULL,
                            n_seasonal=2,
                            vec_earthquakes_index_mjd,
                            vec_earthquakes_relaxation_time) {


  # ensure that vec_earthquakes_relaxation_time is the same length as vec_earthquakes_index_mjd if not NULL
  if(!is.null(vec_earthquakes_relaxation_time)){
    if(length(vec_earthquakes_index_mjd) != length(vec_earthquakes_relaxation_time)){
      stop("Vector of relaxation time not of same length as vector of MJD index of earthquakes")
    }
  }


  # set to default 365.25 days if no earthquakes relaxation time provided
  if(is.null(vec_earthquakes_relaxation_time)){
    vec_earthquakes_relaxation_time = rep(365.25, length(vec_earthquakes_index_mjd))
  }


  # create empty matrix number of columns bias+trend+2*nbr sinusoidal + nbr jumps + number earthquakes
  X = matrix(0, nrow = length(all_mjd_index), ncol = 2 + 2 * n_seasonal +
               length(jumps) + length(vec_earthquakes_index_mjd))

  # add bias, intercept
  X[, 1] = 1

  # add component for trend
  reference_time =  0.5*(all_mjd_index[1]+tail(all_mjd_index,1))
  X[,2] = all_mjd_index - reference_time


  # add seasonal
  if (n_seasonal > 0) {
    for (i in 1:n_seasonal) {
      X[, 2 + (i - 1) * 2 + 1] = sin((all_mjd_index) *  i * 2 * pi/365.25)
      X[, 2 + (i - 1) * 2 + 2] = cos((all_mjd_index) * i * 2 * pi/365.25)
    }
  }

  # add offsets
  if (!is.null(jumps)) {
    for (i in 1:length(jumps)) {
      it = min(which(all_mjd_index > jumps[i] - 1e-06))
      X[, 2 + 2 * n_seasonal + i] = c(rep(0, it - 1), rep(1,
                                                          length(all_mjd_index) - it + 1))
    }
  }


  # exponential decay function for post seismic relaxation
  if(!is.null(vec_earthquakes_index_mjd)){
    for(i in seq_along(vec_earthquakes_index_mjd)){
      tau_i = vec_earthquakes_relaxation_time[i]
      earthquake_mjd_i = vec_earthquakes_index_mjd[i]
      # create vector
      decay_values <- ifelse(all_mjd_index > earthquake_mjd_i,
                             1 - exp(-(all_mjd_index - earthquake_mjd_i) / tau_i),
                             0)

      # create column in matrix
      X[, 2 + 2 * n_seasonal + length(jumps) + i] = decay_values
    }
  }


  # # Slow slip events (tanh)
  # if(!is.null(vec_tanh_mid_point)){
  #   for(i in seq_along(vec_tanh_mid_point)){
  #     # i = 1
  #     t_k = vec_tanh_mid_point[i]
  #     T_k = vec_tanh_length[i]
  #     X[, 2 + 2 * n_seasonal + length(jumps) + i] = 0.5 * (tanh((t_nogap - t_k)/T_k) -1)
  #   }
  # }
  rownames(X) = all_mjd_index
  return(X)
}



res=t(X) %*%X
solve(res)
Matrix::solve(res)
