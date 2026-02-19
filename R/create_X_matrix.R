create_X_matrix <- function(all_mjd_index,
                            jumps,
                            n_seasonal,
                            vec_earthquakes_index_mjd,
                            vec_earthquakes_relaxation_time) {





  # x = download_station_ngl("CHML")
  # # create all jumps by combining jumps due to equipment change and jumps due to earthquakes
  # jumps = c(x$df_equipment_software_changes$modified_julian_date,
  #           x$df_earthquakes$modified_julian_date )
  #
  # # if multiple jumps  to prevent not invertible matrix
  # jumps = unique(jumps)
  # vec_earthquakes_index_mjd = c(x$df_earthquakes$modified_julian_date)
  # # if multiple earthquakes  to prevent not invertible matrix
  # vec_earthquakes_index_mjd = unique(vec_earthquakes_index_mjd)
  # vec_earthquakes_relaxation_time = NULL
  #
  # all_mjd_index = seq(head(x$df_position$modified_julian_day, 1), tail(x$df_position$modified_julian_day, 1))



  # ensure that vec_earthquakes_relaxation_time is the same length as vec_earthquakes_index_mjd if not NULL
  if (!is.null(vec_earthquakes_relaxation_time)) {
    if (length(vec_earthquakes_index_mjd) != length(vec_earthquakes_relaxation_time)) {
      stop("Vector of relaxation time not of same length as vector of MJD index of earthquakes")
    }
  }


  # set to default 365.25 days if no earthquakes relaxation time provided
  if (is.null(vec_earthquakes_relaxation_time)) {
    vec_earthquakes_relaxation_time <- rep(365.25, length(vec_earthquakes_index_mjd))
  }


  # create empty matrix number of columns bias+trend+2*nbr sinusoidal + nbr jumps + number earthquakes
  X <- matrix(0, nrow = length(all_mjd_index), ncol = 2 + 2 * n_seasonal +
                length(jumps) + length(vec_earthquakes_index_mjd))

  # add bias, intercept
  X[, 1] <- 1

  # add component for trend and scale with respect to middle of time axis
  reference_time <- 0.5 * (all_mjd_index[1] + tail(all_mjd_index, 1))
  X[, 2] <- all_mjd_index - reference_time

  # add seasonal
  if (n_seasonal > 0) {
    for (i in 1:n_seasonal) {
      omega_i <- (i / 365.25) * 2 * pi
      X[, 2 + (i - 1) * 2 + 1] <- sin((all_mjd_index) * omega_i)
      X[, 2 + (i - 1) * 2 + 2] <- cos((all_mjd_index) * omega_i)
    }
  }

  # add offsets
  if (!is.null(jumps)) {
    for (i in 1:length(jumps)) {
      it <- min(which(all_mjd_index > jumps[i] - 1e-06))
      X[, 2 + 2 * n_seasonal + i] <- c(rep(0, it - 1), rep(
        1,
        length(all_mjd_index) - it + 1
      ))
    }
  }

  # exponential decay function for post seismic relaxation
  if (!is.null(vec_earthquakes_index_mjd)) {
    for (i in seq_along(vec_earthquakes_index_mjd)) {
      tau_i <- vec_earthquakes_relaxation_time[i]
      earthquake_mjd_i <- vec_earthquakes_index_mjd[i]
      # create vector
      decay_values <- ifelse(all_mjd_index >= earthquake_mjd_i,
                             1 - exp(-(all_mjd_index - earthquake_mjd_i) / tau_i),
                             0
      )

      # create column in matrix
      X[, 2 + 2 * n_seasonal + length(jumps) + i] <- decay_values
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
  rownames(X) <- all_mjd_index
  return(X)
}

