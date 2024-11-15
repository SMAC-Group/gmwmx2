\name{gmwmx2}
\alias{gmwmx2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{~~ Dummy title ~~}
\description{~~ Dummy description ~~}
\usage{
gmwmx2(x, n_seasonal = 2, vec_earthquakes_relaxation_time = NULL, component = "N", toeplitz_approx_var_cov_wv = TRUE, stochastic_model = "wn + fl")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{n_seasonal}{
%%     ~~Describe \code{n_seasonal} here~~
}
  \item{vec_earthquakes_relaxation_time}{
%%     ~~Describe \code{vec_earthquakes_relaxation_time} here~~
}
  \item{component}{
%%     ~~Describe \code{component} here~~
}
  \item{toeplitz_approx_var_cov_wv}{
%%     ~~Describe \code{toeplitz_approx_var_cov_wv} here~~
}
  \item{stochastic_model}{
%%     ~~Describe \code{stochastic_model} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or standard data sets, see data().

## The function is currently defined as
function (x, n_seasonal = 2, vec_earthquakes_relaxation_time = NULL, 
    component = "N", toeplitz_approx_var_cov_wv = TRUE, stochastic_model = "wn + fl") 
{
    if (!stochastic_model \%in\% c("wn + pl", "wn + fl")) {
        stop("Argument stochastic_model should be either 'wn + fl' or 'wn + pl'")
    }
    if (!component \%in\% c("N", "E", "V")) {
        stop("Argument ?´component should be either´ 'N', 'E' or 'V'")
    }
    if (!n_seasonal \%in\% c(1, 2)) {
        stop("Argument `n_seasonal` should take either value `1` or `2`")
    }
    all_mjd_index <- seq(head(x$df_position$modified_julian_day, 
        1), tail(x$df_position$modified_julian_day, 1), by = 1)
    jumps <- c(x$df_equipment_software_changes$modified_julian_date, 
        x$df_earthquakes$modified_julian_date)
    jumps <- unique(jumps)
    vec_earthquakes_index_mjd <- c(x$df_earthquakes$modified_julian_date)
    vec_earthquakes_index_mjd <- unique(vec_earthquakes_index_mjd)
    if (length(jumps) == 0) {
        jumps <- NULL
    }
    X <- create_X_matrix(all_mjd_index = all_mjd_index, jumps = jumps, 
        n_seasonal = n_seasonal, vec_earthquakes_index_mjd = vec_earthquakes_index_mjd, 
        vec_earthquakes_relaxation_time = vec_earthquakes_relaxation_time)
    id_X_sub <- which(rownames(X) \%in\% x$df_position$modified_julian_day)
    X_sub <- X[id_X_sub, ]
    if (component == "N") {
        y <- x$df_position$northings_fractional_portion
    }
    else if (component == "E") {
        y <- x$df_position$eastings_fractional_portion
    }
    else if (component == "V") {
        y <- x$df_position$vertical_fractional_portion
    }
    beta_hat <- .lm.fit(y = y, x = X_sub)$coefficients
    if (n_seasonal == 1) {
        names_beta_hat <- c("Intercept", "Trend", "Sin (Annual)", 
            "Cos (Annual)", if (length(jumps) > 0) paste0("Jump: ", 
                jumps), if (length(vec_earthquakes_index_mjd) > 
                0) paste0("Earthquake: ", vec_earthquakes_index_mjd))
    }
    else if (n_seasonal == 2) {
        names_beta_hat <- c("Intercept", "Trend", "Sin (Annual)", 
            "Cos (Annual)", "Sin (Semi-Annual)", "Cos (Semi-Annual)", 
            if (length(jumps) > 0) paste0("Jump: MJD ", jumps), 
            if (length(vec_earthquakes_index_mjd) > 0) paste0("Earthquake: MJD ", 
                vec_earthquakes_index_mjd))
    }
    names(beta_hat) <- names_beta_hat
    eps_hat_sub <- y - X_sub \%*\% beta_hat
    eps_hat_filled <- vector(mode = "numeric", length = length(all_mjd_index))
    eps_hat_filled[id_X_sub] <- eps_hat_sub
    wv_emp_eps_hat_filled <- wv::wvar(eps_hat_filled)
    vec_omega <- as.numeric(all_mjd_index \%in\% x$df_position$modified_julian_day)
    p_hat <- estimate_p1_p2_mle_cpp(vec_omega)
    pstar_hat <- p_hat[2]/(p_hat[1] + p_hat[2])
    vec_autocov_omega <- create_vec_theo_autocov_omega_cpp(p1 = p_hat[1], 
        p2 = p_hat[2], length(all_mjd_index))
    XtX <- t(X) \%*\% X
    inv_XtX <- Matrix::solve(XtX)
    H <- X \%*\% inv_XtX \%*\% t(X)
    D <- diag(length(all_mjd_index)) - H
    quantities_D <- pre_compute_quantities_on_D_only_required_smarter_cpp(D, 
        approx_type = "3")
    if (all(vec_omega == 1)) {
        no_missing <- TRUE
    }
    else {
        no_missing <- FALSE
    }
    if (stochastic_model == "wn + fl") {
        gamma_init_not_transformed <- find_initial_values_wn_fl(wv_emp_eps_hat_filled)
        gamma_init = log(gamma_init_not_transformed)
    }
    else if (stochastic_model == "wn + pl") {
        gamma_init_not_transformed = find_initial_values_wn_pl(signal = eps_hat_sub, 
            wv_emp = wv_emp_eps_hat_filled)
        gamma_init = c(log(gamma_init_not_transformed[1]), inv_trans_kappa_pl(gamma_init_not_transformed[2]), 
            log(gamma_init_not_transformed[3]))
    }
    res_gmwmx <- optim(par = gamma_init, fn = objective_function_w_missing, 
        wv_obj = wv_emp_eps_hat_filled, n = length(vec_omega), 
        quantities_D = quantities_D, approx_type = "3", vec_autocov_omega = vec_autocov_omega, 
        pstar_hat = pstar_hat, no_missing = no_missing, stochastic_model = stochastic_model)
    if (stochastic_model == "wn + fl") {
        gamma_hat_1 <- exp(res_gmwmx$par)
    }
    else if (stochastic_model == "wn + pl") {
        gamma_hat_1 <- c(exp(res_gmwmx$par[1]), trans_kappa_pl(res_gmwmx$par[2]), 
            exp(res_gmwmx$par[3]))
    }
    if (toeplitz_approx_var_cov_wv) {
        max_J <- wv_emp_eps_hat_filled$J
        min_n_to_get_var_wv <- length(vec_omega) + sum_of_powers_of_2(1, 
            max_J - 1) + 2
        if (stochastic_model == "wn + fl") {
            vec_mean_autocov <- vec_mean_autocov_powerlaw(-1, 
                length(vec_omega)) * gamma_hat_1[2]
            vec_mean_autocov[1] <- vec_mean_autocov[1] + gamma_hat_1[1]
        }
        else if (stochastic_model == "wn + pl") {
            vec_mean_autocov <- powerlaw_autocovariance(kappa = gamma_hat_1[2], 
                sigma2 = gamma_hat_1[3], n = length(vec_omega))
            vec_mean_autocov[1] <- vec_mean_autocov[1] + gamma_hat_1[1]
        }
        vec_mean_autocov_eps_hat <- compute_all_mean_diag_fast_w_linear_interp_only_required_cpp(mat_D_q_term_1 = quantities_D$mat_D_q_term_1, 
            mat_D_q_term_2 = quantities_D$mat_D_q_term_2, sum_on_sub_diag_of_D = quantities_D$sum_on_sub_diag_of_D, 
            vec_autocov = vec_mean_autocov, approx_type = "3")
        if (!no_missing) {
            autocov_wn_fl_times_omega <- vec_mean_autocov_eps_hat * 
                (vec_autocov_omega + pstar_hat^2)
        }
        else {
            autocov_wn_fl_times_omega <- vec_mean_autocov_eps_hat
        }
        autocov_wn_fl_times_omega_w_zeroes <- vector(mode = "numeric", 
            length = min_n_to_get_var_wv)
        autocov_wn_fl_times_omega_w_zeroes[1:(length(vec_omega))] <- autocov_wn_fl_times_omega
        var_cov_omega <- toeplitz(as.vector(vec_autocov_omega))
        Sigma_wv <- get_theo_cov_matrix_wvar_cpp(n = length(vec_omega), 
            autocov_vec_X = autocov_wn_fl_times_omega_w_zeroes)
        inv_var_cov_nu_hat <- Matrix::solve(Sigma_wv)
        if (stochastic_model == "wn + fl") {
            var_cov_mat_wn <- gamma_hat_1[1] * diag(length(vec_omega))
            var_cov_mat_flicker <- var_cov_powerlaw_cpp(sigma2 = gamma_hat_1[2], 
                kappa = -1, n = length(vec_omega))
            var_cov_mat_epsilon <- var_cov_mat_wn + var_cov_mat_flicker
        }
        else if (stochastic_model == "wn + pl") {
            vec_autocov_stationary_powerlaw = powerlaw_autocovariance(kappa = gamma_hat_1[2], 
                sigma2 = gamma_hat_1[3], n = length(vec_omega))
            vec_autocov_stationary_powerlaw[1] <- gamma_hat_1[1] + 
                vec_autocov_stationary_powerlaw[1]
            var_cov_mat_epsilon <- toeplitz(as.vector(vec_autocov_stationary_powerlaw))
        }
    }
    else {
        if (stochastic_model == "wn + fl") {
            var_cov_mat_wn <- gamma_hat_1[1] * diag(length(vec_omega))
            var_cov_mat_flicker <- var_cov_powerlaw_cpp(sigma2 = gamma_hat_1[2], 
                kappa = -1, n = length(vec_omega))
            var_cov_mat_epsilon <- var_cov_mat_wn + var_cov_mat_flicker
        }
        else if (stochastic_model == "wn + pl") {
            vec_autocov_stationary_powerlaw = powerlaw_autocovariance(kappa = gamma_hat_1[2], 
                sigma2 = gamma_hat_1[3], n = length(vec_omega))
            vec_autocov_stationary_powerlaw[1] <- gamma_hat_1[1] + 
                vec_autocov_stationary_powerlaw[1]
            var_cov_mat_epsilon <- toeplitz(as.vector(vec_autocov_stationary_powerlaw))
        }
        var_cov_omega <- toeplitz(as.vector(vec_autocov_omega))
        var_cov_eps_hat_w_missing <- var_cov_mat_epsilon * (var_cov_omega + 
            pstar_hat^2)
        var_cov_nu_hat <- compute_cov_wv_cpp_approx_faster(Sigma_X = var_cov_eps_hat_w_missing)
        inv_var_cov_nu_hat <- Matrix::solve(var_cov_nu_hat)
    }
    res_gmwmx_2 <- optim(par = res_gmwmx$par, fn = objective_function_w_missing, 
        wv_obj = wv_emp_eps_hat_filled, n = length(vec_omega), 
        quantities_D = quantities_D, approx_type = "3", vec_autocov_omega = vec_autocov_omega, 
        omega = inv_var_cov_nu_hat, pstar_hat = pstar_hat, no_missing = no_missing, 
        stochastic_model = stochastic_model)
    if (stochastic_model == "wn + fl") {
        gamma_hat_2 <- exp(res_gmwmx_2$par)
    }
    else if (stochastic_model == "wn + pl") {
        gamma_hat_2 <- c(exp(res_gmwmx_2$par[1]), trans_kappa_pl(res_gmwmx_2$par[2]), 
            exp(res_gmwmx_2$par[3]))
    }
    if (stochastic_model == "wn + fl") {
        vec_mean_autocov <- vec_mean_autocov_powerlaw(kappa = -1, 
            length(vec_omega)) * gamma_hat_2[2]
        vec_mean_autocov[1] <- vec_mean_autocov[1] + gamma_hat_2[1]
    }
    else if (stochastic_model == "wn + pl") {
        vec_mean_autocov <- powerlaw_autocovariance(kappa = gamma_hat_2[2], 
            sigma2 = gamma_hat_2[3], n = length(vec_omega))
        vec_mean_autocov[1] <- vec_mean_autocov[1] + gamma_hat_2[1]
        vec_mean_autocov = as.vector(vec_mean_autocov)
    }
    vec_mean_autocov_eps_hat <- compute_all_mean_diag_fast_w_linear_interp_only_required_cpp(mat_D_q_term_1 = quantities_D$mat_D_q_term_1, 
        mat_D_q_term_2 = quantities_D$mat_D_q_term_2, sum_on_sub_diag_of_D = quantities_D$sum_on_sub_diag_of_D, 
        vec_autocov = vec_mean_autocov, approx_type = "3")
    if (!no_missing) {
        vec_mean_per_diag_w_missing <- vec_mean_autocov_eps_hat * 
            (vec_autocov_omega + pstar_hat^2)
        theo_wv <- autocovariance_to_wv(vec_mean_per_diag_w_missing, 
            tau = wv_emp_eps_hat_filled$scales)
    }
    else {
        theo_wv <- autocovariance_to_wv(vec_mean_autocov_eps_hat, 
            tau = wv_emp_eps_hat_filled$scales)
    }
    if (no_missing) {
        var_cov_beta_hat <- inv_XtX \%*\% t(X) \%*\% var_cov_mat_epsilon \%*\% 
            X \%*\% inv_XtX
    }
    else {
        var_cov_beta_hat <- pstar_hat^(-2) * inv_XtX \%*\% t(X) \%*\% 
            ((var_cov_omega + pstar_hat^2) * var_cov_mat_epsilon) \%*\% 
            X \%*\% inv_XtX
    }
    std_beta_hat_gmwmx_3 <- sqrt(diag(var_cov_beta_hat))
    ret <- list(beta_hat = beta_hat, std_beta_hat = std_beta_hat_gmwmx_3, 
        gamma_hat = gamma_hat_2, vartheta_hat = p_hat, component = component, 
        design_matrix_X = X_sub, y = y, empirical_wvar = wv_emp_eps_hat_filled, 
        theoretical_wvar = theo_wv, df_position = x$df_position, 
        df_earthquakes = x$df_earthquakes, df_equipment_software_changes = x$df_equipment_software_changes, 
        p_hat = p_hat, p_star_hat = pstar_hat, stochastic_model = stochastic_model)
    class(ret) <- "fit_gnss_ts_ngl"
    return(ret)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
