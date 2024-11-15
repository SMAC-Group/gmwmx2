library(gmwmx2)
devtools::load_all()
x = download_station_ngl("RIVO")
plot(x)
plot(x, component = "N")

vec_earthquakes_relaxation_time=NULL
n_seasonal=2
component="N"

# create full index
all_mjd_index <- seq(head(x$df_position$modified_julian_day, 1), tail(x$df_position$modified_julian_day, 1), by = 1)

# create all jumps by combining jumps due to equipment change and jumps due to earthquakes
jumps <- c(
  x$df_equipment_software_changes$modified_julian_date,
  x$df_earthquakes$modified_julian_date
)
jumps

# if multiple jumps  to prevent not invertible matrix
jumps <- unique(jumps)
vec_earthquakes_index_mjd <- c(x$df_earthquakes$modified_julian_date)
vec_earthquakes_index_mjd

# if multiple earthquakes to prevent not invertible matrix
vec_earthquakes_index_mjd <- unique(vec_earthquakes_index_mjd)

if (length(jumps) == 0) {
  jumps <- NULL
}

# create design matrix
X <- create_X_matrix(
  all_mjd_index = all_mjd_index,
  jumps = jumps,
  n_seasonal = n_seasonal,
  vec_earthquakes_index_mjd = vec_earthquakes_index_mjd,
  vec_earthquakes_relaxation_time = vec_earthquakes_relaxation_time
)




dim(X)

# obtain X_sub
id_X_sub <- which(rownames(X) %in% x$df_position$modified_julian_day)
X_sub <- X[id_X_sub, ]

# check individual component of matrix x
plot(x=rownames(X_sub), y=X_sub[,1], type="l")
plot(x=rownames(X_sub), y=X_sub[,2], type="l")
plot(x=rownames(X_sub), y=X_sub[,3], type="l")
plot(x=rownames(X_sub), y=X_sub[,4], type="l")
plot(x=rownames(X_sub), y=X_sub[,5], type="l")
plot(x=rownames(X_sub), y=X_sub[,6], type="l")
plot(x=rownames(X_sub), y=X_sub[,7], type="l")
plot(x=rownames(X_sub), y=X_sub[,8], type="l")
plot(x=rownames(X_sub), y=X_sub[,9], type="l")

# Extract Y given specified component
if (component == "N") {
  y <- x$df_position$northings_fractional_portion
} else if (component == "E") {
  y <- x$df_position$eastings_fractional_portion
} else if (component == "V") {
  y <- x$df_position$vertical_fractional_portion
}

# obtain beta hat
beta_hat <- .lm.fit(y = y, x = X_sub)$coefficients

# create vector of name of parameters
if (n_seasonal == 1) {
  names_beta_hat <- c(
    "Intercept", "Trend", "Sin (Annual)", "Cos (Annual)",
    if (length(jumps) > 0) paste0("Jump: ", jumps),
    if (length(vec_earthquakes_index_mjd) > 0) paste0("Earthquake: ", vec_earthquakes_index_mjd)
  )
} else if (n_seasonal == 2) {
  names_beta_hat <- c(
    "Intercept", "Trend", "Sin (Annual)", "Cos (Annual)",
    "Sin (Semi-Annual)", "Cos (Semi-Annual)",
    if (length(jumps) > 0) paste0("Jump: MJD ", jumps),
    if (length(vec_earthquakes_index_mjd) > 0) paste0("Earthquake: MJD ", vec_earthquakes_index_mjd)
  )
}
# assign names beta hat
names(beta_hat) <- names_beta_hat

# # plot signal and estimated model
plot(x=rownames(X_sub), y=y, type="l")
lines(x=rownames(X_sub), y= X_sub%*% beta_hat, col="red", lwd=2)
beta_hat[2] * 365.25 * 1000

fit = gmwmx2(x, n_seasonal = 2, component = "N", stochastic_model = "wn + fl")
plot(fit)
