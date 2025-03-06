convert_to_mjd_2 <- function(vec_date) {
  dt <- as.POSIXct(vec_date, format = "%y%b%d", tz = "UTC")
  # Calculate Julian Date
  jd <- as.numeric(dt) / 86400 + 2440587.5

  # Calculate Modified Julian Date
  mjd <- jd - 2400000.5
  return(mjd)
}
