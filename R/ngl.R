#' Download GNSS position time series and steps reference from Nevada Geodetic Laboratory with IGS14 reference frame.
#'
#' @importFrom data.table fread
#' @importFrom utils read.table
#' @importFrom dplyr filter select
#' @importFrom  utils download.file
#' @importFrom magrittr %>%
#' @param  station_name A \code{string} specifying the station name.
#' @export
#' @examples
#' station_1LSU <- download_station_ngl("1LSU")
#' attributes(station_1LSU)
#' @return A \code{list} of class \code{gnss_ts_ngl} that contains three \code{data.frame}: The \code{data.frame} \code{df_position} which contains the position time series extracted from the .tenv3 file available from the Nevada Geodetic Laboratory, the
#' \code{data.frame} \code{df_equipment_software_changes} which specify the equipment or software changes for that stations and the \code{data.frame} \code{df_earthquakes} that specify the earthquakes associated with that station.
download_station_ngl <- function(station_name) {


  # station_name="1LSU"

  # download file from  string formatted as "http://geodesy.unr.edu/gps_timeseries/tenv3/IGS14/", station_name, ".tenv3"
  # see help file for .tenv3 file here : http://geodesy.unr.edu/gps_timeseries/README_tenv3.txt

  # check that station is available
  df_all_stations <- download_all_stations_ngl()
  if (!station_name %in% df_all_stations$station_name) {
    stop("Invalid station name")
  }


  # -------------- load .tenv3 file using data.table fread
  df_position = data.table::fread(
    paste0("http://geodesy.unr.edu/gps_timeseries/tenv3/IGS14/", station_name, ".tenv3"),
    header = TRUE
  )

  # rewrite colnames
  colnames(df_position) <- c(
    "station_name",
    "date",
    "decimal_year",
    "modified_julian_day",
    "gps_week",
    "day_of_gps_week",
    "longitude_reference_meridian",
    "eastings_integer_portion",
    "eastings_fractional_portion",
    "northings_integer_portion",
    "northings_fractional_portion",
    "vertical_integer_portion",
    "vertical_fractional_portion",
    "antenna_height",
    "east_sigma",
    "north_sigma",
    "vertical_sigma",
    "east_north_correlation",
    "east_vertical_correlation",
    "north_vertical_correlation",
    "nominal_station_latitude",
    "nominal_station_longitude",
    "nominal_station_height"
  )

  # load steps from file steps and extract all steps associated with the station
  # see README for step file: http://geodesy.unr.edu/NGLStationPages/steps_readme.txt

  # # Using fread for fast reading
  dt <- data.table::fread("http://geodesy.unr.edu/NGLStationPages/steps.txt",
              fill = 7,  # Handle varying number of columns
              header = FALSE,  # Assuming no header
              na.strings = "") # Empty fields become NA

  # filter earthquakes from equipment / software changes
  df_equipment_software_changes = dt %>% dplyr::filter(V3 == 1 | V3 == 3) %>% dplyr::select(c("V1", "V2", "V3", "V4"))
  colnames(df_equipment_software_changes)= c("station_name", "date_YYMMDD", "step_type_code", "type_equipment_change")
  df_earthquakes = dt %>% dplyr::filter(V3 == 2 )
  colnames(df_earthquakes) = c("station_name", "date_YYMMDD", "step_type_code", "treshold_distance_km", "distance_station_to_epicenter_km", "event_magnitude", "usgs_event_id")

  # subset
  df_equipment_software_changes_sub <- df_equipment_software_changes |> dplyr::filter(station_name == !!station_name)
  df_earthquakes_sub <- df_earthquakes |> dplyr::filter(station_name == !!station_name)

  # convert to MJD
  df_equipment_software_changes_sub$modified_julian_date <- convert_to_mjd_2(df_equipment_software_changes_sub$date_YYMMDD)
  # convert to MJD
  df_earthquakes_sub$modified_julian_date <- convert_to_mjd_2(df_earthquakes_sub$date_YYMMDD)

  ret <- list(
    "df_position" = df_position,
    "df_equipment_software_changes" = df_equipment_software_changes,
    "df_earthquakes" = df_earthquakes
  )
  class(ret) <- "gnss_ts_ngl"

  return(ret)
}

convert_to_mjd_2 <- function(vec_date) {
  dt <- as.POSIXct(vec_date, format = "%y%b%d", tz = "UTC")
  # Calculate Julian Date
  jd <- as.numeric(dt) / 86400 + 2440587.5

  # Calculate Modified Julian Date
  mjd <- jd - 2400000.5
  return(mjd)
}



#' Download all stations name and location from Nevada Geodetic Laboratory (NGL)
#' @export
#' @return Return a \code{data.frame} with all stations name, latitude, longitude and heights.
download_all_stations_ngl <- function() {
  # load file from http://geodesy.unr.edu/NGLStationPages/llh.out using data.table fread

  # load all stations
  df_all_stations = data.table::fread(
    "http://geodesy.unr.edu/NGLStationPages/llh.out",
    header = FALSE
  )
  colnames(df_all_stations) <- c("station_name", "latitude", "longitude", "height")
  return(df_all_stations)
}



#' Download estimated velocities provided by the Nevada Geodetic Laboratory (NGL) for all stations.
#' @export
#' @return Return a \code{data.frame} with all stations name, information about the time series for each station, estimated velocities and estimated standard deviation of the estimated velocities.
download_estimated_velocities_ngl <- function() {

  # # load file from http://geodesy.unr.edu/velocities/midas.IGS14.txt
  # README for file available at NGL website: http://geodesy.unr.edu/velocities/midas.readme.txt


  df_estimated_velocities_midas = data.table::fread(
    "http://geodesy.unr.edu/velocities/midas.IGS14.txt",
    header = FALSE
  )

  # subset columns
  df_estimated_velocities_midas <- df_estimated_velocities_midas[, c(1, 2, 5, 9, 10, 11, 12, 13, 14, 25, 26, 27)]

  colnames(df_estimated_velocities_midas) <- c(
    "station_name", # Column 1: 4 character station ID
    "midas_version_label", # Column 2: MIDAS version label
    "time_series_duration_year", # Column 5: Time series duration (years)
    "east_velocity_m_yr", # Column 9: East  velocity (m/yr)
    "north_velocity_m_yr", # Column 10: North  velocity (m/yr)
    "up_velocity_m_yr", # Column 11: Up  velocity (m/yr)
    "east_velocity_unc_m_yr", # Column 12: East mode velocity uncertainty (m/yr)
    "north_velocity_unc_m_yr", # Column 13: North mode velocity uncertainty (m/yr)
    "up_velocity_unc_m_yr", # Column 14: Up mode velocity uncertainty (m/yr)
    "latitude", # Column 25: Latitude (degrees)
    "longitude", # Column 26: Longitude (degrees)
    "height" # Column 27: Height (m) of station
  )
  return(df_estimated_velocities_midas)
}





#' Plot a \code{gnss_ts_ngl} object
#' @param x A \code{gnss_ts_ngl} object.
#' @param component A \code{string} with value either "N", "E" or "V" that specify which component to plot (Northing, Easting or Vertical).
#' @param ... Additional graphical parameters.
#' @importFrom utils head tail
#' @importFrom graphics abline box layout legend par plot.new lines
#' @export
#' @examples
#' station_1LSU <- download_station_ngl("1LSU")
#' plot(station_1LSU)
#' plot(station_1LSU, component = "N")
#' plot(station_1LSU, component = "E")
#' plot(station_1LSU, component = "V")
#' @return No return value. Plot a \code{gnss_ts_ngl} object.
plot.gnss_ts_ngl <- function(x, component = NULL, ...) {
  # compute NA over the time series
  all_mjd <- seq(
    head(x$df_position$modified_julian_day, 1),
    tail(x$df_position$modified_julian_day, 1)
  )
  missing_mjd <- all_mjd[which(!all_mjd %in% x$df_position$modified_julian_day)]

  if (is.null(component)) {
    # Save the current graphical parameters
    old_par <- par(no.readonly = TRUE)

    # set parameters for layout
    mat_layout <- matrix(c(1, 2, 3, 4), ncol = 1, nrow = 4)
    layout(mat_layout, heights = c(.1, 1, 1, 1))
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("center",
      horiz = T,
      legend = c("NA", "Equipment/Software change", "Earthquake"),
      col = c("grey60", "blue", "darkorange"),
      pch = c(15, NA, NA),
      pt.cex = c(2, NA, NA),
      # x.intersp = 0.8,
      text.width = c(.1, .3, .1),
      lty = c(NA, 1, 1), bty = "n"
    )
    par(mar = c(4, 4.1, 2, 2.1))

    # north
    plot(x$df_position$modified_julian_day,
      y = x$df_position$northings_fractional_portion, type = "l",
      xlab = "MJD", ylab = "Northing (m)", las = 1
    )
    # add missing data
    for (i in seq_along(missing_mjd)) {
      abline(v = missing_mjd[i], col = "grey60")
    }

    # add equipment change
    for (i in seq((dim(x$df_equipment_software_changes)[1]))) {
      abline(v = x$df_equipment_software_changes$modified_julian_date, col = "blue")
    }

    # add earthquake
    for (i in seq((dim(x$df_earthquakes)[1]))) {
      abline(v = x$df_earthquakes$modified_julian_date, col = "darkorange")
    }
    box()


    # east
    plot(x$df_position$modified_julian_day,
      y = x$df_position$eastings_fractional_portion, type = "l",
      xlab = "MJD", ylab = "Easting (m)", las = 1
    )

    # missing data
    for (i in seq_along(missing_mjd)) {
      abline(v = missing_mjd[i], col = "grey60")
    }

    # add equipment change
    for (i in seq((dim(x$df_equipment_software_changes)[1]))) {
      abline(v = x$df_equipment_software_changes$modified_julian_date, col = "blue")
    }

    # add earthquake
    for (i in seq((dim(x$df_earthquakes)[1]))) {
      abline(v = x$df_earthquakes$modified_julian_date, col = "darkorange")
    }
    box()

    # UP
    plot(x$df_position$modified_julian_day,
      y = x$df_position$vertical_fractional_portion,
      type = "l", xlab = "MJD", ylab = "Vertical (m)", las = 1
    )

    # missing data
    for (i in seq_along(missing_mjd)) {
      abline(v = missing_mjd[i], col = "grey60")
    }

    # add equipment change
    for (i in seq((dim(x$df_equipment_software_changes)[1]))) {
      abline(v = x$df_equipment_software_changes$modified_julian_date, col = "blue")
    }

    # add earthquake
    for (i in seq((dim(x$df_earthquakes)[1]))) {
      abline(v = x$df_earthquakes$modified_julian_date, col = "darkorange")
    }
    box()

    # Restore the original graphical parameters
    par(old_par)

    # Reset the layout if desired
    layout(1) # Reset to a single plot layout
  } else if (!is.null(component)) {
    if (!component %in% c("N", "E", "V")) {
      stop("Component should be either 'N', 'E', or 'V'")
    }


    # Save the current graphical parameters
    old_par <- par(no.readonly = TRUE)

    # set parameters for layout
    mat_layout <- matrix(c(1, 2), ncol = 1, nrow = 2)
    layout(mat_layout, heights = c(.1, 1, 1, 1))
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend("center",
      horiz = T,
      legend = c("NA", "Equipment/Software change", "Earthquake"),
      col = c("grey60", "blue", "darkorange"),
      pch = c(15, NA, NA),
      pt.cex = c(2, NA, NA),
      text.width = c(.1, .4, .1),
      lty = c(NA, 1, 1), bty = "n"
    )
    par(mar = c(4, 4.1, 2, 2.1))
    if (component == "N") {
      # north
      plot(x$df_position$modified_julian_day,
        y = x$df_position$northings_fractional_portion, type = "l",
        xlab = "MJD", ylab = "Northing (m)", las = 1
      )
      # add missing data
      for (i in seq_along(missing_mjd)) {
        abline(v = missing_mjd[i], col = "grey60")
      }

      # add equipment change
      for (i in seq((dim(x$df_equipment_software_changes)[1]))) {
        abline(v = x$df_equipment_software_changes$modified_julian_date, col = "blue")
      }

      # add earthquake
      for (i in seq((dim(x$df_earthquakes)[1]))) {
        abline(v = x$df_earthquakes$modified_julian_date, col = "darkorange")
      }
      box()
    } else if (component == "E") {
      # east
      plot(x$df_position$modified_julian_day,
        y = x$df_position$eastings_fractional_portion, type = "l",
        xlab = "MJD", ylab = "Easting (m)", las = 1
      )

      # missing data
      for (i in seq_along(missing_mjd)) {
        abline(v = missing_mjd[i], col = "grey60")
      }

      # add equipment change
      for (i in seq((dim(x$df_equipment_software_changes)[1]))) {
        abline(v = x$df_equipment_software_changes$modified_julian_date, col = "blue")
      }

      # add earthquake
      for (i in seq((dim(x$df_earthquakes)[1]))) {
        abline(v = x$df_earthquakes$modified_julian_date, col = "darkorange")
      }
      box()
    } else if (component == "V") {
      # UP
      plot(x$df_position$modified_julian_day,
        y = x$df_position$vertical_fractional_portion,
        type = "l", xlab = "MJD", ylab = "Vertical (m)", las = 1
      )

      # missing data
      for (i in seq_along(missing_mjd)) {
        abline(v = missing_mjd[i], col = "grey60")
      }

      # add equipment change
      for (i in seq((dim(x$df_equipment_software_changes)[1]))) {
        abline(v = x$df_equipment_software_changes$modified_julian_date, col = "blue")
      }

      # add earthquake
      for (i in seq((dim(x$df_earthquakes)[1]))) {
        abline(v = x$df_earthquakes$modified_julian_date, col = "darkorange")
      }
      box()
    }

    # Restore the original graphical parameters
    par(old_par)

    # Reset the layout if desired
    layout(1) # Reset to a single plot layout
  }
}
