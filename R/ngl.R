#'Load .tenv3 GNSS position time series format and steps reference from Nevada Geodetic Laboratory with IGS14  reference frame.
#'
#' @importFrom data.table fread
#' @importFrom utils read.table
#' @importFrom dplyr filter
#' @importFrom magrittr %>%
#' @param  station_name A \code{string} specifying the station name.
#' @export
#' @examples
#' station_1LSU = download_station_ngl("1LSU")
#' str(station_1LSU)
#' @return A \code{list} of class \code{gnss_ts_ngl} that contains three \code{data.frame}: The \code{data.frame} \code{df_position} which contains the position time series extracted from the .tenv3 file available from the Nevada Geodetic Laboratory, the
#' \code{data.frame} \code{df_equipment_software_changes} which specify the equipment or software changes for that stations and the \code{data.frame} \code{df_earthquakes} that specfiy the earthquakes associated with that station.
download_station_ngl = function(station_name) {
  # see help file for .tenv3 file here : http://geodesy.unr.edu/gps_timeseries/README_tenv3.txt
  # station_name="1LSU"
  # check that station is available
  all_stations = download_all_stations_names_ngl()
  if(! station_name %in% all_stations){
    stop("Invalid station name")
  }
  # -------------- load .tenv3 file
  name_tmp_dir = tempdir()
  name_tmp_file = paste0(name_tmp_dir,"/", station_name, ".tenv3")
  cmd = paste0("wget -O ", name_tmp_file," ", "http://geodesy.unr.edu/gps_timeseries/tenv3/IGS14/", station_name, ".tenv3")
  # run command
  system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  # extract file from temporary
  df_position = read.table(name_tmp_file, header = T)
  # rewrite colnames
  colnames(df_position) =  c(
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
  name_tmp_dir = tempdir()
  name_tmp_file = paste0(name_tmp_dir,"/", "steps.txt")
  cmd = paste0("wget -O ", name_tmp_file," ", " http://geodesy.unr.edu/NGLStationPages/steps.txt")

  # run command
  system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

  # read all lines
  all_lines = readLines(name_tmp_file)

  # Count words in each line
  word_counts <- sapply(strsplit(trimws(all_lines), "\\s+"), length)
  id_equipment_software_changes = which(word_counts==4)
  id_earthquakes = which(word_counts==7)

  # extract both dataset
  df_equipment_software_changes <- data.table::fread(text = all_lines[id_equipment_software_changes],
                                header = FALSE,
                                col.names = c("station_name","date_YYMMDD", "step_type_code", "type_equipment_change"))
  # subset
  df_equipment_software_changes = df_equipment_software_changes  |> dplyr::filter(station_name == !!station_name)
  # subset
  df_earthquakes <- data.table::fread(text = all_lines[id_earthquakes],
                                                     header = FALSE,
                                                     col.names = c("station_name","date_YYMMDD", "step_type_code", "treshold_distance_km", "distance_station_to_epicenter_km", "event_magnitude", "usgs_event_id"))
  df_earthquakes = df_earthquakes |> dplyr::filter(station_name == !!station_name)
  ret = list("df_position" = df_position,
             "df_equipment_software_changes" = df_equipment_software_changes,
             "df_earthquakes" = df_earthquakes)
  class(ret) = "gnss_ts_ngl"

  return(ret)
}




#' Download all stations name from Nevada Geodetic Laboratory
#' @export
download_all_stations_names_ngl = function(){
  # load file from http://geodesy.unr.edu/NGLStationPages/llh.out
  name_tmp_dir = tempdir()
  name_tmp_file = paste0(name_tmp_dir,"/", "all_stations.txt")
  cmd = paste0("wget -O ", name_tmp_file," ", " http://geodesy.unr.edu/NGLStationPages/llh.out")

  # run command
  system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)

  # read all lines
  df_all_stations = read.table(name_tmp_file)
  vec_all_stations = df_all_stations[,1]
  return(vec_all_stations)
}




