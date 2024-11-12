rm(list=ls())

library(gmwmx2)
library(dplyr)
library(maps)



test = download_all_stations_ngl()
x = download_station_ngl(test$station_name[56])
plot(x)
fit=gmwmx2(x = x, n_seasonal = 2)
unique(x$df_position$station_name)



all_velocities = download_estimated_velocities_ngl()
name_tmp_dir = tempdir()
name_tmp_file = paste0(name_tmp_dir,"/", "steps.txt")
download.file("http://geodesy.unr.edu/NGLStationPages/steps.txt", name_tmp_file, quiet = TRUE)

# read all lines
all_lines = readLines(name_tmp_file)

# Count words in each line to identify equipment/software changes vs earthquakes
word_counts <- sapply(strsplit(trimws(all_lines), "\\s+"), length)
id_equipment_software_changes = which(word_counts==4)
id_earthquakes = which(word_counts==7)

# extract both dataset
df_equipment_software_changes <- data.table::fread(text = all_lines[id_equipment_software_changes],
                                                   header = FALSE,
                                                   col.names = c("station_name","date_YYMMDD", "step_type_code", "type_equipment_change"))
# subset
df_earthquakes <- data.table::fread(text = all_lines[id_earthquakes],
                                    header = FALSE,
                                    col.names = c("station_name","date_YYMMDD", "step_type_code", "treshold_distance_km", "distance_station_to_epicenter_km", "event_magnitude", "usgs_event_id"))

library(dplyr)
df_count_earthquake = df_earthquakes %>%
  group_by(station_name) %>%
  summarise("n_earthquake" =n())


df_count_equipment_software_changes = df_equipment_software_changes %>%
  group_by(station_name) %>%
  summarise("n_equipment_change" =n())

df_station_time = all_velocities %>% select(station_name, time_series_duration_year)

# merge station duration and number of earthquakes, equipment changes
df1 = dplyr::full_join(df_count_earthquake, df_count_equipment_software_changes)
df2 = dplyr::full_join(df1,df_station_time )

# select only long stations
n_year_min = 5
df3 = df2 %>% filter(time_series_duration_year == n_year_min, n_earthquake==1,
                     n_equipment_change==3)


df3$station_name
station_name_1 = df3$station_name[2]
station_name_1
x = download_station_ngl(station_name_1) # problem here
plot(x)
plot(x, component = "N")
fit = gmwmx2(x, n_seasonal = 2, component = "N")
plot(fit)
# x=download_station_ngl("DARR")
plot(x, component = "E")
fit = gmwmx2(x, n_seasonal = 2, component = "E")
plot(fit)




df_all_station = download_all_stations_ngl()
df_velocities = download_estimated_velocities_ngl()
# get all longest station
df_velocities_long = df_velocities %>% filter(time_series_duration_year>15)

# Plot the world map
# Set map limits (for example, focusing on Europe and North America)
xlim <- c(-190, 190)  # Longitude range
ylim <- c(-70, 70)   # Latitude range

# Plot the world map with zoomed-in limits
map("world", fill = TRUE, col = "lightblue", bg = "lightyellow", lwd = 0.5,
    xlim = xlim, ylim = ylim)


# Add points for latitude and longitude
# points(df_velocities$longitude, df_velocities$latitude, col = "red", pch = 19, cex = 1.5)

# only long
points(df_velocities_long$longitude, df_velocities_long$latitude, col = "red", pch = 19, cex = 1.5)

# check for stations with

df_velocities_long

# Add a title
# check for some names
df_velocities_long$station_name[1:5]
x = download_station_ngl("0ALF")
plot(x)
fit_1 = gmwmx2::gmwmx2(x = x, n_seasonal = 2, component = "V")
