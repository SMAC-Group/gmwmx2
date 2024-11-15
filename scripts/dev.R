library(gmwmx2)
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


#merge
df1 = dplyr::full_join(df_count_earthquake, df_count_equipment_software_changes)
df2 = dplyr::full_join(df1,df_station_time )
df3 = df2 %>% filter(time_series_duration_year <6, n_earthquake>1, n_equipment_change>1)


df3$station_name
x = download_station_ngl("CHML") # problem here
plot(x, component = "N")
fit = gmwmx2(x, n_seasonal = 2, component = "N")
plot(fit)
# x=download_station_ngl("DARR")
plot(x, component = "E")
fit = gmwmx2(x, n_seasonal = 2, component = "E")
plot(fit)
