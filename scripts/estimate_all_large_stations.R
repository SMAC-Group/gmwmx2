rm(list=ls())

library(gmwmx2)
library(dplyr)
library(maps)

# Estimate little network in Switzerland
all_station = download_all_stations_ngl()
head(all_station, 10)

# download all 4 stations
# create matrix
df_network = all_station%>% filter(station_name %in% c("CERN", "JUVI", "AUBO", "SATI"))


df_estimated_velocities = data.frame(matrix(NA, nrow=dim(df_network)[1], ncol = 6))
for(station_index in seq_along(df_network$station_name)){

  station_name = df_network$station_name[station_index]
  # extract station
  station_data = download_station_ngl(station_name = station_name)
  fit_N = gmwmx2(station_data, n_seasonal = 2, component = "N", stochastic_model = "wn + pl")
  fit_E = gmwmx2(station_data, n_seasonal = 2, component = "E", stochastic_model = "wn + pl")
  df_estimated_velocities[station_index, 1] = station_name
  df_estimated_velocities[station_index, 2:6] =   c(fit_N$beta_hat[2], fit_N$std_beta_hat[2],fit_E$beta_hat[2], fit_E$std_beta_hat[2], dim(fit_N$design_matrix_X)[1])
  cat(paste0(station_index ,"/", length(df_network$station_name), "\n"))
}

colnames(df_estimated_velocities) = c("station_name", "estimated_trend_N", "std_estimated_trend_N", "estimated_trend_E", "std_estimated_trend_E", "n_data")
df_estimated_velocities$estimated_trend_N_scaled = df_estimated_velocities$estimated_trend_N * 365.25
df_estimated_velocities$std_estimated_trend_N_scaled = df_estimated_velocities$std_estimated_trend_N * 365.25
df_estimated_velocities$estimated_trend_E_scaled = df_estimated_velocities$estimated_trend_E  * 365.25
df_estimated_velocities$std_estimated_trend_E_scaled =  df_estimated_velocities$std_estimated_trend_E * 365.25

# plot
# library(leaflet)
install.packages("elevatr")
library(elevatr)
library(raster)

# Define the extent (longitude/latitude range) you want for the relief map
xlims <- c(44, 48)
ylims <- c(25, 52)
locs <- df_network%>%dplyr::select(latitude, longitude)
locs <- data.frame(x = c(xlims[1], xlims[2]), y = c(ylims[1], ylims[2]))

# Fetch elevation data
relief <- get_elev_raster(locs, prj = "EPSG:4326", z = 6)


df_estimated_velocities_2 = dfplyr::ful
library(leaflet)
# install.packages("leaflet")
leaflet(data = df_estimated_velocities) %>%
  addTiles() %>%
  setView(lng = 6.1432, lat = 46.2044, zoom = 12) %>% # Center on Geneva

  # Add markers for each location
  addMarkers(~lon, ~latit, popup = ~name) %>%

  # Add a circle around Geneva to show the surrounding region
  addCircles(lng = 6.1432, lat = 46.2044, radius = 10000, # 10 km radius
             color = "blue", fillOpacity = 0.2, weight = 1,
             popup = "Geneva Surrounding Region")
library(rnaturalearth)

# Load Switzerland map
# Load Europe map data
europe <- ne_countries(scale = "medium", returnclass = "sf") %>%
  dplyr::filter(admin %in% c("France", "Switzerland", "Germany", "Italy"))

# Specify Geneva's coordinates
geneva <- data.frame(long = 6.1423, lat = 46.2044)

# Arrow length and direction (adjust end longitude and latitude)
arrow_start <- c(6.1423, 46.2044)  # Geneva coordinates
arrow_end <- c(6.2423, 46.2044)
library("rnaturalearth")
library("rnaturalearthdata")
world <- ne_countries(scale = "medium", returnclass = "sf")
ggplot(data = world) +
  geom_sf() +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("World map", subtitle = paste0("(", length(unique(world$NAME)), " countries)"))

library(ggplot2)
ggplot(data = europe) +
  geom_sf() +
  coord_sf(xlim = c(5, 10), ylim = c(45, 48), expand = FALSE) +
  geom_point(aes(x = arrow_start[1], y = arrow_start[2]), color = "red", size = 2) +
  geom_segment(aes(x = arrow_start[1], y = arrow_start[2], xend = arrow_end[1], yend = arrow_end[2]),
               arrow = arrow(length = unit(0.3, "cm")), color = "blue", size = 1) +
  theme_minimal() +
  ggtitle("Map of Geneva and Surrounding Region with Arrow")


fit1 = gmwmx2(x = x, n_seasonal = 2, component = "N", stochastic_model = "wn + pl")
plot(fit1)
summary(fit1)
fit2 = gmwmx2(x = x, n_seasonal = 2, component = "N", stochastic_model = "wn + fl")
plot(fit2)
summary(fit2)
x=download_station_ngl("0ABI")
fit3 = gmwmx2(x = x, n_seasonal = 2, component = "N", stochastic_model = "wn + pl")
plot(fit3)
summary(fit3)
fit4 = gmwmx2(x = x, n_seasonal = 2, component = "N", stochastic_model = "wn + fl")
plot(fit4)
summary(fit4)
x=download_station_ngl("ZWOL")
plot(x)
fit5 = gmwmx2(x = x, n_seasonal = 2, component = "N", stochastic_model = "wn + pl")
plot(fit5)
summary(fit5)
fit6 = gmwmx2(x = x, n_seasonal = 2, component = "N", stochastic_model = "wn + fl")
plot(fit6)
summary(fit6)




# Plot the world map
# Set map limits (for example, focusing on Europe and North America)
xlim <- c(-190, 190)  # Longitude range
ylim <- c(-70, 70)   # Latitude range

# Plot the world map with zoomed-in limits
map("world", fill = TRUE, col = "lightblue", bg = "lightyellow", lwd = 0.5,
    xlim = xlim, ylim = ylim)
range(all_station$longitude)
range(all_station$latitude)

# Fix longitudes that are less than -180 by adding 360
all_station$longitude2 <- ifelse(all_station$longitude < -180,
                                all_station$longitude + 360,
                                all_station$longitude)
range(all_station$longitude2)
# only long
points(all_station$longitude2, all_station$latitude, col = "red", pch = 19, cex = 1.5)















all_station_processed = data.table::fread("http://geodesy.unr.edu/NGLStationPages/GlobalStationList", fill=Inf)

# Load HTML content
html_file <- "path/to/your/large_html_file.html"
page <- read_html(html_file)

# Extract station names using CSS selector for links within the table
station_codes <- page %>%
  html_nodes("a") %>%               # Select all anchor tags
  html_text() %>%                   # Get text from each anchor tag (station code)
  str_trim() %>%                    # Trim whitespace
  str_subset("^(00NA|01NA|02NA)")   # Keep only codes starting with the specified patterns

# View the extracted station codes
print(station_codes)
#


test = download_all_stations_ngl()
x = download_station_ngl(test$station_name[56])
plot(x)
fit=gmwmx2(x = x, n_seasonal = 2)
unique(x$df_position$station_name)





# identify station of interest
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
n_year_max = 10
df3 = df2 %>% filter(time_series_duration_year < n_year_max, n_earthquake==1,
                     n_equipment_change==1)



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
