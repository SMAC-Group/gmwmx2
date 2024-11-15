

load("scripts/estimate_large_network_hpc/mat_result_simulation_2024-11-15_16-57-06.rda")


# load estimated velocities mida

colnames(df_all_results)=c("station_name", "estimated_trend_N", "std_estimated_trend_N", "estimated_trend_E", "std_estimated_trend_E", "length_signal")

df_all_results = df_all_results %>%
  dplyr::mutate(
    estimated_trend_N = as.numeric(estimated_trend_N),
    std_estimated_trend_N = as.numeric(std_estimated_trend_N),
    estimated_trend_E = as.numeric(estimated_trend_E),
    std_estimated_trend_E = as.numeric(std_estimated_trend_E),
    length_signal = as.numeric(length_signal)
  ) %>%

  dplyr::mutate(

    estimated_trend_N_scaled = estimated_trend_N * 365.25,
    std_estimated_trend_N_scaled = std_estimated_trend_N * 365.25,
    estimated_trend_E_scaled = estimated_trend_E * 365.25,
    std_estimated_trend_E_scaled = std_estimated_trend_E * 365.25
  ) %>% filter(station_name!="KIEL") %>% filter(station_name != "YSDN")


# add latitude and longitude
estimated_velocities_midas = download_estimated_velocities_ngl()

estimated_velocities_midas_sub = estimated_velocities_midas %>% filter(station_name %in% df_all_results$station_name) %>% select(station_name, latitude, longitude)

estimated_velocities_midas_sub$longitude <- ifelse(estimated_velocities_midas_sub$longitude < -180,
                                                   estimated_velocities_midas_sub$longitude + 360,
                                                   estimated_velocities_midas_sub$longitude)

#merge
df_all_results = dplyr::left_join(df_all_results, estimated_velocities_midas_sub)

df_estimated_velocities_gmwmx = df_all_results

df_estimated_velocities_gmwmx = df_estimated_velocities_gmwmx %>% filter(length_signal > 3650)
usethis::use_data(df_estimated_velocities_gmwmx, overwrite = T)


