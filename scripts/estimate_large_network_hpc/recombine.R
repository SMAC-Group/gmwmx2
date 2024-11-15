# define path
folder = "estimate_station_ngl/"
path = "estimate_station_ngl/data_temp"

# list files
all_files = list.files(path = path)

# load first file
load(paste0(path, "/", all_files[1]))
ncol_file = length(vec_to_save)

# create df to save
df_all_results = data.frame(matrix(NA, ncol=ncol_file))
colnames(df_all_results) = names(vec_to_save)

# for all files load and bind
for(file_index in seq_along(all_files)){
  file_i = all_files[file_index]
  file_name = paste0(path,"/",file_i)
  load(file_name)
  df_all_results = rbind(df_all_results, vec_to_save)

}

colnames(df_all_results) = names(vec_to_save)


df_all_results = df_all_results[-1,]


# save matrix of results
time = Sys.time()
time_2 = gsub(" ", "_", time)
time_3 = gsub(":", "-", time_2)
file_name_to_save = paste0(paste("estimate_station_ngl/mat_result_simulation", time_3, sep="_"),
                           ".rda")
print(file_name_to_save)
save(df_all_results, file=file_name_to_save)


