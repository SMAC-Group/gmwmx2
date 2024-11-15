rm(list=ls())
library(gmwmx2)
df_all_stations = download_all_stations_ngl()
df_all_stations$longitude2 <- ifelse(df_all_stations$longitude < -180,
                                     df_all_stations$longitude + 360,
                                     df_all_stations$longitude)
library("rnaturalearth")
library("rnaturalearthdata")


# define function to make color transparent
make_transparent <- function(colors, alpha = 0.5) {
  # Ensure alpha is between 0 and 1
  if(alpha < 0 || alpha > 1) {
    stop("Alpha value must be between 0 and 1")
  }

  # Convert colors to RGB and add alpha
  transparent_colors <- sapply(colors, function(col) {
    rgb_val <- col2rgb(col) / 255
    rgb(rgb_val[1], rgb_val[2], rgb_val[3], alpha = alpha)
  })

  return(transparent_colors)
}

world <- ne_countries(scale = "medium", returnclass = "sf")
library(ggplot2)
my_col = c("#e96bff")

my_col_trans = make_transparent(my_col,alpha = .25)





library(tikzDevice)
tikz(file  = "scripts/graphs/map_all_stations.tex", width = 7, height = 6, standAlone = T)

ggplot(data = world) +
  geom_sf() +
  xlab("Longitude") + ylab("Latitude") +
  ggtitle("NGL GNSS Stations")+
  geom_point(data = df_all_stations, aes(x = longitude2, y = latitude), color =my_col_trans, size = 1)+
  coord_sf(expand = FALSE) +  # Removes extra space around the plot
  theme_minimal() +  # Minimal theme
  theme(
    plot.margin = unit(c(0, 0, 0, 0), "cm"),  # Set plot margin to zero
    axis.text.x = element_blank(),  # Remove x-axis labels
    axis.text.y = element_blank(),   # Remove y-axis labels
    axis.ticks = element_blank() ,
    axis.title = element_blank()
  )

dev.off()
