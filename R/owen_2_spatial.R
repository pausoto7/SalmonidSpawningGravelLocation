

# SPATIAL COMPONENT ---------------------

library(terra)
library(sf)      # vectors / shapefiles
library(tmap)

load("tabular data/owen_crk/thalwag_data.RData")

dem <- rast("spatial data/owen_crk/DEM.tif")


dem_utm <- project(dem, "EPSG:3156")  # NAD83(CSRS) / UTM zone 9N
writeRaster(dem_utm, "spatial data/owen_crk/dem_utm.tif", overwrite = TRUE)

# your thalwag points in UTM (you basically have this already)
pts_utm <- st_as_sf(
  thalwag_data,
  coords = c("POINT_X", "POINT_Y"),
  crs = 3156
)
st_write(pts_utm, "spatial data/owen_crk/pour_points.shp", append = FALSE)


tmap_mode("view")  # static map

tm_shape(dem_utm) +
  tm_raster(style = "cont", palette = "-viridis", title = "Elevation") +
  tm_shape(pts_utm) +
  tm_dots(col = "OBJ_ORDER", size = 0.1, palette = "Set1", title = "X-section") +
  tm_layout(legend.outside = TRUE)



# BUILD FLOW DIRECTION AND ACCUMULATION --------------------------------------------

library(whitebox)
wbt_init()  # optional but nice if you haven’t run whitebox yet


source("R/make_watershed.R")

river_name <- "owen_crk"

make_watershed(river_name, "spatial data/owen_crk/dem_utm.tif", dem_utm, "spatial data/owen_crk/pour_points.shp")



