

# SPATIAL COMPONENT ---------------------

library(terra)
library(sf)      # vectors / shapefiles
library(tmap)
library(sf)
library(dplyr)
library(zoo)
library(whitebox)
wbt_init()  # optional but nice if you haven’t run whitebox yet

load("tabular data/mamquam_river/thalwag_data.RData")

dem_mamq <- rast("spatial data/mamquam_river/DEM.tif")


# dem_utm_mamq <- project(dem_mamq, "EPSG:3156")  # NAD83(CSRS) / UTM zone 9N
# writeRaster(dem_utm_mamq, "spatial data/mamquam_river/dem_utm.tif", overwrite = TRUE)
# 
# # your thalwag points in UTM (you basically have this already)
# pts_utm <- st_as_sf(
#   thalwag_data_m,
#   coords = c("POINT_X", "POINT_Y"),
#   crs = 3156
# )
# st_write(pts_utm, "spatial data/mamquam_river/pour_points.shp", append = FALSE)
# 

# tmap_mode("view")  # static map
# 
# tm_shape(dem_utm_mamq) +
#   tm_raster(style = "cont", palette = "-viridis", title = "Elevation") +
#   tm_shape(pts_utm) +
#   tm_dots(col = "OBJ_ORDER", size = 0.1, palette = "Set1", title = "X-section") +
#   tm_layout(legend.outside = TRUE)
# 


# BUILD FLOW DIRECTION AND ACCUMULATION --------------------------------------------



#source("R/make_watershed.R")

river_name <- "mamquam_river"

path_start <- sprintf("spatial data/%s/watershed_files_v2", river_name)


watershed_polygon_mamq <- make_watershed(path_start, "spatial data/mamquam_river/dem_utm.tif", dem_utm_mamq, "spatial data/mamquam_river/pour_points.shp")

watershed_df_mamq <- st_drop_geometry(watershed_polygon_mamq)


#### MATH _--------------------------------------
load(file = "tabular data/mamquam_river/owen_stream_properties.RData")


stream_properties_all_mamq <- stream_properties_mamq %>%
  right_join(watershed_df_mamq, by = join_by(OBJ_ORDER == watersheds))


# math

stream_properties_calc_mamq <- stream_properties_all_mamq %>%
  mutate(logh = log10(h), 
         logA = log10(ws_up_A_m2 ))


stream_properties_calc_mamq <- stream_properties_calc_mamq %>%
  filter(OBJ_ORDER <4)

fit_m <- lm(logh ~ logA, data=stream_properties_calc_mamq)
alpha_m <- 10^(coef(fit_m)[1])
beta_m <- coef(fit_m)[2]
summary(fit_m) 







# get final watershed area for input points --------------------------------------------------------



path_start_m <- sprintf("spatial data/%s/ws_files_AOI",  "mamquam_river")


watershed_polygon_m <- make_watershed(path_start_m,
                                      dem_utm_path = "spatial data/mamquam_river/dem_utm.tif",
                                      dem_utm_mamq,
                                      "spatial data/mamquam_river/mamq_pour_points_AOI.shp")

watershed_polygon_m <- sf::st_read("spatial data/mamquam_river/mamq_pour_points_AOI.shp", fid_column_name = "FID")


watershed_AOI_df_m <- st_drop_geometry(watershed_polygon_m)

watershed_area_pnts <- st_read("spatial data/mamquam_river/mamq_pour_points_AOI.shp")


# Paths – edit these:
path_start <- "spatial data/mamquam_river/ws_files_AOI/"


# Make sure thalweg is in the same CRS as the DEM
thalweg  <- st_transform(watershed_area_pnts, st_crs(crs(dem_utm_mamq)))
thalweg  <- st_zm(thalweg, drop = TRUE, what = "ZM")  # ensure 2D points

zvals <- terra::extract(dem_utm_mamq, vect(thalweg))
thalweg$elev_m <- zvals[, 2]   # 2nd column is the raster value




## ---- Paths ----
path_start  <- "spatial data/mamquam_river/ws_files_AOI/"
gpkg_path   <- file.path(path_start, "thalwag_with_drainage_area.gpkg")
dem_path    <- "spatial data/mamquam_river/dem_utm.tif"   # <- make sure this is right

## ---- 1. Read DEM and thalweg+area points ----
dem <- rast(dem_path)

# See what layer name is inside the GPKG
st_layers(gpkg_path)

thalweg_raw <- st_read(gpkg_path, layer = "thalwag_with_drainage_area")

# Make sure CRS matches DEM and drop Z/M
thalweg_transformed <- st_transform(thalweg_raw, st_crs(crs(dem)))
thalweg <- st_zm(thalweg_transformed, drop = TRUE, what = "ZM")

## ---- 2. Extract elevation from DEM ----
zvals <- terra::extract(dem, vect(thalweg))
thalweg$elev_m <- zvals[, 2]  # 2nd column is DEM value

## ---- 3. Add coordinates + along-channel distance ----
coords <- st_coordinates(thalweg)

thalweg <- thalweg %>%
  mutate(
    X = coords[, "X"],
    Y = coords[, "Y"]
  ) %>%
  # Order from upstream to downstream: smaller drainage area first
  arrange(drainage_area_km2) %>%
  mutate(
    dx     = c(0, sqrt(diff(X)^2 + diff(Y)^2)),  # segment length
    dist_m = cumsum(dx)                          # cumulative distance
  )

## ---- 4. Compute local + smoothed slope ----
thalweg <- thalweg %>%
  mutate(
    dz     = c(NA, diff(elev_m)),
    S_seg  = dz / dx
  )

# Rolling-mean slope over ~5 points (adjust 'k' if spacing is very small/large)
k <- 3
thalweg <- thalweg %>%
  mutate(
    S_roll = zoo::rollapply(S_seg, width = k, mean, fill = NA, align = "center"),
    S_use  = ifelse(is.na(S_roll), S_seg, S_roll)  # use smoothed where available
  )

## ---- 5. Compute D50 from hydraulic geometry + slope ----
# >>>> FILL THESE IN from your h ~ A^beta fit for Mamquam <<<<

rho   <- 1000    # water density (kg/m3)
rhos  <- 2650    # sediment density (kg/m3)
tau_c <- 0.045   # critical Shields (try 0.035–0.06 in sensitivity later)

thalweg <- thalweg %>%
  mutate(
    h_m    = alpha_m * (drainage_area_km2^beta_m),         # bankfull depth from HG
    D50_m  = (rho * h_m * S_use) / ((rhos - rho) * tau_c),
    D50_mm = D50_m * 1000,
    spawn_class = case_when(
      is.na(D50_mm)        ~ NA_character_,
      D50_mm < 7           ~ "too fine",
      D50_mm <= 47         ~ "spawning gravel",
      TRUE                 ~ "too coarse"
    )
  )

## ---- 6. Save out for mapping in ArcGIS/QGIS ----
st_write(thalweg,
         file.path(path_start, "mamquam_thalweg_D50.gpkg"),
         delete_dsn = TRUE)

# ----------------------------------------------

summary_data_raw <- sf::st_read("spatial data/mamquam_river/ws_files_AOI/mamquam_thalweg_D50.gpkg", layer = "mamquam_thalweg_D50")

summary_data_df <- sf::st_drop_geometry(summary_data_raw)

# Summary of spawning classes
summary_table <- summary_data_df %>%
  filter()
  group_by(spawn_class) %>%
  summarise(
    n = n(),
    pct = round(n / nrow(summary_data_df) * 100, 1)
  )

summary_table



# Estimate segment lengths as distance differences
summary_data_arr <- summary_data_df %>%
  arrange(dist_m) %>%
  mutate(seg_length_m = c(diff(dist_m), 0))

habitat_length <- summary_data_arr %>%
  group_by(spawn_class) %>%
  summarise(
    total_length_m = sum(seg_length_m, na.rm = TRUE),
    total_length_km = total_length_m / 1000
  )

habitat_length

ggplot(summary_data_df, aes(x = S_use, y = D50_mm, color = spawn_class)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(
    title = "Relationship Between Slope and Predicted D50",
    x = "Channel Slope (S_use)",
    y = "Predicted Median Grain Size D50 (mm)",
    color = "Substrate Class"
  ) +
  theme_minimal(base_size = 13)




