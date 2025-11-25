
# SPATIAL COMPONENT ---------------------


library(terra)
library(sf)      # vectors / shapefiles
library(tmap)
library(sf)
library(dplyr)
library(zoo)
library(whitebox)
library(mapview)
wbt_init()  # optional but nice if you haven’t run whitebox yet


load("tabular data/owen_crk/thalwag_data.RData")



dem <- rast("spatial data/owen_crk/DEM.tif")

dem_utm <- project(dem, "EPSG:3156")  # NAD83(CSRS) / UTM zone 9N
writeRaster(dem_utm, "spatial data/owen_crk/dem_utm.tif", overwrite = TRUE)



## ===========================
## 3. Calibration points: Owen cross-sections -> pour points
## ===========================
# thalwag_data_m: one point per Owen/Mamquam XS with X/Y, OBJ_ORDER, etc.
load("tabular data/owen_creek/thalwag_data.RData")

pts_utm <- st_as_sf(
  thalwag_data,
  coords = c("POINT_X", "POINT_Y"),
  crs = 3156
)

st_write(pts_utm, "spatial data/owen_crk/pour_points.shp", append = FALSE)

## ===========================
## 4. Watershed + area for calibration points, then HG fit
## ===========================
source("R/make_watershed_2.R")
load("tabular data/owen_crk/owen_stream_properties.RData")
load("tabular data/owen_crk/owen_watershed_df.RData")

path_start_o <- sprintf("spatial data/%s/watershed_files/", "owen_crk") 


# This should create dem_breached, d8_pntr, fac_area, watersheds, and
# thalwag_with_drainage_area.gpkg in path_calib
watershed_polygon_owen <- make_watershed(
  path_start   = path_start_o,
  dem_utm_path = "spatial data/owen_crk/dem_utm.tif",
  dem_utm      = dem_utm,
  pour_points  = "spatial data/owen_crk/pour_points.shp"
)

# Drop geometry to get watershed areas (ws_area_m2, ws_area_km2)
watershed_df_o <- st_drop_geometry(watershed_polygon_owen)

# Load cross-section hydraulic properties (h, w) for calibration sites
load("tabular data/owen_crk/owen_stream_properties.RData")
# stream_properties_mamq: columns OBJ_ORDER, h, w

# Join XS properties to watershed area by cross-section order
stream_properties_all_o <- stream_properties %>%
  right_join(watershed_df_o,
             by = join_by(OBJ_ORDER == watersheds))

# Use km2 consistently in HG
stream_properties_calc <- stream_properties_all_o %>%
  mutate(
    logh = log10(h),
    logA = log10(ws_cum_up_m2)
  )

write.csv(stream_properties_calc, "tabular data/owen_crk/owen_stream_properties.csv")

fit <- lm(logh ~ logA, data=stream_properties_calc)
alpha <- 10^(coef(fit)[1])
beta <- coef(fit)[2]
summary(fit) 



## ===========================
## 5. AOI pour points + watershed (Mamquam mainstem)
## ===========================
path_calib <- sprintf("spatial data/%s/ws_files_AOI_2",  "owen_crk")

watershed_polygon_m <- make_watershed(path_calib,
                                      dem_utm_path = "spatial data/owen_crk/dem_utm.tif",
                                      dem_utm,
                                      "spatial data/owen_crk/owen_pour_points_AOI.shp")


# Preserve original FID when you read AOI pour points
watershed_area_pnts <- st_read(
  "spatial data/owen_crk/owen_xsec_AOI.shp",
  fid_column_name = "FID"
)


## ===========================
## 6. thalweg_o points with drainage area (from make_watershed)
## ===========================
# make_watershed() wrote the snapped AOI pour points with drainage area to:
gpkg_path <- file.path("spatial data/owen_crk/ws_files_AOI/", "thalwag_with_drainage_area.gpkg")
dem_path  <- "spatial data/owen_crk/dem_utm.tif"

dem <- rast(dem_path)

# See layers for info / debugging
st_layers(gpkg_path)

# Read thalweg_o points + drainage area; keep FID column from GPKG
thalweg_o_raw <- st_read(
  gpkg_path,
  layer = "thalwag_with_drainage_area",
  fid_column_name = "FID"
)

# CRS -> DEM CRS, drop Z/M
thalweg_o <- thalweg_o_raw %>%
  st_transform(st_crs(crs(dem))) %>%
  st_zm(drop = TRUE, what = "ZM")



## ===========================
## 7. Attach elevation, compute distance & slope
## ===========================
# Elevation from DEM
zvals <- terra::extract(dem, vect(thalweg_o))
thalweg_o$elev_m <- zvals[, 2]

# Coordinates
coords <- st_coordinates(thalweg_o)

thalweg_o <- thalweg_o %>%
  mutate(
    X = coords[, "X"],
    Y = coords[, "Y"]
  ) %>%
  # Use FID to preserve original creation order (upstream -> downstream)
  arrange(FID) %>%
  mutate(
    dx     = c(0, sqrt(diff(X)^2 + diff(Y)^2)),  # segment length (m)
    dist_m = cumsum(dx),                         # cumulative distance (m)
    dz     = c(NA, diff(elev_m)),
    # downstream-positive slope: if Z decreases downstream, dz is negative
    S_seg  = -dz / dx
  )

# Clean up slope edge cases
thalweg_o$S_seg[!is.finite(thalweg_o$S_seg)] <- NA_real_

# Rolling mean slope (k = 3 points)
k <- 3
thalweg_o <- thalweg_o %>%
  mutate(
    S_roll = zoo::rollapply(
      S_seg, width = k, mean, na.rm = TRUE,
      fill = NA_real_, align = "center"
    ),
    S_use = ifelse(is.na(S_roll), S_seg, S_roll),
    S_use = abs(S_use)  # belt-and-suspenders: ensure slope > 0
  )



## ===========================
## 8. Compute D50 from HG + slope
## ===========================
rho   <- 1000   # kg/m3
rhos  <- 2650   # kg/m3
tau_c <- 0.047  # critical Shields parameter

thalweg_o <- thalweg_o %>%
  mutate(
    # h(A) with area in km² (consistent with HG fit)
    h_m   = alpha_m * (drainage_area_km2^beta_m),
    D50_m = (rho * h_m * S_use) / ((rhos - rho) * tau_c),
    D50_mm = D50_m * 1000,
    spawn_class = case_when(
      is.na(D50_mm)        ~ NA_character_,
      D50_mm < 7           ~ "too fine",
      D50_mm <= 47         ~ "spawning gravel",
      TRUE                 ~ "too coarse"
    )
  )

# # Save for mapping
# st_write(
#   thalweg_o,
#   file.path(path_aoi, "mamquam_thalweg_o_D50_v3.gpkg"),
#   delete_dsn = TRUE, 
#   delete_layer = TRUE
# )

## ===========================
## 9. Summaries & diagnostics
## ===========================
# summary_data_raw <- st_read(
#   file.path(path_aoi, "mamquam_thalweg_o_D50.gpkg"),
#   layer = "mamquam_thalweg_o_D50",
#   fid_column_name = "FID"
# )

summary_data_df_o_o <- st_drop_geometry(thalweg_o_o)


# Spawning class counts and percentages
summary_table <- summary_data_df_o %>%
  group_by(spawn_class) %>%
  summarise(
    n   = n(),
    pct = round(n / nrow(summary_data_df_o) * 100, 1),
    .groups = "drop"
  )

print(summary_table)

# Segment lengths (approximate)
summary_data_arr <- summary_data_df_o %>%
  arrange(FID) %>%  # same order as slope calcs
  mutate(seg_length_m = c(diff(dist_m), 0))

habitat_length <- summary_data_arr %>%
  group_by(spawn_class) %>%
  summarise(
    total_length_m  = sum(seg_length_m, na.rm = TRUE),
    total_length_km = total_length_m / 1000,
    .groups = "drop"
  )

print(habitat_length)

# D50 vs slope plot
ggplot(summary_data_df_o,
       aes(x = S_use, y = D50_mm, color = spawn_class)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(
    title = "Relationship Between Slope and Predicted D50",
    x = "Channel Slope (S_use)",
    y = "Predicted Median Grain Size D50 (mm)",
    color = "Substrate Class"
  ) +
  theme_minimal(base_size = 13)



# # Make sure thalweg_o is in the same CRS as the DEM
# thalweg_o  <- st_transform(watershed_area_pnts, st_crs(crs(dem_utm)))
# thalweg_o  <- st_zm(thalweg_o, drop = TRUE, what = "ZM")  # ensure 2D points
# 
# zvals <- terra::extract(dem_utm, vect(thalweg_o))
# thalweg_o$elev_m <- zvals[, 2]   # 2nd column is the raster value
# 
# 
# 
# 
# 
# ## ---- 4. Compute local + smoothed slope ----
# thalweg_o <- thalweg_o %>%
#   mutate(
#     dz     = c(NA, diff(elev_m)),
#     S_seg  = dz / dx
#   )
# 
# # Rolling-mean slope over ~5 points (adjust 'k' if spacing is very small/large)
# k <- 3
# thalweg_o <- thalweg_o %>%
#   mutate(
#     S_roll = zoo::rollapply(S_seg, width = k, mean, fill = NA, align = "center"),
#     S_use  = ifelse(is.na(S_roll), S_seg, S_roll)  # use smoothed where available
#   )
# 
# ## ---- 5. Compute D50 from hydraulic geometry + slope ----
# # >>>> FILL THESE IN from your h ~ A^beta fit for owen <<<<
# 
# rho   <- 1000    # water density (kg/m3)
# rhos  <- 2650    # sediment density (kg/m3)
# tau_c <- 0.045   # critical Shields (try 0.035–0.06 in sensitivity later)
# 
# thalweg_o <- thalweg_o %>%
#   mutate(
#     h_m    = alpha * (drainage_area_km2^beta),         # bankfull depth from HG
#     D50_m  = (rho * h_m * S_use) / ((rhos - rho) * tau_c),
#     D50_mm = D50_m * 1000,
#     spawn_class = case_when(
#       is.na(D50_mm)        ~ NA_character_,
#       D50_mm < 7           ~ "too fine",
#       D50_mm <= 47         ~ "spawning gravel",
#       TRUE                 ~ "too coarse"
#     )
#   )
# 
# ## ---- 6. Save out for mapping in ArcGIS/QGIS ----
# st_write(thalweg_o,
#          file.path(path_start, "owen_thalweg_o_D50.gpkg"),
#          delete_dsn = TRUE)
# 
# 
# 
# 
