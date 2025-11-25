#mamquam_2_spatial_v2

## ===========================
## 1. Libraries & setup
## ===========================
library(terra)
library(sf)
library(tmap)
library(dplyr)
library(zoo)
library(whitebox)
wbt_init()

# Your custom functions
source("R/make_watershed.R")

## ===========================
## 2. DEM prep (Mamquam)
## ===========================
dem_mamq <- rast("spatial data/mamquam_river/DEM.tif")

# Reproject DEM to UTM (same CRS you use everywhere else)
dem_utm_mamq <- project(dem_mamq, "EPSG:3156")  # NAD83(CSRS) / UTM 9N
writeRaster(dem_utm_mamq,
            "spatial data/mamquam_river/dem_utm.tif",
            overwrite = TRUE)

## ===========================
## 3. Calibration points: Owen cross-sections -> pour points
## ===========================
# thalwag_data_m: one point per Owen/Mamquam XS with X/Y, OBJ_ORDER, etc.
load("tabular data/mamquam_river/thalwag_data.RData")

pts_utm <- st_as_sf(
  thalwag_data_m,
  coords = c("POINT_X", "POINT_Y"),
  crs = 3156
)

st_write(pts_utm,
         "spatial data/mamquam_river/pour_points.shp",
         append = FALSE)

## ===========================
## 4. Watershed + area for calibration points, then HG fit
## ===========================
river_name <- "mamquam_river"
path_calib <- sprintf("spatial data/%s/watershed_files/", river_name)

# This should create dem_breached, d8_pntr, fac_area, watersheds, and
# thalwag_with_drainage_area.gpkg in path_calib
watershed_polygon_mamq <- make_watershed(
  path_start   = path_calib,
  dem_utm_path = "spatial data/mamquam_river/dem_utm.tif",
  dem_utm      = dem_utm_mamq,
  pour_points  = "spatial data/mamquam_river/pour_points.shp"
)

# Drop geometry to get watershed areas (ws_area_m2, ws_area_km2)
watershed_df_mamq <- st_drop_geometry(watershed_polygon_mamq)

# Load cross-section hydraulic properties (h, w) for calibration sites
load("tabular data/mamquam_river/mamquam_stream_properties.RData")
# stream_properties_mamq: columns OBJ_ORDER, h, w

# Join XS properties to watershed area by cross-section order
stream_properties_all_mamq <- stream_properties_mamq %>%
  right_join(watershed_df_mamq,
             by = join_by(OBJ_ORDER == watersheds))

# Use km2 consistently in HG
stream_properties_calc_mamq <- stream_properties_all_mamq %>%
  mutate(
    logh = log10(h),
    logA = log10(ws_area_km2)
  )

fit_m  <- lm(logh ~ logA, data = stream_properties_calc_mamq)
alpha_m <- 10^(coef(fit_m)[1])
beta_m  <- coef(fit_m)[2]
summary(fit_m)  # keep for QC

## ===========================
## 5. AOI pour points + watershed (Mamquam mainstem)
## ===========================
path_aoi <- sprintf("spatial data/%s/ws_files_AOI", river_name)

watershed_polygon_m <- make_watershed(
  path_start   = path_aoi,
  dem_utm_path = "spatial data/mamquam_river/dem_utm.tif",
  dem_utm      = dem_utm_mamq,
  pour_points  = "spatial data/mamquam_river/mamq_pour_points_AOI.shp"
)

# Preserve original FID when you read AOI pour points
watershed_area_pnts <- st_read(
  "spatial data/mamquam_river/mamq_pour_points_AOI.shp",
  fid_column_name = "FID"
)

## ===========================
## 6. Thalweg points with drainage area (from make_watershed)
## ===========================
# make_watershed() wrote the snapped AOI pour points with drainage area to:
gpkg_path <- file.path(path_aoi, "thalwag_with_drainage_area.gpkg")
dem_path  <- "spatial data/mamquam_river/dem_utm.tif"

dem <- rast(dem_path)

# See layers for info / debugging
st_layers(gpkg_path)

# Read thalweg points + drainage area; keep FID column from GPKG
thalweg_raw <- st_read(
  gpkg_path,
  layer = "thalwag_with_drainage_area",
  fid_column_name = "FID"
)

# CRS -> DEM CRS, drop Z/M
thalweg <- thalweg_raw %>%
  st_transform(st_crs(crs(dem))) %>%
  st_zm(drop = TRUE, what = "ZM")

## ===========================
## 7. Attach elevation, compute distance & slope
## ===========================
# Elevation from DEM
zvals <- terra::extract(dem, vect(thalweg))
thalweg$elev_m <- zvals[, 2]

# Coordinates
coords <- st_coordinates(thalweg)

thalweg <- thalweg %>%
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
thalweg$S_seg[!is.finite(thalweg$S_seg)] <- NA_real_

# Rolling mean slope (k = 3 points)
k <- 3
thalweg <- thalweg %>%
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

thalweg <- thalweg %>%
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
#   thalweg,
#   file.path(path_aoi, "mamquam_thalweg_D50_v3.gpkg"),
#   delete_dsn = TRUE, 
#   delete_layer = TRUE
# )

## ===========================
## 9. Summaries & diagnostics
## ===========================
# summary_data_raw <- st_read(
#   file.path(path_aoi, "mamquam_thalweg_D50.gpkg"),
#   layer = "mamquam_thalweg_D50",
#   fid_column_name = "FID"
# )

summary_data_df <- st_drop_geometry(thalweg)


# Spawning class counts and percentages
summary_table <- summary_data_df %>%
  group_by(spawn_class) %>%
  summarise(
    n   = n(),
    pct = round(n / nrow(summary_data_df) * 100, 1),
    .groups = "drop"
  )

print(summary_table)

# Segment lengths (approximate)
summary_data_arr <- summary_data_df %>%
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
ggplot(summary_data_df,
       aes(x = S_use, y = D50_mm, color = spawn_class)) +
  geom_point(size = 2, alpha = 0.8) +
  labs(
    title = "Relationship Between Slope and Predicted D50",
    x = "Channel Slope (S_use)",
    y = "Predicted Median Grain Size D50 (mm)",
    color = "Substrate Class"
  ) +
  theme_minimal(base_size = 13)
