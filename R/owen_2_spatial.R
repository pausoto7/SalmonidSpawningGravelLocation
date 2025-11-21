

# SPATIAL COMPONENT ---------------------


library(terra)
library(sf)      # vectors / shapefiles
library(tmap)
library(sf)
library(dplyr)
library(zoo)
library(whitebox)
wbt_init()  # optional but nice if you haven’t run whitebox yet


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

# 
# tmap_mode("view")  # static map
# 
# tm_shape(dem_utm) +
#   tm_raster(style = "cont", palette = "-viridis", title = "Elevation") +
#   tm_shape(pts_utm) +
#   tm_dots(col = "OBJ_ORDER", size = 0.1, palette = "Set1", title = "X-section") +
#   tm_layout(legend.outside = TRUE)
# 


# BUILD FLOW DIRECTION AND ACCUMULATION --------------------------------------------


source("R/make_watershed.R")


path_start_o <- sprintf("spatial data/%s/ws_files_AOI",  "owen_crk")



watershed_polygon <- make_watershed(path_start_o, 
                                    "spatial data/owen_crk/dem_utm.tif",
                                    dem_utm, 
                                    "spatial data/owen_crk/owen_xsec_AOI.shp")

watershed_df <- st_drop_geometry(watershed_polygon)




save(watershed_df, file = "tabular data/owen_crk/owen_watershed_df.RData")


# -----------------------------------------------------------------



load("tabular data/owen_crk/owen_stream_properties.RData")
load("tabular data/owen_crk/owen_watershed_df.RData")

stream_properties_all <- stream_properties %>%
  right_join(watershed_df, by = join_by(OBJ_ORDER == watersheds))


# math

stream_properties_calc <- stream_properties_all %>%
  mutate(logh = log10(h), 
         logA = log10(ws_area_m2 ))


fit <- lm(logh ~ logA, data=stream_properties_calc)
alpha <- 10^(coef(fit)[1])
beta <- coef(fit)[2]
summary(fit) 



# get final watershed area for input points --------------------------------------------------------



path_start_m <- sprintf("spatial data/%s/ws_files_AOI/",  "owen_crk")


watershed_polygon_m <- make_watershed(path_start_o,
                                      dem_utm_path = "spatial data/owen_crk/dem_utm.tif",
                                      dem_utm,
                                      "spatial data/owen_crk/owen_pour_points_AOI.shp")



watershed_AOI_df_m <- st_drop_geometry(watershed_polygon_m)

watershed_area_pnts <- st_read("spatial data/owen_crk/owen_pour_points_AOI.shp")


# Paths – edit these:
path_start <- "spatial data/owen_crk/ws_files_AOI/"


# Make sure thalweg is in the same CRS as the DEM
thalweg  <- st_transform(watershed_area_pnts, st_crs(crs(dem_utm)))
thalweg  <- st_zm(thalweg, drop = TRUE, what = "ZM")  # ensure 2D points

zvals <- terra::extract(dem_utm, vect(thalweg))
thalweg$elev_m <- zvals[, 2]   # 2nd column is the raster value




## ---- Paths ----
path_start  <- "spatial data/owen_crk/ws_files_AOI/"
gpkg_path   <- file.path(path_start, "thalwag_with_drainage_area.gpkg")
dem_path    <- "spatial data/owen_crk/dem_utm.tif"   # <- make sure this is right

## ---- 1. Read DEM and thalweg+area points ----
dem <- rast(dem_path)

# See what layer name is inside the GPKG
st_layers(gpkg_path)

thalweg <- st_read(gpkg_path, layer = "thalwag_with_drainage_area")

# Make sure CRS matches DEM and drop Z/M
thalweg <- st_transform(thalweg, st_crs(crs(dem)))
thalweg <- st_zm(thalweg, drop = TRUE, what = "ZM")

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
# >>>> FILL THESE IN from your h ~ A^beta fit for owen <<<<

rho   <- 1000    # water density (kg/m3)
rhos  <- 2650    # sediment density (kg/m3)
tau_c <- 0.045   # critical Shields (try 0.035–0.06 in sensitivity later)

thalweg <- thalweg %>%
  mutate(
    h_m    = alpha * (drainage_area_km2^beta),         # bankfull depth from HG
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
         file.path(path_start, "owen_thalweg_D50.gpkg"),
         delete_dsn = TRUE)




