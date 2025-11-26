
source("R/make_watershed_3.R")


library(tidyverse)
library(terra)
library(sf)      # vectors / shapefiles
library(tmap)
library(sf)
library(zoo)
library(whitebox)
library(mapview)

dem <- rast("spatial data/owen_crk/DEM.tif")

dem_utm <- project(dem, "EPSG:3156")  # NAD83(CSRS) / UTM zone 9N
# This should create dem_breached, d8_pntr, fac_area, watersheds, and
# thalwag_with_drainage_area.gpkg in path_calib
watershed_polygon_owen_new <- make_watershed(
  path_start   = "spatial data/owen_crk/new_workflow",
  dem_utm_path = "spatial data/owen_crk/dem_utm.tif",
  dem_utm      = dem_utm,
  pour_points  = "spatial data/owen_crk/Owen Crk Pour Points.shp"
)

drainage_df_owen <- st_drop_geometry(watershed_polygon_owen_new)

new_pp_owen <- st_read("spatial data/owen_crk/Owen Crk Pour Points.shp")

# Make sure CRS matches DEM and drop Z/M
pp_transformed <- st_transform(new_pp_owen, st_crs(crs(dem)))
pp_owen_transformed <- st_zm(pp_transformed, drop = TRUE, what = "ZM")


zvals <- terra::extract(dem, vect(pp_owen_transformed))
pp_owen_transformed$elev_m <- zvals[, 2]  # 2nd column is DEM value


all_new_owen <- pp_owen_transformed %>%
  mutate(Name = as.integer(Name)) %>%
  full_join(drainage_df_owen, by = join_by(Name == watersheds) ) 


# extract coordinates in Name order
coords <- all_new_owen %>%
  arrange(Name) %>%
  st_coordinates()

# compute geodesic distances (meters)
dists_m <- geosphere::distHaversine(
  coords[-nrow(coords), ], 
  coords[-1, ]
)

# assemble output table
distance_table <- all_new_owen %>%
  arrange(as.numeric(Name)) %>%
  slice(1:(n() - 1)) %>%
  transmute(
    from = Name,
    to   = lead(Name),
    distance_m = round(dists_m, 1),
    distance_km = round(dists_m / 1000, 4)
  ) %>%
  st_drop_geometry() %>%
  full_join(all_new_owen, by = join_by(from == Name))


slope_table <- distance_table %>%
  arrange(as.numeric(from)) %>%
  mutate(
    elev_this = elev_m,
    elev_next = lead(elev_m),
    slope_m_m = (elev_next - elev_this) / distance_m  ) %>%
    rename(stream_width = PopupInfo) %>%
  select(from, to, distance_m, elev_this, elev_next, slope_m_m, stream_width, ws_cum_up_km2)


# CALCULATE DRAINAGE AREA  ----------------------------------------------------------
# Derive regional Q–A relationship from similar watersheds
# CSV must have columns: station, area (km2), flow (m3/s)

qa_regional <- read.csv("spatial data/owen_crk/owen_similar_watersheds.csv",
                        stringsAsFactors = FALSE)

qa_regional_clean <- qa_regional %>%
  dplyr::rename(
    area_km2 = area,
    Q_m3s    = flow
  ) %>%
  dplyr::filter(
    !is.na(area_km2), !is.na(Q_m3s),
    area_km2 > 0, Q_m3s > 0
  )

# Fit Q = a * A^b  →  ln(Q) = ln(a) + b * ln(A)
fit_QA <- lm(log(Q_m3s) ~ log(area_km2), data = qa_regional_clean)
coef_QA <- coef(fit_QA)

a <- exp(coef_QA[1])  # alpha (units depend on area units; here A in km², Q in m³/s)
b <- coef_QA[2]       # beta (dimensionless)

# Optional: quick check in console
message("Regional Q–A fit for Owen: Q = ", round(a, 3), " * A^", round(b, 3))

# Hydraulic constants
k_s <- 0.3  # roughness coefficient
g   <- 9.81 # gravity, m/s^2

# Use regional Q–A curve to estimate Q and unit discharge q at each cross-section
owen_w_discharge <- slope_table %>%
  mutate(
    stream_width = as.numeric(stream_width),
    Q            = a * ws_cum_up_km2^b,  # Q in m3/s, A in km2
    q            = Q / stream_width      # unit discharge (m2/s)
  )

owen_hydro <- owen_w_discharge %>%
  mutate(
    height = (
      (k_s^(1/3) * q^2 * stream_width) /
        (64 * g * slope_m_m)
    )^(3/10)
  )

hg_fit <- owen_hydro %>%
  filter(
    !is.na(ws_cum_up_km2),
    !is.na(height),
    ws_cum_up_km2 > 0,
    height > 0
  ) %>%
  lm(log(height) ~ log(ws_cum_up_km2), data = .)

coef_hg <- coef(hg_fit)
alpha   <- exp(coef_hg[1])   # alpha in h = alpha * A^beta
beta    <- coef_hg[2]        # beta (dimensionless)

message("Hydraulic geometry: h = ",
        round(alpha, 3), " * A^", round(beta, 3),
        "   (A in km2, h in m)")



write.csv(owen_hydro, "tabular data/owen_crk/owen_crk_hydro_output.csv")




# NEW DATA PORTION ----------------------------------------------------------------------

# Preserve original FID when you read AOI pour points
watershed_area_pnts <- st_read(
  "spatial data/owen_crk/owen_xsec_AOI.shp",
  fid_column_name = "FID"
) 

watershed_area_pnts <- watershed_area_pnts %>%
  filter(as.numeric(FID) %% 2 ==1)

st_write(
  watershed_area_pnts,
  "spatial data/owen_crk/watershed_area_pnts_clean.shp"
)

watershed_area_pnts <- st_read(
  "spatial data/owen_crk/watershed_area_pnts_clean5.shp", 
  fid_column_name = "FID"
)

watershed_AOI <- make_watershed(
  path_start   = "spatial data/owen_crk/new_workflow",
  dem_utm_path = "spatial data/owen_crk/dem_utm.tif",
  dem_utm      = dem_utm,
  pour_points  = "spatial data/owen_crk/watershed_area_pnts_clean5.shp"
)

ws_AOI_df <- watershed_AOI %>%
  st_drop_geometry() %>%
  dplyr::rename(
    FID = watersheds,                # adjust if make_watershed uses a different id
    ws_cum_up_km2 = ws_cum_up_km2    # or whatever column name you chose
  )

# join drainage area onto your AOI thalweg points
thalweg_o <- thalweg_o %>%
  mutate(FID = as.integer(FID)) %>%
  left_join(ws_AOI_df, by = "FID")




# Make sure thalweg_o is in the same CRS as the DEM
thalweg_o  <- st_transform(watershed_area_pnts, st_crs(crs(dem_utm)))
thalweg_o  <- st_zm(thalweg_o, drop = TRUE, what = "ZM")  # ensure 2D points

zvals <- terra::extract(dem_utm, vect(thalweg_o))
thalweg_o$elev_m <- zvals[, 2]   # 2nd column is the raster value


## ---- 4. Compute local + smoothed slope ----
thalweg_o <- thalweg_o %>%
  mutate(
    dz     = c(NA, diff(elev_m)),
    S_seg  = dz / dx
  )

# Rolling-mean slope over ~5 points (adjust 'k' if spacing is very small/large)
k <- 3
thalweg_o <- thalweg_o %>%
  mutate(
    S_roll = zoo::rollapply(S_seg, width = k, mean, fill = NA, align = "center"),
    S_use  = ifelse(is.na(S_roll), S_seg, S_roll)  # use smoothed where available
  )

## ---- 5. Compute D50 from hydraulic geometry + slope ----
# >>>> FILL THESE IN from your h ~ A^beta fit for owen <<<<

rho   <- 1000    # water density (kg/m3)
rhos  <- 2650    # sediment density (kg/m3)
tau_c <- 0.045   # critical Shields (try 0.035–0.06 in sensitivity later)

thalweg_o <- thalweg_o %>%
  mutate(
    h_m    = alpha * (ws_cum_up_km2^beta),         # bankfull depth from HG
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
st_write(thalweg_o,
         file.path(path_start, "owen_thalweg_o_D50.gpkg"),
         delete_dsn = TRUE)





