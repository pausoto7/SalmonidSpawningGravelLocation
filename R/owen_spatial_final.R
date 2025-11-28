# ============================================================
# PREDICTING D50 FOR OWEN CREEK USING DEM + HYDRAULIC GEOMETRY
# ============================================================

# This script does two main things:
#  1) Uses calibration points with known widths + DEM + regional Q–A
#     to estimate depth via a resistance equation, then fit
#     hydraulic geometry h = alpha * A^beta.
#  2) Applies that hydraulic geometry to a denser set of AOI
#     cross-section points to compute D50 and classify spawning habitat.
#
# Make sure "R/make_watershed_3.R" is in your working directory
# and that all file paths below exist.

# ------------------------------------------------------------
# 0. SETUP
# ------------------------------------------------------------

source("R/make_watershed_3.R")

library(tidyverse)
library(terra)       # rasters
library(sf)          # vectors / shapefiles
library(tmap)
library(zoo)         # rollapply for smoothing slopes
library(whitebox)    # watershed / hydrologic tools
library(mapview)
library(geosphere)   # distHaversine for distance in lat/long

# Base path for outputs from make_watershed and final GPKG
path_start <- "spatial data/owen_crk/new_workflow"


# ------------------------------------------------------------
# 1. LOAD DEM AND BUILD WATERSHEDS FOR CALIBRATION POINTS
# ------------------------------------------------------------

# 1.1 Load original DEM (likely in geographic coordinates)
dem <- rast("spatial data/owen_crk/DEM.tif")

# 1.2 Project DEM to UTM for WhiteboxTools (hydrologic processing)
dem_utm <- project(dem, "EPSG:3156")  # NAD83(CSRS) / UTM zone 9N

# 1.3 Save UTM DEM to disk if not already present (make_watershed uses path)
writeRaster(
  dem_utm,
  "spatial data/owen_crk/dem_utm.tif",
  overwrite = TRUE
)

# 1.4 Run custom watershed function for calibration pour points
#     This should create hydrologic derivatives and return polygons
#     with a field like "watersheds" and cumulative area ws_cum_up_km2.
watershed_polygon_owen_new <- make_watershed(
  path_start   = path_start,
  dem_utm_path = "spatial data/owen_crk/dem_utm.tif",
  dem_utm      = dem_utm,
  pour_points  = "spatial data/owen_crk/Owen Crk Pour Points.shp"
)

# 1.5 Drop geometry to retain just attributes (e.g., watersheds, ws_cum_up_km2)
drainage_df_owen <- st_drop_geometry(watershed_polygon_owen_new)


# ------------------------------------------------------------
# 2. CALIBRATION POINTS: ADD ELEVATION AND DISTANCES
# ------------------------------------------------------------

# 2.1 Load calibration pour points (these have widths in PopupInfo)
new_pp_owen <- st_read("spatial data/owen_crk/Owen Crk Pour Points.shp")

# 2.2 Reproject pour points to match the original DEM CRS (for elevation extraction)
pp_transformed <- st_transform(new_pp_owen, st_crs(crs(dem)))

# 2.3 Drop any Z/M dimension to keep simple 2D points
pp_owen_transformed <- st_zm(pp_transformed, drop = TRUE, what = "ZM")

# 2.4 Extract elevation from original DEM at each pour point
#     terra::extract returns a data.frame with point ID + raster values.
zvals <- terra::extract(dem, vect(pp_owen_transformed))
pp_owen_transformed$elev_m <- zvals[, 2]  # second column is DEM value

# 2.5 Join pour points to watershed attributes (cumulative area, etc.)
#     Name is the point ID; watersheds is polygon ID from make_watershed.
all_new_owen <- pp_owen_transformed %>%
  mutate(Name = as.integer(Name)) %>%
  full_join(drainage_df_owen, by = join_by(Name == watersheds))


# ------------------------------------------------------------
# 3. DISTANCES AND SLOPES AT CALIBRATION POINTS
# ------------------------------------------------------------

# 3.1 Extract coordinates in Name order (assumes Name follows channel order)
coords <- all_new_owen %>%
  arrange(Name) %>%
  st_coordinates()

# 3.2 Compute geodesic distances between successive points (in meters)
#     distHaversine expects lon/lat, so this assumes DEM / points are in geographic CRS.
dists_m <- geosphere::distHaversine(
  coords[-nrow(coords), ],
  coords[-1, ]
)

# 3.3 Build a table with from→to segments, distances, and join back to attributes
distance_table <- all_new_owen %>%
  arrange(as.numeric(Name)) %>%
  slice(1:(n() - 1)) %>%   # drop last point (no "next" point)
  transmute(
    from        = Name,
    to          = lead(Name),
    distance_m  = round(dists_m, 1),
    distance_km = round(dists_m / 1000, 4)
  ) %>%
  st_drop_geometry() %>%
  full_join(all_new_owen, by = join_by(from == Name))

# 3.4 Compute slope for each segment using elevations and distances
slope_table <- distance_table %>%
  arrange(as.numeric(from)) %>%
  mutate(
    elev_this = elev_m,
    elev_next = lead(elev_m),
    # This assumes point order is chosen so this difference is positive.
    slope_m_m = (elev_next - elev_this) / distance_m
  ) %>%
  # Rename PopupInfo (from Google Earth) to stream_width
  rename(stream_width = PopupInfo) %>%
  # Keep only fields we actually need downstream
  select(from, to, distance_m, elev_this, elev_next,
         slope_m_m, stream_width, ws_cum_up_km2)


# ------------------------------------------------------------
# 4. REGIONAL Q–A RELATION AND CALIBRATION DEPTH
# ------------------------------------------------------------

# 4.1 Load regional area–discharge data (area in km2, flow in m3/s)
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

# 4.2 Fit regional Q = a * A^b  via log–log regression
#     ln(Q) = ln(a) + b * ln(A)
fit_QA <- lm(log(Q_m3s) ~ log(area_km2), data = qa_regional_clean)
coef_QA <- coef(fit_QA)

a <- exp(coef_QA[1])  # a (depends on units: Q in m3/s, A in km2 here)
b <- coef_QA[2]       # b (dimensionless)

message("Regional Q–A fit for Owen: Q = ", round(a, 3), " * A^", round(b, 3))

# 4.3 Hydraulic constants for resistance equation
k_s <- 0.25   # roughness height (m) – assumption, can be tuned
g   <- 9.81  # gravity (m/s^2)

# 4.4 Use Q–A to estimate Q and unit discharge q at each calibration cross-section
owen_w_discharge <- slope_table %>%
  mutate(
    stream_width = as.numeric(stream_width), # convert width to numeric
    Q            = a * ws_cum_up_km2^b,      # discharge in m3/s
    q            = Q / stream_width          # unit discharge (m2/s)
  )

# 4.5 Compute depth at calibration sections using the resistance equation from notes:
#     h = [ (k_s^(1/3) * q^2) / (64 * g * S) ]^(3/10)
#     where:
#       h = flow depth (m)
#       q = unit discharge (m2/s)
#       S = slope (m/m)
#
#     Note: here we do a small safety check to avoid S <= 0.
owen_hydro <- owen_w_discharge %>%
  mutate(
    S_use  = ifelse(slope_m_m > 0, slope_m_m, NA_real_),   # avoid zero/negative slopes
    height = (
      (k_s^(1/3) * q^2) /
        (64 * g * S_use)
    )^(3/10)
  )

# 4.6 Fit hydraulic geometry: h = alpha * A^beta
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



owen_hg_data <- owen_hydro %>%
  filter(
    !is.na(ws_cum_up_km2),
    !is.na(height),
    ws_cum_up_km2 > 0,
    height > 0
  )

r2_text <- round(summary(hg_fit)$r.squared, 3)

owen_hg_plot <- ggplot(owen_hg_data,
                       aes(x = ws_cum_up_km2, y = height)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Drainage area (km², log scale)",
    y = "Depth (m, log scale)",
    title = paste("Hydraulic geometry: h = αA^β (R² =", r2_text, ")")
  ) + 
  theme_bw()

print(owen_hg_plot)



message("Hydraulic geometry: h = ",
        round(alpha, 3), " * A^", round(beta, 3),
        "   (A in km2, h in m)")

# 4.7 Save calibration table for inspection / QC
write.csv(owen_hydro, sprintf("tabular data/owen_crk/owen_crk_hydro_output_%s.csv", gsub("\\.", "", as.character(k_s))),
          row.names = FALSE)


# ------------------------------------------------------------
# 5. NEW DATA: AOI CROSS-SECTIONS (DENSE THALWEG POINTS)
# ------------------------------------------------------------

# GOAL: Use alpha, beta + DEM slope + drainage area to estimate h and D50
# at a denser set of AOI cross-section points (owen_xsec_AOI.shp).

# # 5.1 Read AOI points with FID preserved
# watershed_area_pnts <- st_read(
#   "spatial data/owen_crk/owen_xsec_AOI.shp",
#   fid_column_name = "FID"
# )
# 
# # 5.2 Optional: filter every other point if you only want odd FIDs
# watershed_area_pnts <- watershed_area_pnts %>%
#   filter(as.numeric(FID) %% 2 == 1) %>%
#   select(FID, geometry) 
# 
# # 5.3 Save a cleaned AOI shapefile (for reproducibility)
# st_write(
#   watershed_area_pnts,
#   "spatial data/owen_crk/watershed_area_pnts_clean2.shp",
#   delete_dsn = TRUE
# )

# 5.4 Re-read cleaned AOI points (now the canonical input for watersheds)
watershed_area_pnts <- st_read(
  "spatial data/owen_crk/watershed_area_pnts_clean5.shp"
)

# 5.5 Use make_watershed to delineate drainage areas for AOI points
watershed_AOI <- make_watershed(
  path_start   = path_start,
  dem_utm_path = "spatial data/owen_crk/dem_utm.tif",
  dem_utm      = dem_utm,
  pour_points  = "spatial data/owen_crk/watershed_area_pnts_clean5.shp"
)

# SIDE QUEST -----------------------


pts_path   <- file.path(path_start, "pour_points_snapped_raw.shp")
poly_path  <- file.path(path_start, "watersheds_polygons_o.gpkg")

out_path   <- file.path(path_start, "pour_points_with_ws_meta.gpkg")


# 1. Read points and polygons --------------------------------

# Points: snapped pour points + drainage area from FAC
pts <- st_read(pts_path, quiet = TRUE)

# Polygons: watershed polygons + local area
polys <- st_read(poly_path, quiet = TRUE)

# Check CRS and align if needed
if (st_crs(pts) != st_crs(polys)) {
  polys <- st_transform(polys, st_crs(pts))
}


# 2. Choose which polygon attributes to copy to points -------

# Keep just the watershed-related attributes you care about
# (add/remove fields as needed)
polys_keep <- polys %>%
  dplyr::select(
    watersheds,
    ws_cum_up_km2
  )


# 3. Spatial join: point-in-polygon (INTERSECT) ---------------

# Each point gets the attributes of the polygon it falls inside
pts_joined <- st_join(
  pts,
  polys_keep,
  join = st_intersects,  # point-in-polygon
  left = TRUE            # keep all points even if no polygon match
)



# END QUEST -------------------------

# ------------------------------------------------------------
# 5. SPATIAL JOIN: ATTACH POLYGON (WATERSHED) ATTRIBUTES TO AOI POINTS
# ------------------------------------------------------------

# 5.1 Make sure polygons and points are in the same CRS
watershed_AOI <- st_transform(watershed_AOI, st_crs(watershed_area_pnts))

# 5.2 Keep only the polygon fields you actually need on the points
#     (add/remove fields from this select() as needed)
polys_keep <- watershed_AOI %>%
  dplyr::select(
    ws_cum_up_km2,   # drainage area (km2) or whatever field you want for A
    watersheds       # polygon ID (optional but nice to keep)
  )

# 5.3 Spatial join: each AOI point gets the attributes of the polygon it falls in
thalweg_o <- st_join(
  watershed_area_pnts,
  polys_keep,
  join = st_intersects,  # point-in-polygon
  left = TRUE            # keep all points even if no polygon match
)

# At this point, thalweg_o has:
#  - original point attributes (e.g. FID, etc.)
#  - ws_cum_up_km2 from polygons


# ------------------------------------------------------------
# 6. BUILD THALWEG OBJECT WITH DRAINAGE AREA, ELEVATION, SLOPE
# ------------------------------------------------------------

# 6.1 Put thalweg points into UTM CRS for distances + DEM_UTM sampling
thalweg_o <- thalweg_o %>%
  st_transform(st_crs(dem_utm)) %>%
  st_zm(drop = TRUE, what = "ZM")   # ensure 2D

# 6.2 Extract elevation from UTM DEM
zvals_AOI <- terra::extract(dem_utm, vect(thalweg_o))
thalweg_o$elev_m <- zvals_AOI[, 2]

# 6.3 Compute along-channel distances between successive points
#     Here we use Euclidean distances in UTM coordinates (meters).
#     Sort points along-channel by whatever index makes sense; if FID is ordered,
#     use that. Otherwise, replace FID with your own ordering field.
thalweg_o <- thalweg_o %>%
  arrange(FID)   # or some other along-channel index

coords_o <- st_coordinates(thalweg_o)

dx_vec <- sqrt(diff(coords_o[, 1])^2 + diff(coords_o[, 2])^2)

# 6.4 Add dx (segment length) and slope S_seg using dz = elev_next - elev_this
thalweg_o <- thalweg_o %>%
  mutate(
    dx        = c(dx_vec, NA_real_),  # last point has no downstream neighbor
    elev_this = elev_m,
    elev_next = lead(elev_m),
    dz        = elev_next - elev_this,
    S_seg     = dz / dx               # raw local segment slope
  )

# 6.5 Smooth slope and keep magnitude
k <- 30  # window size; tweak based on point spacing
thalweg_o_sloped <- thalweg_o %>%
  mutate(
    S_seg_abs = abs(S_seg),
    S_roll    = zoo::rollapply(
      S_seg_abs,
      width = k,
      FUN   = mean,
      fill  = NA,
      align = "center"
    ),
    S_use     = ifelse(is.na(S_roll), S_seg_abs, S_roll)  # fallback to raw |slope|
  )


# ------------------------------------------------------------
# 7. COMPUTE D50 USING HYDRAULIC GEOMETRY + SHIELDS EQUATION
# ------------------------------------------------------------

rho   <- 1000   # water density (kg/m3)
rhos  <- 2650   # sediment density (kg/m3)
tau_c <- 0.045  # critical Shields (dimensionless)

# h = alpha * A^beta   (A = ws_cum_up_km2, h in m)
# τ = ρ g h S
# D50 = τ / ((ρs - ρ) g τ_c) = (ρ h S) / ((ρs - ρ) τ_c)

thalweg_o_final <- thalweg_o_sloped %>%
  mutate(
    h_m   = alpha * (ws_cum_up_km2^beta),   # bankfull depth from hydraulic geometry
    D50_m = (rho * h_m * S_use) / ((rhos - rho) * tau_c),
    D50_mm = D50_m * 1000,
    spawn_class = case_when(
      is.na(D50_mm)        ~ NA_character_,
      D50_mm < 7           ~ "too fine",
      D50_mm <= 42         ~ "spawning gravel",
      TRUE                 ~ "too coarse"
    )
  ) %>%
  arrange(as.numeric(FID)) %>%
  rename(FID_ORIG = FID)

thalweg_o_final_df <- thalweg_o_final %>%
  st_drop_geometry() 

thalweg_o_final_df_count <- thalweg_o_final_df %>%
  count(spawn_class) %>%
  mutate(percent = round(100*n/sum(n), 2))
  





summarize_d50(thalweg_o_final_df)


owen <- thalweg_o_final %>%
  filter(!is.na(D50_mm)) %>%
  arrange(D50_mm) %>%
  mutate(
    cumulative = cumsum(rep(1, n())) / n()
  )

n_owen <- nrow(owen)

# ---- 2. Plot S-curve ----
ggplot(owen, aes(x = D50_mm, y = cumulative)) +
  geom_line(color = "#2C7BB6", size = 1.2) +
  
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(trans = "log10") +
  
  # Gravel thresholds (Mamquam = 7–47 mm)
  geom_vline(xintercept = 7,  color = "grey40", linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = 47, color = "firebrick", linetype = "dotted", linewidth = 0.8) +
  
  labs(
    title = "Owen Creek – Cumulative Grain Size Distribution",
    #subtitle = paste0("n = ", n_mamq, " predicted bed-material samples"),
    x = "Predicted median grain size D50 (mm, log10 scale)",
    y = "% finer"
  ) + theme_bw()

# ------------------------------------------------------------
# 8. SAVE OUTPUT FOR MAPPING / FURTHER ANALYSIS
# ------------------------------------------------------------

write.csv(thalweg_o_final_df, sprintf("tabular data/owen_crk/owen_properties_output_k%s.csv", gsub("\\.", "", as.character(k_s))))
write.csv(thalweg_o_final_df_count, sprintf("tabular data/owen_crk/owen_properties_output_proportion_k%s.csv", gsub("\\.", "", as.character(k_s))))


st_write(
  thalweg_o_final,
  file.path(path_start, sprintf("owen_thalweg_o_D50_%s.gpkg", gsub("\\.", "", as.character(k_s)))),
  delete_dsn = TRUE
)
