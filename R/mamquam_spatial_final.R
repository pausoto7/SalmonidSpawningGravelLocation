# ============================================================
# PREDICTING D50 FOR MAMQUAM RIVER USING DEM + HYDRAULIC GEOMETRY
# ============================================================

# This script does two main things:
#  1) Uses calibration points with known widths + DEM + regional Q–A
#     to estimate depth via a resistance equation, then fit
#     hydraulic geometry h = alpha * A^beta.
#  2) Applies that hydraulic geometry to a denser set of AOI
#     cross-section points to compute D50 and classify spawning habitat.
#
# Requires:
#   - R/make_watershed_3.R (your existing function)
#   - Mamquam dem_utm.tif (UTM DEM)
#   - Mamquam calibration pour points (with width field)
#   - Mamquam AOI thalweg points for dense D50 mapping
#
# ------------------------------------------------------------
# 0. SETUP
# ------------------------------------------------------------

source("R/make_watershed_3.R")

library(tidyverse)
library(terra)       # rasters
library(sf)          # vectors / shapefiles
library(zoo)         # rollapply for smoothing slopes
library(whitebox)    # watershed / hydrologic tools
library(mapview)

# ------------------------------------------------------------
# 0.1 USER CONFIG: CHECK / EDIT THESE PATHS + FIELD NAMES
# ------------------------------------------------------------

# Base path for Mamquam outputs
path_start <- "spatial data/mamquam_river/new_workflow"

# UTM DEM (this one you *have* used before)
dem_utm_path <- "spatial data/mamquam_river/dem_utm.tif"

# Calibration pour points (Mamquam). From earlier code:
#   "spatial data/mamquam_river/mamq_pour_points_AOI.shp"
mamq_calib_points_path <- "spatial data/mamquam_river/mamquam_with_widths.shp"

# AOI / dense thalweg points for Mamquam (UPDATE to your actual AOI file)
# e.g., "spatial data/mamquam_river/mamq_xsec_AOI.shp" or similar
mamq_aoi_points_path <- "spatial data/mamquam_river/mamq_100m_pnts_aoi.shp"  
# ^ if you have a denser AOI file, replace this with that path.

# Regional Q–A CSV for Mamquam (UPDATE this to your actual file)
mamq_qa_csv_path <- "spatial data/mamquam_river/mamquam_similar_watershed.csv"
# (leave as-is if you haven’t created it yet; script will error if the file
#  doesn’t exist, which is better than silently using a fake path)

# Column in calibration pour points that holds bankfull width
# For Owen you used PopupInfo (from Google Earth). If Mamquam uses a different
# name (e.g., "width_m"), change it here.
width_col <- "PopupInfo"

# Roughness height for resistance equation (m)
k_s <- 0.25
# Slope smoothing window (number of points)
k_slope <- 3

# Critical Shields parameter (dimensionless)
tau_c <- 0.045

# Output tabular path (make sure this folder exists)
tabular_path <- "tabular data/mamquam_river"


# ------------------------------------------------------------
# 1. LOAD DEM AND BUILD WATERSHEDS FOR CALIBRATION POINTS
# ------------------------------------------------------------

# 1.1 Load UTM DEM directly
dem_utm <- rast(dem_utm_path)

# 1.2 Run custom watershed function for Mamquam calibration pour points
watershed_polygon_mamq <- make_watershed(
  path_start   = path_start,
  dem_utm_path = dem_utm_path,
  dem_utm      = dem_utm,
  pour_points  = mamq_calib_points_path, 
  poly_name = "SAMPLEDATA"
)

# 1.3 Drop geometry to retain polygon attributes (watersheds, ws_cum_up_km2, etc.)
drainage_df_mamq <- st_drop_geometry(watershed_polygon_mamq)


# ------------------------------------------------------------
# 2. CALIBRATION POINTS: ADD ELEVATION AND DISTANCES
# ------------------------------------------------------------

# 2.1 Load calibration pour points
mamq_pp <- st_read(mamq_calib_points_path)

# 2.2 Reproject to UTM DEM CRS (for elevation + distance)
mamq_pp_utm <- mamq_pp %>%
  st_transform(st_crs(dem_utm)) %>%
  st_zm(drop = TRUE, what = "ZM")

# 2.3 Extract elevation from UTM DEM at each pour point
zvals_mamq <- terra::extract(dem_utm, vect(mamq_pp_utm))
mamq_pp_utm$elev_m <- zvals_mamq[, 2]

# 2.4 Join pour points to watershed attributes
# NOTE: THIS ASSUMES your pour point layer has an ID field that matches
#       "watersheds" from make_watershed (for Owen it was Name == watersheds).
#       If your Mamquam points use a different ID field, change "Name" below.
mamq_all <- mamq_pp_utm %>%
  mutate(Name = as.integer(Name)) %>%  # << change 'Name' if needed
  right_join(drainage_df_mamq, by = join_by(Name == watersheds))%>%
  filter(!is.na(FolderPath)) %>%
  arrange(as.numeric(Name))
  


# ------------------------------------------------------------
# 3. DISTANCES AND SLOPES AT CALIBRATION POINTS
# ------------------------------------------------------------

# 3.1 Extract coordinates in Name order (assumes Name follows channel order)
coords_mamq <- mamq_all %>%
  st_coordinates()

# 3.2 Compute straight-line distances between successive points (in meters)
dx_m <- sqrt(diff(coords_mamq[, 1])^2 + diff(coords_mamq[, 2])^2)

# 3.3 Build distance + attribute table
mamq_distance_table <- mamq_all %>%
  arrange(as.numeric(Name)) %>%
  slice(1:(n() - 1)) %>%   # last point has no downstream neighbor
  transmute(
    from        = Name,
    to          = lead(Name),
    distance_m  = round(dx_m, 1),
    distance_km = round(dx_m / 1000, 4)
  ) %>%
  st_drop_geometry() %>%
  full_join(mamq_all, by = join_by(from == Name))

# 3.4 Compute slope using elevations and distances
mamq_slope_table <- mamq_distance_table %>%
  arrange(as.numeric(from)) %>%
  mutate(
    elev_this = elev_m,
    elev_next = lead(elev_m),
    slope_m_m = (elev_next - elev_this) / distance_m
  ) %>%
  # rename width field to stream_width
  rename(stream_width = !!sym(width_col)) %>%
  select(from, to, distance_m, elev_this, elev_next,
         slope_m_m, stream_width, ws_cum_up_km2)


# ------------------------------------------------------------
# 4. REGIONAL Q–A RELATION AND CALIBRATION DEPTH
# ------------------------------------------------------------

# 4.1 Load regional area–discharge data for Mamquam
mamq_qa_regional <- read.csv(mamq_qa_csv_path, stringsAsFactors = FALSE)

mamq_qa_clean <- mamq_qa_regional %>%
  dplyr::rename(
    area_km2 = area,
    Q_m3s    = flow
  ) %>%
  dplyr::filter(
    !is.na(area_km2), !is.na(Q_m3s),
    area_km2 > 0, Q_m3s > 0
  )

# 4.2 Fit Q = a * A^b
mamq_fit_QA <- lm(log(Q_m3s) ~ log(area_km2), data = mamq_qa_clean)
mamq_coef_QA <- coef(mamq_fit_QA)

a_mamq <- exp(mamq_coef_QA[1])
b_mamq <- mamq_coef_QA[2]

message("Regional Q–A fit for Mamquam: Q = ",
        round(a_mamq, 3), " * A^", round(b_mamq, 3))

g <- 9.81  # gravity

# 4.3 Estimate Q and unit discharge q at each calibration cross-section
mamq_w_discharge <- mamq_slope_table %>%
  mutate(
    stream_width = as.numeric(stream_width),
    Q            = a_mamq * ws_cum_up_km2^b_mamq,
    q            = Q / stream_width
  )

# 4.4 Compute depth via resistance equation:
#     h = [ (k_s^(1/3) * q^2) / (64 * g * S) ]^(3/10)
mamq_hydro <- mamq_w_discharge %>%
  mutate(
    S_use  = ifelse(slope_m_m > 0, slope_m_m, NA_real_),
    height = (
      (k_s^(1/3) * q^2) /
        (64 * g * S_use)
    )^(3/10)
  )

# 4.5 Fit hydraulic geometry: h = alpha * A^beta
mamq_hg_fit <- mamq_hydro %>%
  filter(
    !is.na(ws_cum_up_km2),
    !is.na(height),
    ws_cum_up_km2 > 0,
    height > 0
  ) %>%
  lm(log(height) ~ log(ws_cum_up_km2), data = .)

mamq_coef_hg <- coef(mamq_hg_fit)
alpha_mamq   <- exp(mamq_coef_hg[1])
beta_mamq    <- mamq_coef_hg[2]

mamq_hg_data <- mamq_hydro %>%
  filter(
    !is.na(ws_cum_up_km2),
    !is.na(height),
    ws_cum_up_km2 > 0,
    height > 0
  )

r2_text <- round(summary(mamq_hg_fit)$r.squared, 3)

mamq_hg_plot <- ggplot(mamq_hg_data,
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

print(mamq_hg_plot)




message("Mamquam hydraulic geometry: h = ",
        round(alpha_mamq, 3), " * A^", round(beta_mamq, 3),
        "   (A in km2, h in m)")

# 4.6 Save calibration table
if (!dir.exists(tabular_path)) dir.create(tabular_path, recursive = TRUE)

write.csv(
  mamq_hydro,
  file.path(tabular_path,
            sprintf("mamq_hydro_output_ks%s.csv", gsub("\\.", "", as.character(k_s)))),
  row.names = FALSE
)


# ------------------------------------------------------------
# 5. NEW DATA: AOI CROSS-SECTIONS (DENSE THALWEG POINTS)
# ------------------------------------------------------------

# 5.1 Read AOI points (thalweg / X-sections)
mamq_aoi <- st_read(mamq_aoi_points_path, quiet = TRUE,
                    fid_column_name = "FID")

# 5.2 Use make_watershed to delineate drainage areas for AOI points
watershed_AOI_mamq <- make_watershed(
  path_start   = path_start,
  dem_utm_path = dem_utm_path,
  dem_utm      = dem_utm,
  pour_points  = mamq_aoi_points_path, 
  poly_name = "100m"
)

# 5.3 Spatial join: attach ws_cum_up_km2 to AOI points
watershed_AOI_mamq <- st_transform(watershed_AOI_mamq, st_crs(mamq_aoi))

polys_keep_mamq <- watershed_AOI_mamq %>%
  dplyr::select(
    ws_cum_up_km2,
    watersheds
  )

thalweg_m <- st_join(
  mamq_aoi,
  polys_keep_mamq,
  join = st_intersects,
  left = TRUE
)


# ------------------------------------------------------------
# 6. BUILD THALWEG OBJECT WITH DRAINAGE AREA, ELEVATION, SLOPE
# ------------------------------------------------------------

thalweg_m <- thalweg_m %>%
  st_transform(st_crs(dem_utm)) %>%
  st_zm(drop = TRUE, what = "ZM")

zvals_m_AOI <- terra::extract(dem_utm, vect(thalweg_m))
thalweg_m$elev_m <- zvals_m_AOI[, 2]

# Assume FID is your along-channel index; change if needed
thalweg_m <- thalweg_m %>%
  arrange(FID)

coords_m_thal <- st_coordinates(thalweg_m)
dx_vec <- sqrt(diff(coords_m_thal[, 1])^2 + diff(coords_m_thal[, 2])^2)

thalweg_m <- thalweg_m %>%
  mutate(
    dx        = c(dx_vec, NA_real_),
    elev_this = elev_m,
    elev_next = lead(elev_m),
    dz        = elev_next - elev_this,
    S_seg     = dz / dx
  )

# Slope smoothing
thalweg_m_sloped <- thalweg_m %>%
  mutate(
    S_seg_abs = abs(S_seg),
    S_roll    = zoo::rollapply(
      S_seg_abs,
      width = k_slope,
      FUN   = mean,
      fill  = NA,
      align = "center"
    ),
    S_use     = ifelse(is.na(S_roll), S_seg_abs, S_roll)
  )


# ------------------------------------------------------------
# 7. COMPUTE D50 USING HYDRAULIC GEOMETRY + SHIELDS EQUATION
# ------------------------------------------------------------

rho  <- 1000
rhos <- 2650

thalweg_m_final <- thalweg_m_sloped %>%
  mutate(
    h_m   = alpha_mamq * (ws_cum_up_km2^beta_mamq),
    D50_m = (rho * h_m * S_use) / ((rhos - rho) * tau_c),
    D50_mm = D50_m * 1000,
    spawn_class = case_when(
      is.na(D50_mm)        ~ NA_character_,
      D50_mm < 7           ~ "too fine",
      D50_mm <= 47         ~ "spawning gravel",
      TRUE                 ~ "too coarse"
    )
  ) %>%
  arrange(as.numeric(FID)) %>%
  rename(FID_ORIG = FID)

thalweg_m_final_df <- thalweg_m_final %>%
  st_drop_geometry()

thalweg_m_final_df_count <- thalweg_m_final_df %>%
  count(spawn_class) %>%
  mutate(percent = round(100 * n / sum(n), 2))

thalweg_m_final_df_count
# ------------------------------------------------------------
# 8. SAVE OUTPUT FOR MAPPING / FURTHER ANALYSIS
# ------------------------------------------------------------

write.csv(
  thalweg_m_final_df,
  file.path(tabular_path,
            sprintf("mamq_properties_output_ks%s.csv", gsub("\\.", "", as.character(k_s)))),
  row.names = FALSE
)

write.csv(
  thalweg_m_final_df_count,
  file.path(tabular_path,
            sprintf("mamq_properties_output_proportion_ks%s.csv", gsub("\\.", "", as.character(k_s)))),
  row.names = FALSE
)

st_write(
  thalweg_m_final,
  file.path(path_start,
            sprintf("mamq_thalweg_D50_ks%s.gpkg", gsub("\\.", "", as.character(k_s)))),
  delete_dsn = TRUE
)

summarize_d50(thalweg_m_final_df)


mamq <- thalweg_m_final_df %>%
  filter(!is.na(D50_mm)) %>%
  arrange(D50_mm) %>%
  mutate(
    cumulative = cumsum(rep(1, n())) / n()
  )

n_mamq <- nrow(mamq)

# ---- 2. Plot S-curve ----
ggplot(mamq, aes(x = D50_mm, y = cumulative)) +
  geom_line(color = "#2C7BB6", size = 1.2) +
  
  scale_y_continuous(labels = scales::percent_format()) +
  scale_x_continuous(trans = "log10") +
  
  # Gravel thresholds (Mamquam = 7–47 mm)
  geom_vline(xintercept = 7,  color = "grey40", linetype = "dashed", linewidth = 0.8) +
  geom_vline(xintercept = 47, color = "firebrick", linetype = "dotted", linewidth = 0.8) +
  
  labs(
      title = "Mamquam River – Cumulative Grain Size Distribution",
    #subtitle = paste0("n = ", n_mamq, " predicted bed-material samples"),
    x = "Predicted median grain size D50 (mm, log10 scale)",
    y = "% finer"
  ) + theme_bw()
  

summarize_d50 <- function(df) {
  df %>%
    summarise(
      n_points     = n(),
      min_mm       = round(min(D50_mm, na.rm = TRUE), 2),
      max_mm       = round(max(D50_mm, na.rm = TRUE), 2),
      median_mm    = round(median(D50_mm, na.rm = TRUE), 2),
      mean_mm      = round(mean(D50_mm, na.rm = TRUE), 2),
      q25_mm       = round(quantile(D50_mm, 0.25, na.rm = TRUE), 2),
      q75_mm       = round(quantile(D50_mm, 0.75, na.rm = TRUE), 2)
    )
}

