
source("R/make_watershed_3.R")

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

# a <- 23.753
# b <- 0.7808
a <- 0.228
b <- 1.109
k_s <- 0.3 #roughness coefficient
g   <- 9.81          # gravity, m/s^2

  
owen_w_discharge <- slope_table %>%
  mutate(stream_width = as.numeric(stream_width), 
         Q = a*ws_cum_up_km2^b, 
         q = Q/stream_width)

owen_hydro <- owen_w_discharge %>%
  mutate(
    height = (
      (k_s^(1/3) * q^2 * stream_width) /
        (64 * g * slope_m_m)
    )^(3/10)
  )

fit_hA <- lm(log(height) ~ log(ws_cum_up_km2), data = owen_hydro)
coef_hA <- coef(fit_hA)
alpha_h <- exp(coef_hA[1])
beta_h  <- coef_hA[2]




library(ggplot2)

ggplot(owen_hydro, aes(x = ws_cum_up_km2, y = height)) +
  geom_point() +
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Drainage area (km²)", y = "Depth h (m)",
       title = "h–A hydraulic geometry (Owen)")

library(ggrepel)

ggplot(owen_hydro, aes(x = ws_cum_up_km2, y = stream_width)) +
  geom_point() +
  geom_text_repel(aes(label = from), size = 3) +
  
  scale_x_log10() +
  scale_y_log10() +
  labs(x = "Drainage area (km²)", y = "Width w (m)",
       title = "w–A hydraulic geometry (Owen)")



