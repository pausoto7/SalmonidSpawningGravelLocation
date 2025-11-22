# library(sf)
# library(terra)
# library(whitebox)   # if you’re using whiteboxr
# library(dplyr)
# library(mapview)
# 
# make_watershed <- function(path_start, dem_utm_path, dem_utm, pour_points){
#   
#   message("breach/fill")
#   wbt_breach_depressions_least_cost(
#     dem    = dem_utm_path,
#     output = file.path(path_start, "dem_breached.tif"),
#     dist   = 30
#   )
#   
#   message("pointer")
#   wbt_d8_pointer(
#     dem    = file.path(path_start, "dem_breached.tif"),
#     output = file.path(path_start, "d8_pntr.tif")
#   )
#   
#   message("flow accumulation (catchment area)")
#   wbt_d8_flow_accumulation(
#     input    = file.path(path_start, "d8_pntr.tif"),
#     output   = file.path(path_start, "fac_area.tif"),
#     out_type = "catchment area",
#     pntr     = TRUE
#   )
#   
#   # ---- FIX: adopt the raster CRS and apply it to vectors up front ----
#   r_fac  <- rast(file.path(path_start, "fac_area.tif"))
#   r_pntr <- rast(file.path(path_start, "d8_pntr.tif"))
#   target_crs <- crs(r_fac)  # terra CRS string
#   
#   # Read raw pour points, transform to raster CRS, and write a temp shapefile
#   pp_raw <- st_read(pour_points, quiet = TRUE)
#   pp_ras <- st_transform(pp_raw, st_crs(target_crs))
#   pp_path <- file.path(path_start, "pour_points_in_raster_crs.shp")
#   st_write(pp_ras, pp_path, delete_dsn = TRUE, quiet = TRUE)
#   
#   message("snap pour points")
#   wbt_snap_pour_points(
#     pour_pts   = pp_path,                               # <- in raster CRS
#     flow_accum = file.path(path_start, "fac_area.tif"),
#     output     = file.path(path_start, "pour_points_snapped.shp"),
#     snap_dist  = 100
#   )
#   
#   pour_snapped <- st_read(file.path(path_start, "pour_points_snapped.shp"), quiet = TRUE)
#   # ---- FIX: ensure snapped points are also exactly in raster CRS (2D) ----
#   pour_snapped <- st_transform(pour_snapped, st_crs(target_crs))
#   pour_snapped <- st_zm(pour_snapped, drop = TRUE, what = "ZM")  # drop Z/M if present
#   
#   # Attach drainage area (m^2) from FAC to the snapped points
#   da_vals <- terra::extract(r_fac, vect(pour_snapped))
#   pour_snapped$drainage_area_m2  <- da_vals[, 2]
#   pour_snapped$drainage_area_km2 <- pour_snapped$drainage_area_m2 / 1e6
#   
#   st_write(pour_snapped, file.path(path_start, "thalwag_with_drainage_area.gpkg"),
#            delete_dsn = TRUE, quiet = TRUE)
#   
#   message("watershed")
#   wbt_watershed(
#     d8_pntr  = file.path(path_start, "d8_pntr.tif"),
#     pour_pts = file.path(path_start, "pour_points_snapped.shp"),
#     output   = file.path(path_start, "watersheds.tif")
#   )
#   
#   ws_rast <- rast(file.path(path_start, "watersheds.tif"))
#   
#   # ---- FIX: drop NoData and polygonize robustly ----
#   # (terra::as.polygons handles NA if mask=TRUE or we clamp first)
#   ws_rast <- classify(ws_rast, cbind(NA, NA, NA))  # ensure NA remains NA
#   ws_poly <- as.polygons(ws_rast, values = TRUE, dissolve = TRUE, na.rm = TRUE)
#   ws_poly_sf <- st_as_sf(ws_poly)
#   ws_poly_sf <- st_make_valid(ws_poly_sf)
#   ws_poly_sf <- st_transform(ws_poly_sf, st_crs(target_crs))
#   ws_poly_sf <- st_zm(ws_poly_sf, drop = TRUE, what = "ZM")
#   
#   # Guard: ensure we actually have polygons
#   if (nrow(ws_poly_sf) == 0) stop("Watershed polygonization produced 0 features.")
#   
#   # Add polygon areas (units from CRS; UTM -> m^2)
#   ws_poly_sf <- ws_poly_sf %>%
#     mutate(
#       ws_area_m2  = as.numeric(st_area(geometry)),
#       ws_area_km2 = ws_area_m2 / 1e6
#     )
#   
#   st_write(ws_poly_sf, file.path(path_start, "watersheds_polygons.shp"),
#            delete_dsn = TRUE, quiet = TRUE)
#   
#   # ---- FIX: make CRS identical before spatial join & use intersects ----
#   pour_snapped <- st_transform(pour_snapped, st_crs(ws_poly_sf))
#   ws_area_only <- ws_poly_sf %>% select(ws_area_m2, ws_area_km2)
#   
#   # Use st_intersects to be robust to points on boundaries; pick first match
#   idx <- st_intersects(pour_snapped, ws_area_only)
#   got_match <- lengths(idx) > 0
#   join_rows <- ifelse(got_match, vapply(idx, `[`, integer(1), 1), NA_integer_)
#   
#   pour_with_da <- bind_cols(
#     pour_snapped,
#     ws_area_only[join_rows, c("ws_area_m2", "ws_area_km2")] |> st_drop_geometry()
#   )
#   
#   # Write out joined points (overwriting earlier GPKG on purpose)
#   st_write(pour_with_da, file.path(path_start, "thalwag_with_drainage_area.gpkg"),
#            delete_dsn = TRUE, quiet = TRUE)
#   
#   # Quick map (optional)
#   suppressWarnings({
#     mapview(dem_utm) +
#       mapview(ws_poly_sf, alpha.regions = 0.2, col.regions = NA, color = "black") +
#       mapview(pour_with_da, zcol = "ws_area_km2")
#   })
#   
#   return(ws_poly_sf)
# }
