library(sf)
library(terra)
library(whitebox)
library(dplyr)
library(mapview)

make_watershed <- function(path_start, dem_utm_path, dem_utm, pour_points) {
  
  # ---------------------------------------------------------------------------
  # 1. Pre-process DEM: breach, pointer, flow accumulation (catchment area)
  # ---------------------------------------------------------------------------
  message("breach/fill")
  wbt_breach_depressions_least_cost(
    dem    = dem_utm_path,
    output = file.path(path_start, "dem_breached.tif"),
    dist   = 30
  )
  
  message("pointer")
  wbt_d8_pointer(
    dem    = file.path(path_start, "dem_breached.tif"),
    output = file.path(path_start, "d8_pntr.tif")
  )
  
  message("flow accumulation (catchment area)")
  wbt_d8_flow_accumulation(
    input    = file.path(path_start, "d8_pntr.tif"),
    output   = file.path(path_start, "fac_area.tif"),
    out_type = "catchment area",
    pntr     = TRUE
  )
  
  r_fac       <- rast(file.path(path_start, "fac_area.tif"))
  target_crs  <- crs(r_fac)
  
  # ---------------------------------------------------------------------------
  # 2. Read pour points, transform to raster CRS, keep FID if available
  # ---------------------------------------------------------------------------
  pp_raw <- tryCatch(
    st_read(pour_points, quiet = TRUE, fid_column_name = "FID"),
    error = function(e) st_read(pour_points, quiet = TRUE)
  )
  
  pp_ras <- pp_raw %>%
    st_transform(st_crs(target_crs)) %>%
    st_zm(drop = TRUE, what = "ZM")
  
  pp_path <- file.path(path_start, "pour_points_in_raster_crs.shp")
  st_write(pp_ras, pp_path, delete_dsn = TRUE, quiet = TRUE)
  
  # ---------------------------------------------------------------------------
  # 3. Snap pour points and attach cumulative upstream drainage area
  # ---------------------------------------------------------------------------
  message("snap pour points")
  wbt_snap_pour_points(
    pour_pts   = pp_path,
    flow_accum = file.path(path_start, "fac_area.tif"),
    output     = file.path(path_start, "pour_points_snapped.shp"),
    snap_dist  = 100
  )
  
  pour_snapped <- st_read(file.path(path_start, "pour_points_snapped.shp"),
                          quiet = TRUE) %>%
    st_transform(st_crs(target_crs)) %>%
    st_zm(drop = TRUE, what = "ZM")
  
  # extract FAC value at each snapped point: this is TOTAL upstream area
  da_vals <- terra::extract(r_fac, vect(pour_snapped))
  pour_snapped$drainage_area_m2  <- da_vals[, 2]
  pour_snapped$drainage_area_km2 <- pour_snapped$drainage_area_m2 / 1e6
  
  st_write(
    pour_snapped,
    file.path(path_start, "POI_with_drainage_area.shp"),
    delete_layer = TRUE,
    quiet = TRUE
  )
  
  # ---------------------------------------------------------------------------
  # 4. Delineate watersheds to polygons
  # ---------------------------------------------------------------------------
  message("watershed")
  wbt_watershed(
    d8_pntr  = file.path(path_start, "d8_pntr.tif"),
    pour_pts = file.path(path_start, "pour_points_snapped.shp"),
    output   = file.path(path_start, "watersheds.tif")
  )
  
  ws_rast <- rast(file.path(path_start, "watersheds.tif"))
  
  ws_poly <- as.polygons(ws_rast, values = TRUE, dissolve = TRUE, na.rm = TRUE)
  ws_poly_sf <- st_as_sf(ws_poly) %>%
    st_make_valid() %>%
    st_transform(st_crs(target_crs)) %>%
    st_zm(drop = TRUE, what = "ZM")
  
  if (nrow(ws_poly_sf) == 0) {
    stop("Watershed polygonization produced 0 features.")
  }
  
  # local polygon area (for reference / QC)
  ws_poly_sf <- ws_poly_sf %>%
    mutate(
      ws_loc_A_m2  = as.numeric(st_area(geometry)),
      ws_loc_A_km2 = ws_loc_A_m2 / 1e6
    )
  
  # ---------------------------------------------------------------------------
  # 5. Copy cumulative upstream area from pour_snapped into polygons
  #    (each polygon should have exactly one pour point)
  # ---------------------------------------------------------------------------
  pour_snapped <- st_transform(pour_snapped, st_crs(ws_poly_sf))
  
  # spatial join: for each polygon, attach attributes of point it contains
  ws_poly_sf <- st_join(
    ws_poly_sf,
    pour_snapped %>% select(drainage_area_m2, drainage_area_km2),
    join = st_contains,  # polygon contains point
    left = TRUE
  )
  
  ws_poly_sf <- ws_poly_sf %>%
    rename(
      ws_up_A_m2  = drainage_area_m2,
      ws_up_A_km2 = drainage_area_km2
    )
  
  

  st_write(
    ws_poly_sf,
    file.path(path_start, "watersheds_polygons.gpkg"),
    delete_dsn = TRUE,
    quiet = TRUE
  )
  
  
  # ---------------------------------------------------------------------------
  # 6. Quick visual QC (optional)
  # ---------------------------------------------------------------------------
  suppressWarnings({
    mapview(dem_utm) +
      mapview(ws_poly_sf, alpha.regions = 0.2,
              col.regions = NA, color = "black") +
      mapview(pour_snapped, zcol = "drainage_area_km2")
  })
  
  return(ws_poly_sf)
}
