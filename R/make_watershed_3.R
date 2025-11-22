library(sf)
library(terra)
library(whitebox)
library(dplyr)
library(mapview)

# path_start: directory where intermediates + outputs will live (must exist)
# dem_utm_path: path to DEM in same CRS as dem_utm (e.g., "…/dem_utm.tif")
# dem_utm: terra rast object of same DEM (used only for mapview/QC)
# pour_points: path to *POINT* layer with pour points (XS locations etc.)
make_watershed <- function(path_start, dem_utm_path, dem_utm, pour_points) {
  
  # ---------------------------------------------------------------------------
  # 1. Pre-process DEM: breach, flow direction, flow accumulation (catchment area)
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
  
  r_fac      <- rast(file.path(path_start, "fac_area.tif"))
  ws_pntr    <- rast(file.path(path_start, "d8_pntr.tif"))
  target_crs <- crs(r_fac)
  
  # ---------------------------------------------------------------------------
  # 2. Read pour points, force them to 2D POINT in raster CRS
  # ---------------------------------------------------------------------------
  message("read + standardize pour points")
  pp_raw <- tryCatch(
    st_read(pour_points, quiet = TRUE, fid_column_name = "FID"),
    error = function(e) st_read(pour_points, quiet = TRUE)
  )
  
  # Drop Z/M and cast to POINT to keep Whitebox happy
  pp_ras <- pp_raw %>%
    st_zm(drop = TRUE, what = "ZM") %>%
    st_cast("POINT") %>%
    st_transform(st_crs(target_crs))
  
  # Whitebox currently prefers shapefile input; we treat this as temporary
  pp_tmp <- file.path(path_start, "pour_points_in_raster_crs.shp")
  st_write(pp_ras, pp_tmp, delete_dsn = TRUE, quiet = TRUE)
  
  # ---------------------------------------------------------------------------
  # 3. Snap pour points to high-FAC cells and extract cumulative drainage area
  # ---------------------------------------------------------------------------
  message("snap pour points")
  wbt_snap_pour_points(
    pour_pts   = pp_tmp,
    flow_accum = file.path(path_start, "fac_area.tif"),
    output     = file.path(path_start, "pour_points_snapped.shp"),
    snap_dist  = 100
  )
  
  pour_snapped <- st_read(file.path(path_start, "pour_points_snapped.shp"),
                          quiet = TRUE) %>%
    st_zm(drop = TRUE, what = "ZM") %>%
    st_transform(st_crs(target_crs))
  
  # Extract FAC (catchment area) at each snapped point = TOTAL upstream area
  da_vals <- terra::extract(r_fac, vect(pour_snapped))
  # 2nd col = FAC value
  pour_snapped$drainage_area_m2  <- da_vals[, 2]
  pour_snapped$drainage_area_km2 <- pour_snapped$drainage_area_m2 / 1e6
  
  # Save points + drainage area as GPKG for analysis & mapping
  st_write(
    pour_snapped,
    file.path(path_start, "thalwag_with_drainage_area.gpkg"),
    delete_dsn = TRUE,
    quiet = TRUE
  )
  
  # ---------------------------------------------------------------------------
  # 4. Delineate watersheds (raster) and convert to polygons
  # ---------------------------------------------------------------------------
  message("watershed")
  wbt_watershed(
    d8_pntr  = file.path(path_start, "d8_pntr.tif"),
    pour_pts = file.path(path_start, "pour_points_snapped.shp"),
    output   = file.path(path_start, "watersheds.tif")
  )
  
  ws_rast <- rast(file.path(path_start, "watersheds.tif"))
  
  # Polygonize watershed raster; 'watersheds' is the ID field
  ws_poly <- as.polygons(ws_rast, values = TRUE, dissolve = TRUE, na.rm = TRUE)
  ws_poly_sf <- st_as_sf(ws_poly) %>%
    st_make_valid() %>%
    st_transform(st_crs(target_crs)) %>%
    st_zm(drop = TRUE, what = "ZM")
  
  if (nrow(ws_poly_sf) == 0) {
    stop("Watershed polygonization produced 0 features.")
  }
  
  # Local polygon area (geometric area of each watershed patch)
  ws_poly_sf <- ws_poly_sf %>%
    mutate(
      ws_local_area_m2  = as.numeric(st_area(geometry)),
      ws_local_area_km2 = ws_local_area_m2 / 1e6
    )
  
  # ---------------------------------------------------------------------------
  # 5. Attach cumulative upstream area from FAC to polygons by watershed ID
  # ---------------------------------------------------------------------------
  message("attach upstream area to polygons")
  
  # Extract watershed ID at each pour point from the watershed raster
  ws_id_vals <- terra::extract(ws_rast, vect(pour_snapped))[, 2]
  pour_snapped$ws_id <- ws_id_vals
  
  # Build lookup table: one row per watershed ID with its cumulative area
  da_lookup <- pour_snapped %>%
    st_drop_geometry() %>%
    select(ws_id, drainage_area_m2, drainage_area_km2) %>%
    distinct(ws_id, .keep_all = TRUE)
  
  # Attach to polygons using watershed ID
  ws_poly_sf <- ws_poly_sf %>%
    rename(ws_id = watersheds) %>%
    left_join(da_lookup, by = "ws_id") %>%
    rename(
      ws_upstream_area_m2  = drainage_area_m2,
      ws_upstream_area_km2 = drainage_area_km2
    )
  
  # Write polygons to GPKG
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
  
  # Return SF polygons with both local + upstream area
  return(ws_poly_sf)
}
