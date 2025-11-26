library(sf)
library(terra)
library(whitebox)
library(dplyr)
library(mapview)

# ---------------------------------------------------------------------------
# make_watershed()
#
# Arguments:
#   path_start   : directory where intermediates + outputs will be written
#   dem_utm_path : path to DEM (in metres, projected; e.g., UTM)
#   dem_utm      : terra::rast of same DEM (used for QC if needed)
#   pour_points  : path to *POINT* layer (AOI / XS / pour points)
#
# Returns:
#   sf polygon layer with:
#     - watersheds          : watershed ID (from raster)
#     - ws_local_area_m2    : geometric polygon area (m2)
#     - ws_local_area_km2   : geometric polygon area (km2)
#     - ws_cum_up_m2        : FAC-based upstream area at pour point (m2)
#     - ws_cum_up_km2       : FAC-based upstream area at pour point (km2)
#
# Side products in path_start:
#   - dem_breached.tif
#   - d8_pntr.tif
#   - fac_area.tif          : D8 catchment area (m2)
#   - pour_points_in_raster_crs.shp
#   - pour_points_snapped.shp
#   - AOI_with_drainage_area.gpkg
#   - watersheds.tif
#   - watersheds_polygons_o.gpkg
# ---------------------------------------------------------------------------

# make_watershed_new <- function(path_start, dem_utm_path, dem_utm, pour_points) {
#   
#   # -------------------------------------------------------------------------
#   # 1. DEM pre-processing: breach, D8 pointer, D8 flow accumulation (FAC)
#   # -------------------------------------------------------------------------
#   message("1) DEM pre-processing: breach / pointer / flow accumulation")
#   
#   message("  - breach_depressions_least_cost")
#   wbt_breach_depressions_least_cost(
#     dem    = dem_utm_path,
#     output = file.path(path_start, "dem_breached.tif"),
#     dist   = 30
#   )
#   
#   message("  - d8_pointer")
#   wbt_d8_pointer(
#     dem    = file.path(path_start, "dem_breached.tif"),
#     output = file.path(path_start, "d8_pntr.tif")
#   )
#   
#   message("  - d8_flow_accumulation (catchment area)")
#   wbt_d8_flow_accumulation(
#     input    = file.path(path_start, "d8_pntr.tif"),
#     output   = file.path(path_start, "fac_area.tif"),
#     out_type = "catchment area",  # outputs area, not cell count
#     pntr     = TRUE
#   )
#   
#   # Load FAC and pointer rasters
#   r_fac      <- rast(file.path(path_start, "fac_area.tif"))   # m2
#   ws_pntr    <- rast(file.path(path_start, "d8_pntr.tif"))
#   target_crs <- crs(r_fac)  # CRS used for all raster-based processing
#   
#   
#   # -------------------------------------------------------------------------
#   # 2. Read pour points, standardize geometry & CRS, write temporary shapefile
#   # -------------------------------------------------------------------------
#   message("2) Read + standardize pour points")
#   
#   # Read original pour points (keep FID if present)
#   pp_raw <- tryCatch(
#     st_read(pour_points, quiet = TRUE, fid_column_name = "FID"),
#     error = function(e) st_read(pour_points, quiet = TRUE)
#   )
#   
#   # Standardize:
#   #   - drop Z/M
#   #   - cast to POINT
#   #   - transform to raster CRS
#   #   - create a simple integer ID (pp_id)
#   #   - select only pp_id so shapefile field names are safe
#   pp_ras <- pp_raw %>%
#     st_zm(drop = TRUE, what = "ZM") %>%
#     st_cast("POINT") %>%
#     st_transform(st_crs(target_crs)) %>%
#     mutate(pp_id = row_number()) %>%
#     dplyr::select(pp_id)  # only one attribute + geometry
#   
#   # Write temporary shapefile in raster CRS for Whitebox
#   pp_tmp <- file.path(path_start, "pour_points_in_raster_crs.shp")
#   st_write(pp_ras, pp_tmp, delete_dsn = TRUE, quiet = TRUE)
#   
#   
#   # -------------------------------------------------------------------------
#   # 3. Snap pour points to high-FAC cells and extract upstream area (FAC)
#   # -------------------------------------------------------------------------
#   message("3) Snap pour points to FAC maxima and extract drainage area")
#   
#   # Snap points to nearby high-flow-accumulation cells
#   wbt_snap_pour_points(
#     pour_pts   = pp_tmp,
#     flow_accum = file.path(path_start, "fac_area.tif"),
#     output     = file.path(path_start, "pour_points_snapped.shp"),
#     snap_dist  = 10  # search radius in map units (m if DEM in metres)
#   )
#   
#   # Read snapped pour points back in
#   pour_snapped <- st_read(
#     file.path(path_start, "pour_points_snapped.shp"),
#     quiet = TRUE
#   ) %>%
#     st_zm(drop = TRUE, what = "ZM") %>%
#     st_transform(st_crs(target_crs))
#   
#   # Extract FAC (catchment area in m2) at each snapped point
#   da_vals <- terra::extract(r_fac, vect(pour_snapped))
#   # 2nd column is FAC value (1st is cell ID)
#   pour_snapped$drainage_area_m2  <- da_vals[, 2]
#   pour_snapped$drainage_area_km2 <- pour_snapped$drainage_area_m2 / 1e6
#   
#   # Save snapped pour points with drainage area as GPKG for QC / mapping
#   st_write(
#     pour_snapped,
#     file.path(path_start, "AOI_with_drainage.gpkg"),
#     delete_dsn = TRUE,
#     quiet = TRUE
#   )
#   
#   
#   # -------------------------------------------------------------------------
#   # 4. Delineate watersheds (raster) from snapped pour points
#   # -------------------------------------------------------------------------
#   message("4) Delineate watersheds from snapped pour points")
#   
#   wbt_watershed(
#     d8_pntr  = file.path(path_start, "d8_pntr.tif"),
#     pour_pts = file.path(path_start, "pour_points_snapped.shp"),
#     output   = file.path(path_start, "watersheds.tif")
#   )
#   
#   ws_rast <- rast(file.path(path_start, "watersheds.tif"))
#   
#   
#   # -------------------------------------------------------------------------
#   # 5. Convert watershed raster to polygons + compute local (geometric) area
#   # -------------------------------------------------------------------------
#   message("5) Polygonize watershed raster and compute local areas")
#   
#   # Raster → polygons; 'watersheds' is the ID of each basin
#   ws_poly <- as.polygons(ws_rast, values = TRUE, dissolve = TRUE, na.rm = TRUE)
#   
#   ws_poly_sf <- st_as_sf(ws_poly) %>%
#     st_make_valid() %>%
#     st_transform(st_crs(target_crs)) %>%
#     st_zm(drop = TRUE, what = "ZM")
#   
#   if (nrow(ws_poly_sf) == 0) {
#     stop("Watershed polygonization produced 0 features.")
#   }
#   
#   # Geometric area of each watershed polygon
#   ws_poly_sf <- ws_poly_sf %>%
#     mutate(
#       ws_local_area_m2  = as.numeric(st_area(geometry)),
#       ws_local_area_km2 = ws_local_area_m2 / 1e6
#     )
#   
#   
#   # -------------------------------------------------------------------------
#   # 6. Attach FAC-based upstream area to polygons using watershed ID
#   # -------------------------------------------------------------------------
#   message("6) Attach cumulative upstream area from FAC to polygons")
#   
#   # 6.1 Get watershed ID for each snapped pour point from ws_rast
#   ws_id_vals <- terra::extract(ws_rast, vect(pour_snapped))
#   # 2nd column is watershed ID
#   pour_snapped$watersheds <- ws_id_vals[, 2]
#   
#   # 6.2 Build lookup of upstream area per watershed ID from FAC values
#   da_lookup <- pour_snapped %>%
#     st_drop_geometry() %>%
#     dplyr::select(watersheds, drainage_area_m2, drainage_area_km2) %>%
#     dplyr::filter(!is.na(watersheds)) %>%
#     dplyr::distinct(watersheds, .keep_all = TRUE)
#   
#   # 6.3 Join FAC-based upstream area onto watershed polygons
#   ws_poly_sf_upstream <- ws_poly_sf %>%
#     left_join(da_lookup, by = "watersheds") %>%
#     rename(
#       ws_cum_up_m2  = drainage_area_m2,
#       ws_cum_up_km2 = drainage_area_km2
#     )
#   
#   # NOTE:
#   #   - ws_local_area_* = geometric area of each watershed polygon
#   #   - ws_cum_up_*     = D8 FAC upstream area at pour point (more accurate
#   #                       representation of contributing area)
#   
#   
#   # -------------------------------------------------------------------------
#   # 7. Write output polygons and return
#   # -------------------------------------------------------------------------
#   message("7) Writing watershed polygons to GPKG")
#   
#   st_write(
#     ws_poly_sf_upstream,
#     file.path(path_start, "watersheds_polygons_o.gpkg"),
#     delete_dsn = TRUE,
#     quiet = TRUE
#   )
#   
#   return(ws_poly_sf_upstream)
# }


library(sf)
library(terra)
library(whitebox)
library(dplyr)

# ---------------------------------------------------------------------------
# make_watershed_simple()
#
# Arguments:
#   path_start   : directory where intermediates + outputs will be written
#   dem_utm_path : path to DEM (projected, metres; e.g., UTM)
#   dem_utm      : terra::rast of same DEM (used only for reprojection/QC)
#   pour_points  : path to *POINT* layer (AOI / XS / pour points)
#
# Outputs written into path_start:
#   - dem_breached.tif
#   - d8_pntr.tif
#   - fac_area.tif
#   - pour_points_in_raster_crs.shp
#   - pour_points_snapped.shp
#   - pour_points_snapped_with_drainage.gpkg   (snapped points + DA)
#   - watersheds.tif
#   - watersheds_polygons_simple.gpkg          (watershed polygons + local area)
#
# Return value:
#   sf polygon layer of watersheds with:
#     - watershed_id       : integer ID from watershed raster
#     - ws_local_area_m2   : polygon area (m²)
#     - ws_local_area_km2  : polygon area (km²)
#
# You can then:
#   - use pour_points_snapped_with_drainage.gpkg to get drainage area per point
#   - join / QC however you like downstream.
# ---------------------------------------------------------------------------

library(sf)
library(terra)
library(whitebox)
library(dplyr)

# ---------------------------------------------------------------------------
# make_watershed_simple()
#
# Minimal helper:
#   - preprocess DEM
#   - compute FAC
#   - snap pour points
#   - get drainage area per snapped point
#   - delineate watersheds and polygonize
#
# No FID logic, no cumulative area – you handle that downstream.
# ---------------------------------------------------------------------------

make_watershed_simple <- function(path_start, dem_utm_path, dem_utm, pour_points) {
  
  # -------------------------------------------------------------------------
  # 1. DEM pre-processing: breach, D8 pointer, D8 flow accumulation (FAC)
  # -------------------------------------------------------------------------
  message("1) DEM pre-processing: breach / pointer / flow accumulation")
  
  # 1.1 Breach depressions
  wbt_breach_depressions_least_cost(
    dem    = dem_utm_path,
    output = file.path(path_start, "dem_breached.tif"),
    dist   = 30
  )
  
  # 1.2 D8 flow direction
  wbt_d8_pointer(
    dem    = file.path(path_start, "dem_breached.tif"),
    output = file.path(path_start, "d8_pntr.tif")
  )
  
  # 1.3 D8 flow accumulation as catchment area (m²)
  wbt_d8_flow_accumulation(
    input    = file.path(path_start, "d8_pntr.tif"),
    output   = file.path(path_start, "fac_area.tif"),
    out_type = "catchment area",
    pntr     = TRUE
  )
  
  r_fac      <- rast(file.path(path_start, "fac_area.tif"))   # m²
  target_crs <- crs(r_fac)
  
  
  # -------------------------------------------------------------------------
  # 2. Read pour points, standardize geometry/CRS, write temp shapefile
  # -------------------------------------------------------------------------
  message("2) Read + standardize pour points")
  
  pp_raw <- st_read(pour_points, quiet = TRUE)
  
  pp_std <- pp_raw %>%
    st_zm(drop = TRUE, what = "ZM") %>%
    st_cast("POINT") %>%
    st_transform(st_crs(target_crs))
  
  # geometry-only shapefile to keep Whitebox happy
  pp_tmp <- file.path(path_start, "pour_points_in_raster_crs.shp")
  st_write(pp_std["geometry"], pp_tmp, delete_dsn = TRUE, quiet = TRUE)
  
  
  # -------------------------------------------------------------------------
  # 3. Snap pour points to high-FAC cells and extract drainage area
  # -------------------------------------------------------------------------
  message("3) Snap pour points to FAC maxima and extract drainage area")
  
  wbt_snap_pour_points(
    pour_pts   = pp_tmp,
    flow_accum = file.path(path_start, "fac_area.tif"),
    output     = file.path(path_start, "pour_points_snapped.shp"),
    snap_dist  = 10
  )
  
  pour_snapped <- st_read(
    file.path(path_start, "pour_points_snapped.shp"),
    quiet = TRUE
  ) %>%
    st_zm(drop = TRUE, what = "ZM") %>%
    st_transform(st_crs(target_crs))
  
  # IMPORTANT: Whitebox often creates a non-integer FID field.
  # Rename or drop it so GDAL doesn't choke when writing to GPKG.
  if ("FID" %in% names(pour_snapped)) {
    pour_snapped <- pour_snapped %>%
      dplyr::rename(fid_wbt = FID)
  }
  
  # Extract FAC (catchment area m²) at each snapped point
  da_vals <- terra::extract(r_fac, vect(pour_snapped))
  fac_col_name <- setdiff(names(da_vals), "ID")
  
  pour_snapped$drainage_area_m2  <- da_vals[[fac_col_name]]
  pour_snapped$drainage_area_km2 <- pour_snapped$drainage_area_m2 / 1e6
  
  
  out_pts_gpkg <- file.path(path_start, "pour_points_snapped_with_drainage.gpkg")
  
  # If a file already exists, remove it first (silently)
  if (file.exists(out_pts_gpkg)) {
    unlink(out_pts_gpkg, recursive = TRUE, force = TRUE)
  }
  st_write(
    pour_snapped,
    file.path(path_start, "pour_points_snapped_with_drainage.gpkg"),
    delete_dsn = TRUE,
    delete_layer = TRUE,
    quiet = TRUE
  )
  
  
  # -------------------------------------------------------------------------
  # 4. Delineate watersheds (raster) from snapped pour points
  # -------------------------------------------------------------------------
  message("4) Delineate watersheds from snapped pour points")
  
  wbt_watershed(
    d8_pntr  = file.path(path_start, "d8_pntr.tif"),
    pour_pts = file.path(path_start, "pour_points_snapped.shp"),
    output   = file.path(path_start, "watersheds.tif")
  )
  
  ws_rast <- rast(file.path(path_start, "watersheds.tif"))
  names(ws_rast) <- "watershed_id"
  
  
  # -------------------------------------------------------------------------
  # 5. Polygonize watershed raster + compute local polygon areas
  # -------------------------------------------------------------------------
  message("5) Polygonize watershed raster and compute local areas")
  
  ws_poly <- as.polygons(ws_rast, values = TRUE, dissolve = TRUE, na.rm = TRUE)
  
  ws_poly_sf <- st_as_sf(ws_poly) %>%
    st_make_valid() %>%
    st_transform(st_crs(dem_utm)) %>%  # for consistent area in m²
    st_zm(drop = TRUE, what = "ZM")
  
  if (nrow(ws_poly_sf) == 0) {
    stop("Watershed polygonization produced 0 features.")
  }
  
  # Make sure watershed ID column is named watershed_id
  non_geom_cols <- setdiff(names(ws_poly_sf), attr(ws_poly_sf, "sf_column"))
  if (!"watershed_id" %in% non_geom_cols) {
    ws_poly_sf <- ws_poly_sf %>%
      rename(watershed_id = all_of(non_geom_cols[1]))
  }
  
  ws_poly_sf <- ws_poly_sf %>%
    mutate(
      ws_local_area_m2  = as.numeric(st_area(geometry)),
      ws_local_area_km2 = ws_local_area_m2 / 1e6
    )
  
  st_write(
    ws_poly_sf,
    file.path(path_start, "watersheds_polygons_simple.gpkg"),
    delete_dsn = TRUE,
    quiet = TRUE
  )
  
  return(ws_poly_sf)
}

