library(mapview)
library(tidyverse)


make_watershed <- function(river_name, dem_utm_path, dem_utm, pour_points){
  
  
  path_start <- sprintf("spatial data/%s/watershed_files/", river_name)

  
  
  
  print("breach/fill")
  # 1a. Breach/fill depressions so flow is continuous
  wbt_breach_depressions_least_cost(
    dem   = dem_utm_path,
    output = paste0(path_start, "dem_breached.tif"), 
    dist = 25                       # MAY NEED TO BE MODIFIED LATER?
  )
  
  print("pointer")
  # 1b. D8 flow direction (pointer)
  wbt_d8_pointer(
    dem    = paste0(path_start, "dem_breached.tif"),
    output = paste0(path_start, "d8_pntr.tif")
  )
  
  print("flow accumulation")
  # 1c. D8 flow accumulation, outputting *catchment area* (m²)
  wbt_d8_flow_accumulation(
    input  = paste0(path_start, "d8_pntr.tif"),
    output = paste0(path_start,"fac_area.tif"),
    out_type = "catchment area",
    pntr   = TRUE
  )
  
  print("snap pour points")
  # snap to highest flow-accum cell within, say, 100 m
  wbt_snap_pour_points(
    pour_pts   = pour_points,
    flow_accum = paste0(path_start,"fac_area.tif"),
    output     = paste0(path_start, "pour_points_snapped.shp"),
    snap_dist  = 100  # metres
  )
  
  pour_snapped <- st_read(paste0(path_start, "pour_points_snapped.shp"))
  
  
  fac <- rast(paste0(path_start,"fac_area.tif"))
  
  # terra::extract returns a data frame; second column is raster value
  da_vals <- terra::extract(fac, vect(pour_snapped))
  
  pour_snapped$drainage_area_m2  <- da_vals[, 2]
  pour_snapped$drainage_area_km2 <- pour_snapped$drainage_area_m2 / 1e6
  
  st_write(pour_snapped, paste0(path_start, "thalwag_with_drainage_area.gpkg"), delete_dsn = TRUE)
  
  
  wbt_watershed(
    d8_pntr = paste0(path_start, "d8_pntr.tif"),
    pour_pts = paste0(path_start, "pour_points_snapped.shp"),
    output = paste0(path_start,"watersheds.tif")
  )
  
  ws_rast  <- rast(paste0(path_start,"watersheds.tif"))
  ws_poly  <- as.polygons(ws_rast, dissolve = TRUE, values = TRUE)
  ws_poly_sf <- st_as_sf(ws_poly)
  
  
  
  
  # Add area columns to watershed polygons
  ws_poly_sf <- ws_poly_sf %>%
    mutate(
      ws_area_m2  = as.numeric(st_area(geometry)),   # area in m²
      ws_area_km2 = ws_area_m2 / 1e6                 # area in km²
    )
  
  
  
  st_write(ws_poly_sf,
           paste0(path_start, "watersheds_polygons.shp"),
           delete_dsn = TRUE)
  
  
  
  ws_poly_sf %>% st_geometry_type() %>% table()
  summary(ws_poly_sf$ws_area_km2)
  
  
  
  
  
  
  # Keep only the area columns (and any ID, if present)
  ws_area_only <- ws_poly_sf %>%
    select(ws_area_m2, ws_area_km2)
  
  # Spatial join: attach watershed area to each pour point
  pour_with_da <- st_join(
    pour_snapped,     # points
    ws_area_only,     # polygons with area
    join = st_within, # each point should be within its watershed polygon
    left = TRUE
  )
  
  # Inspect result
  head(st_drop_geometry(pour_with_da))
  
  st_write(pour_with_da, paste0(path_start, "thalwag_with_drainage_area.gpkg"), delete_dsn = TRUE)
  
  
  mapview(dem_utm) +
    mapview(ws_poly_sf, alpha.regions = 0.2, col.regions = NA, color = "black") +
    mapview(pour_with_da, zcol = "ws_area_km2")
  
  
  return(ws_poly_sf)
}
