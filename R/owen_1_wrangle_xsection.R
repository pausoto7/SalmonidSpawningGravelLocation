

library(readxl)
library(tidyverse)

xsection_raw <- read_xls("tabular data/owen_crk/points_with_elev_BWK_v2.xls") 

xsection <- xsection_raw %>%
  arrange(OBJECTID) %>%
  group_by(ORIG_FID) %>%
  mutate(distance = 0.5*row_number()-0.5) %>%
  ungroup() %>%
  mutate(OBJ_ORDER = recode(ORIG_FID, 
                            `6` = 1, 
                            `5` = 2, 
                            `4`= 3,
                            `3`= 4, 
                            `2`= 5, 
                            `1`= 6, 
                            .default = ORIG_FID))


ggplot(xsection) + 
  geom_point(aes(x= distance, y = RASTERVALU,  group = OBJ_ORDER, color = as.factor(OBJ_ORDER))) + 
  geom_line(aes(x= distance, y = RASTERVALU,  group = OBJ_ORDER, color = as.factor(OBJ_ORDER))) + 
  theme_bw() +
  facet_wrap(.~OBJ_ORDER)



source("R/trim_to_bankfull.R")


clean_bankfull <- trim_to_bankfull(xsection, edge_prop = 0.4)

ggplot(clean_bankfull ) + 
  geom_point(aes(x= distance, y = RASTERVALU,  group = OBJ_ORDER, color = as.factor(OBJ_ORDER))) + 
  geom_line(aes(x= distance, y = RASTERVALU,  group = OBJ_ORDER, color = as.factor(OBJ_ORDER))) + 
  theme_bw() +
  facet_wrap(.~OBJ_ORDER)



thalwag_data <- clean_bankfull %>%
  group_by(OBJ_ORDER) %>%
  mutate(deepest_point = min(RASTERVALU)) %>%
  distinct(OBJ_ORDER, deepest_point, .keep_all = TRUE)

save(thalwag_data, file = "tabular data/owen_crk/thalwag_data.RData")


# # 1a. Breach/fill depressions so flow is continuous
# wbt_breach_depressions_least_cost(
#   dem   = "dem_utm.tif",
#   output = "dem_breached.tif", 
#   dist = 25                       # MAY NEED TO BE MODIFIED LATER?
# )
# 
# # 1b. D8 flow direction (pointer)
# wbt_d8_pointer(
#   dem    = "dem_breached.tif",
#   output = "d8_pntr.tif"
# )
# 
# # 1c. D8 flow accumulation, outputting *catchment area* (m²)
# wbt_d8_flow_accumulation(
#   input  = "d8_pntr.tif",
#   output = "fac_area.tif",
#   out_type = "catchment area",
#   pntr   = TRUE
# )
# 
# # snap to highest flow-accum cell within, say, 100 m
# wbt_snap_pour_points(
#   pour_pts   = "pour_points.shp",
#   flow_accum = "fac_area.tif",
#   output     = "pour_points_snapped.shp",
#   snap_dist  = 100  # metres
# )
# 
# pour_snapped <- st_read("pour_points_snapped.shp")
# 
# 
# fac <- rast("fac_area.tif")
# 
# # terra::extract returns a data frame; second column is raster value
# da_vals <- terra::extract(fac, vect(pour_snapped))
# 
# pour_snapped$drainage_area_m2  <- da_vals[, 2]
# pour_snapped$drainage_area_km2 <- pour_snapped$drainage_area_m2 / 1e6
# 
# st_write(pour_snapped, "thalwag_with_drainage_area.gpkg", delete_dsn = TRUE)
# 
# 
# wbt_watershed(
#   d8_pntr = "d8_pntr.tif",
#   pour_pts = "pour_points_snapped.shp",
#   output = "watersheds.tif"
# )
# 
# ws_rast  <- rast("watersheds.tif")
# ws_poly  <- as.polygons(ws_rast, dissolve = TRUE, values = TRUE)
# ws_poly_sf <- st_as_sf(ws_poly)
# 
# 
# 
# 
# # Add area columns to watershed polygons
# ws_poly_sf <- ws_poly_sf %>%
#   mutate(
#     ws_area_m2  = as.numeric(st_area(geometry)),   # area in m²
#     ws_area_km2 = ws_area_m2 / 1e6                 # area in km²
#   )
# 
# 
# ws_poly_sf %>% st_geometry_type() %>% table()
# summary(ws_poly_sf$ws_area_km2)
# 
# 
# 
# # Keep only the area columns (and any ID, if present)
# ws_area_only <- ws_poly_sf %>%
#   select(ws_area_m2, ws_area_km2)
# 
# # Spatial join: attach watershed area to each pour point
# pour_with_da <- st_join(
#   pour_snapped,     # points
#   ws_area_only,     # polygons with area
#   join = st_within, # each point should be within its watershed polygon
#   left = TRUE
# )
# 
# # Inspect result
# head(st_drop_geometry(pour_with_da))
# 
# st_write(pour_with_da, "thalwag_with_drainage_area.gpkg", delete_dsn = TRUE)
# 
# 
# library(mapview)
# 
# 
# mapview(dem_utm) +
#   mapview(ws_poly_sf, alpha.regions = 0.2, col.regions = NA, color = "black") +
#   mapview(pour_with_da, zcol = "ws_area_km2")
