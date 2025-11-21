

library(readxl)
library(tidyverse)

xsection_raw <- read_xls("tabular data/mamquam_river/mamquam_elevations3.xls") 

xsection <- xsection_raw %>%
  arrange(OBJECTID) %>%
  group_by(ORIG_FID) %>%
  mutate(distance = 0.5*row_number()-0.5) %>%
  ungroup() %>%
  mutate(OBJ_ORDER = ORIG_FID +1)


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



thalwag_data_m <- clean_bankfull %>%
  group_by(OBJ_ORDER) %>%
  mutate(deepest_point = min(RASTERVALU)) %>%
  distinct(OBJ_ORDER, deepest_point, .keep_all = TRUE)

save(thalwag_data_m, file = "tabular data/mamquam_river/thalwag_data.RData")



# ------------------------

stream_properties_mamq <- clean_bankfull %>%
  group_by(OBJ_ORDER) %>%
  summarise(h = max(RASTERVALU) - min(RASTERVALU), 
            w = max(distance) - min(distance)) 


save(stream_properties_mamq, file = "tabular data/mamquam_river/owen_stream_properties.RData")

print(stream_properties_mamq)


