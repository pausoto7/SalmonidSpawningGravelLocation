

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


# ------------------------

stream_properties <- clean_bankfull %>%
  group_by(OBJ_ORDER) %>%
  summarise(h = max(RASTERVALU) - min(RASTERVALU), 
            w = max(distance) - min(distance)) 
  









