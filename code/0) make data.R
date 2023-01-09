library(tidyverse)
library(janitor)

gpp = readRDS("data/gpp_means.rds") %>% clean_names() %>% 
  rename(gpp = mean, gpp_sd = sd) %>% 
  mutate(log_gpp = log(gpp),
         log_gpp_s = scale(log_gpp, center = T, scale = T),
         log_gpp_s = as.numeric(log_gpp_s),
         gpp_s = scale(gpp, center = T, scale = T)) 

macro_fish_mat = readRDS("data/macro_fish_mat.rds")

macro_fish_mat_siteminmax <- macro_fish_mat %>% # this is really site max and global min
  left_join(gpp) %>% 
  filter(!is.na(log_gpp_s)) %>% 
  group_by(mat_s, log_gpp_s, gpp_s, gpp, log_gpp, site_id, year, ID, mat_site, sdat_site, dw, x) %>% 
  summarize(counts = sum(counts)) %>% 
  ungroup %>% 
  mutate(sample_id_int = as.integer(as.factor(ID))) %>% 
  ungroup %>%  
  mutate(xmin = min(dw))%>% 
  group_by(site_id) %>%
  mutate(xmax = max(dw)) %>% 
  ungroup %>% 
  mutate(id = 1:nrow(.),
         year_id = as.integer(as.factor(year)),
         site_id_int = as.integer(as_factor(site_id)))

saveRDS(macro_fish_mat_siteminmax, file = "data/macro_fish_mat_siteminmax.rds")



# biomass -----------------------------------------------------------------

biomass = macro_fish_mat %>% # this is really site max and global min
  left_join(gpp) %>% 
  filter(!is.na(log_gpp_s)) %>% 
  filter(animal_type == "macroinvertebrates") %>% 
  mutate(dw_total = x*counts) %>% 
  group_by(mat_s, log_gpp_s, gpp_s, gpp, log_gpp, site_id, year, ID, mat_site, sdat_site) %>% 
  summarize(sample_biomass = sum(dw_total)/1000)

saveRDS(biomass, file = "data/biomass.rds")


