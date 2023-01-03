library(tidyverse)
library(janitor)

gpp = readRDS("data/gpp_means.rds") %>% clean_names() %>% 
  rename(gpp = mean, gpp_sd = sd) %>% 
  mutate(log_gpp = log(gpp),
         gpp_s = scale(gpp, center = T, scale = T),
         gpp_s = as.numeric(gpp_s)) 

macro_fish_mat_siteminmax <- readRDS("data/macro_fish_mat.rds") %>% # this is really site max and global min
  left_join(gpp) %>% 
  filter(!is.na(gpp_s)) %>% 
  group_by(mat_s, gpp_s, gpp, site_id, year, ID, mat_site, sdat_site, dw, x) %>% 
  summarize(counts = sum(counts)) %>% 
  ungroup %>% 
  mutate(site_id_int = as.integer(as.factor(ID))) %>% 
  ungroup %>%  
  mutate(xmin = min(dw))%>% 
  group_by(site_id) %>%
  mutate(xmax = max(dw)) %>% 
  ungroup

saveRDS(macro_fish_mat_siteminmax, file = "data/macro_fish_mat_siteminmax.rds")

length(unique(macro_fish_mat_siteminmax$ID))


