library(tidyverse)
library(janitor)

macro_fish_dw <- readRDS("C:/Users/Jeff.Wesner/OneDrive - The University of South Dakota/USD/Github Projects/neon_size_spectra/data/derived_data/macro_fish_dw.rds")

macro_fish_dw %>% 
  ungroup %>% 
  mutate(mat_s = (mat_site - mean(mat_site)/sd(mat_site))) %>% 
  ungroup() %>% 
  filter(mat_s == min(mat_s)) %>% 
  select(mat_s)

mle_mat <- read_csv("C:/Users/Jeff.Wesner/OneDrive - The University of South Dakota/USD/Github Projects/neon_size_spectra/data/derived_data/mle_mat.csv")

macro_fish_mat = macro_fish_dw %>% left_join(mle_mat) %>% 
  group_by(ID) %>% 
  mutate(xmin = min(dw),
         xmax = max(dw),
         x = dw,
         counts = no_m2) %>% 
  ungroup() %>%
  mutate(site_no = as.numeric(factor(site_id)),
         mat_s = (mat_site - mean(mat_site))/sd(mat_site)) %>% 
  mutate(site_id_int = as.integer(as_factor(site_id))) %>%
  mutate(IDname = paste0(site_id, year_month),
         ID = as.integer(as.factor(IDname))) 

macro_fish_mat_globalminmax <- readRDS("data/macro_fish_mat.rds") %>% 
  mutate(site_id_int = as.integer(as_factor(site_id))) %>% 
  ungroup() %>% 
  mutate(xmin = min(dw),
         xmax = max(dw)) %>% 
  ungroup


macro_fish_mat_siteminmax <- readRDS("data/macro_fish_mat.rds") %>% 
  mutate(site_id_int = as.integer(as_factor(site_id))) %>% 
  ungroup() %>%  
  mutate(xmin = min(dw))%>% 
  group_by(site_id) %>%
  mutate(xmax = max(dw)) %>% 
  ungroup

saveRDS(macro_fish_mat_siteminmax, file = "data/macro_fish_mat_siteminmax.rds")
saveRDS(macro_fish_mat_globalminmax , file = "data/macro_fish_mat_globalminmax.rds")
saveRDS(macro_fish_mat, file = "data/macro_fish_mat.rds")




