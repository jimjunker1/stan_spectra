library(rstan)
library(tidyverse)
library(janitor)
library(tidybayes)
library(brms)

interaction_posts_wrangled = readRDS(file = "posteriors/interaction_posts_wrangled.rds")
temp_posts_wrangled = readRDS(file = "posteriors/temp_posts_wrangled.rds")
dat = readRDS(file = "data/macro_fish_mat_siteminmax.rds") 

sim_postpreds = function(model, sample_n = 10000){
  interaction_posts_wrangled %>% 
    filter(.draw <= 10) %>% 
    right_join(dat) %>% 
    mutate(u = runif(nrow(.), 0, 1)) %>% # uniform draw
    mutate(x = (u*xmax^(lambda+1) +  (1-u) * xmin^(lambda+1) ) ^ (1/(lambda+1))) %>% 
    mutate(data = "y_rep") %>% 
    bind_rows(dat %>% 
                mutate(data = "y_raw") %>% 
                mutate(.draw = 0)) %>% 
    rename(sim = x) %>% 
    group_by(.draw) %>% 
    sample_n(sample_n, weight = counts, replace = T) %>% 
    select(sim, .draw, data, site_id, year, sample_id_int)
  
}

sim_postpreds(interaction_posts_wrangled) %>%
  filter(sample_id_int == 5) %>% 
  ggplot(aes(x = .draw, y = sim, color = data, group = .draw)) + 
  geom_violin() +
  scale_y_log10() +
  NULL

sim_postpreds(temp_posts_wrangled) %>% 
  filter(sample_id_int <= 40) %>% 
  ggplot(aes(x = sim, color = data, group = .draw)) +
  geom_density() +
  scale_x_log10() +
  facet_wrap(~sample_id_int, scales = "free_y") +
  NULL


