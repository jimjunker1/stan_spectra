library(tidyverse)
library(brms)

macro_fish_thin <- readRDS("data/macro_fish_thin.rds") 

# values to simulate from
pred_vals = macro_fish_thin %>% ungroup %>% distinct(mat_s, year, site_id)

# varying intercepts (non-centered)
year_offset = macro_fish_thin %>% ungroup %>% distinct(year) %>% 
  group_by(year) %>% 
  mutate(sigma_year = rexp(1, 9),
         alpha_year = rnorm(1, 0, 1),
         offset_year = sigma_year*alpha_year)

site_offset = macro_fish_thin %>% ungroup %>% distinct(site_id) %>% 
  group_by(site_id) %>% 
  mutate(sigma_site = rexp(1, 9),
         alpha_site = rnorm(1, 0, 1),
         offset_site = sigma_site*alpha_site)

# simulate from prior distributions
n_sims = 300

prior_sims = tibble(a = rnorm(n_sims, -1.5, 0.2),
       beta = rnorm(n_sims, 0, 0.1),
       alpha_raw_year = rnorm(n_sims, 0, 1),
       alpha_raw_site = rnorm(n_sims, 0, 1),
       sims = 1:n_sims)


# combine and simulate b values
b_prior_preds = pred_vals %>% 
  expand_grid(prior_sims) %>% 
  left_join(year_offset) %>% 
  left_join(site_offset) %>% 
  mutate(b_only = a + beta*mat_s)

prior_pred_b_vs_mats = b_prior_preds %>% 
  ggplot(aes(x = mat_s, y = b_only, group = sims)) + 
  geom_line(alpha = 0.2) +
  theme_default()
  
saveRDS(prior_pred_b_vs_mats, file = "plots/prior_pred_b_vs_mats.rds")      
