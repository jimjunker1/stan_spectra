library(tidyverse)
library(brms)
library(rstan)
library(tidybayes)
library(janitor)
library(ggridges)

# Does model fit improve when using a global xmin and xmax or with sample specific xmin and xmax?
# varying xmins and xmax ----------------------------------------
macro_fish_thin = readRDS("data/macro_fish_thin.rds") %>% ungroup
test <- readRDS("data/macro_fish_thin.rds") %>% ungroup %>%  
  filter(site_id_int <= 5) %>% 
  sample_n(1000)

# 2) make stan data
stan_data = list(N = nrow(test),
                 mat_s = test$mat_s,
                 year = as.integer(as.factor(test$year)),
                 n_years = length(unique(test$year)),
                 n_sites = length(unique(test$site_id_int)),
                 site = test$site_id_int,
                 counts = test$no_m2,
                 x = test$dw,
                 xmin = test$xmin,
                 xmax = test$xmax)


mod_varxminxmax = stan(file = "models/stan_spectra_mod.stan", 
                   data = stan_data,
                   iter = 300, chains = 1, cores = 4)




# 1) global xmin and xmax ----------------------------------------
test_global <- macro_fish_thin %>% 
  mutate(xmin = min(dw),
         xmax = max(dw))

# 2) make stan data
stan_data_global = list(N = nrow(test_global),
                 mat_s = test_global$mat_s,
                 year = as.integer(as.factor(test_global$year)),
                 n_years = length(unique(test_global$year)),
                 n_sites = length(unique(test_global$site_id_int)),
                 site = test_global$site_id_int,
                 counts = test_global$no_m2,
                 x = test_global$dw,
                 xmin = test_global$xmin,
                 xmax = test_global$xmax)


mod_globalxminxmax = stan(file = "models/stan_spectra_mod.stan", 
                       data = stan_data_global,
                       iter = 400, chains = 1, cores = 4)



mod_globalxminxmax


bs = rPLB(n)

test_global


test = test_global %>% 
  mutate(prop = no_m2/max(no_m2)*10000) %>% 
  arrange(-prop)


test_global %>% 
  ggplot(aes(x = dw*no_m2)) + 
  geom_histogram(bins = 100) + 
  # scale_x_log10() + 
  NULL

sims = tibble(dw = rPLB(10000, b = -1.26, xmin = min(test_global$dw), xmax = max(test_global$dw)))
sims %>% 
  ggplot(aes(x = dw)) + 
  geom_freqpoly(bins = 100) + 
  # scale_x_log10() +
  NULL
