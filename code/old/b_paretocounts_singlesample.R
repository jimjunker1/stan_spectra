library(tidyverse)
library(brms)
library(rstan)
library(tidybayes)
library(sizeSpectra)

# Fit to simulated data ----------------------------------------
N_ind = 2000
xmin = 0.01
xmax = 1000
b = -0.5
alpha = -b-1
x_ind = rPLB(n = N_ind, b = b, xmin = xmin, xmax = xmax) # simulate individual masses from bounded power law


sim_data <- tibble(x = round(x_ind,2)) %>% 
  count(x, name = "counts") %>% 
  mutate(xmin = xmin,
         xmax = xmax,
         n = sum(counts))

fit_b_pareto_sim <- stan(file = "models/b_paretocounts_singlesample.stan",
                        data = list(x = sim_data$x,
                                    N = nrow(sim_data),
                                    counts = sim_data$counts,
                                    xmax = sim_data$xmax,
                                    xmin = sim_data$xmin),
                        iter = 100,
                        chains = 1) 

fit_b_pareto_sim

fit_b_pareto_sim %>% as_draws_df()

# Fit to real data ----------------------------------------
macro_fish_dw <- readRDS("C:/Users/Jeff.Wesner/OneDrive - The University of South Dakota/USD/Github Projects/neon_size_spectra/data/derived_data/macro_fish_dw.rds")
mle_mat <- mle_mat <- read_csv("C:/Users/Jeff.Wesner/OneDrive - The University of South Dakota/USD/Github Projects/neon_size_spectra/data/derived_data/mle_mat.csv") %>% 
  select(-ID)

macro_fish_mat = macro_fish_dw %>% left_join(mle_mat) %>% 
  group_by(ID) %>% 
  mutate(xmin = min(dw),
         xmax = max(dw),
         x = dw,
         counts = no_m2)

macro_fish_test <- macro_fish_mat %>% filter(site_id == "OKSR")

fit_pareto_neon <- stan(file = "models/b_paretocounts_singlesample.stan",
                        data = list(x = macro_fish_test$x,
                                    counts = macro_fish_test$counts,
                                    N = nrow(macro_fish_test),
                                    # sum_counts = fake$sum_counts,
                                    xmax = macro_fish_test$xmax,
                                    xmin = macro_fish_test$xmin),
                        iter = 100,
                        chains = 1) 

fit_pareto_neon

# compare to sizeSpectra
fit_pareto_neon %>% as_draws_df() %>% 
  ungroup() %>% 
  median_qi(b_exp) %>% 
  mutate(model = "Stan") %>% 
  add_row(b_exp = unique(macro_fish_test$b),
          .lower = unique(macro_fish_test$confMin),
          .upper = unique(macro_fish_test$confMax),
          model = "eightMethods.count")

