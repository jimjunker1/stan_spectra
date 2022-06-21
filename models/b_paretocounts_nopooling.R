library(tidyverse)
library(brms)
library(rstan)
library(tidybayes)
library(sizeSpectra)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Fit to real data ----------------------------------------
macro_fish_dw <- readRDS("C:/Users/Jeff.Wesner/OneDrive - The University of South Dakota/USD/Github Projects/neon_size_spectra/data/derived_data/macro_fish_dw.rds")
mle_mat <- mle_mat <- read_csv("C:/Users/Jeff.Wesner/OneDrive - The University of South Dakota/USD/Github Projects/neon_size_spectra/data/derived_data/mle_mat.csv") %>% 
  select(-ID)
macro_fish_mat = macro_fish_dw %>% left_join(mle_mat) %>% 
  group_by(ID) %>% 
  mutate(xmin = min(dw),
         xmax = max(dw),
         x = dw,
         counts = no_m2) %>% 
  ungroup() %>% 
  mutate(site_no = as.numeric(factor(site_id)))

macro_fish_test <- macro_fish_mat

dat <- list(x = macro_fish_test$x,
            counts = macro_fish_test$counts,
            N = nrow(macro_fish_test),
            n_groups = length(unique(macro_fish_test$site_no)),
            site = macro_fish_test$site_no,
            # sum_counts = fake$sum_counts,
            xmax = macro_fish_test$xmax,
            xmin = macro_fish_test$xmin)

rt <- stanc(file="models/b_paretocounts_nopooling.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)      # compile

fit_pareto_neon_nopool <- sampling(sm, data=dat,iter = 1000, chains = 1)

# saveRDS(fit_pareto_neon_nopool, file = "models/fit_pareto_neon_nopool.rds")
fit_pareto_neon_nopool <- readRDS(file = "models/fit_pareto_neon_nopool.rds")


# compare to sizeSpectra
compare_neon_eight <- fit_pareto_neon_varint %>% as_draws_df() %>%
  select(starts_with("b_exp")) %>% 
  ungroup() %>% 
  pivot_longer(starts_with("b_exp")) %>%  
  mutate(site_no = parse_number(name)) %>% 
  group_by(site_no) %>% 
  rename(b = value) %>% 
  median_qi(b) %>% 
  mutate(model = "Stan") %>% 
  add_row(b = unique(macro_fish_test$b),
          .lower = unique(macro_fish_test$confMin),
          .upper = unique(macro_fish_test$confMax),
          model = "eightMethods.count",
          site_no = macro_fish_test %>% distinct(site_id, ID, site_no) %>% pull(site_no)) %>% 
  left_join(macro_fish_test %>% ungroup() %>% distinct(site_no, site_id))


compare_neon_eight %>%
  ggplot(aes(x = site_id, y = b, ymin = .lower, ymax = .upper)) + 
  geom_pointrange(aes(color = model), size = 0.1,position = position_dodge(width = 0.1))



# Check chains ------------------------------------------------------------

traceplot(fit_pareto_neon_nopool)







