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
  mutate(site_no = as.numeric(factor(site_id)),
         mat_s = (mat_site - mean(mat_site))/sd(mat_site))

macro_fish_test <- macro_fish_mat %>% filter(ID <= 30)

dat <- list(x = macro_fish_test$x,
            counts = macro_fish_test$counts,
            N = nrow(macro_fish_test),
            n_groups = length(unique(macro_fish_test$year)),
            year = as.integer(as.factor(macro_fish_test$year)),
            # sum_counts = fake$sum_counts,
            xmax = macro_fish_test$xmax,
            xmin = macro_fish_test$xmin,
            mat_s = macro_fish_test$mat_s,
            site_no = as.integer(as.factor(macro_fish_test$site_id)),
            n_sites = length(unique(macro_fish_test$site_id)))

rt <- stanc(file="models/b_paretocounts_varyingintercepts_regression.stan")
sm <- stan_model(stanc_ret = rt, verbose=FALSE)      # compile

fit_pareto_neon_varint_regression  <- sampling(sm, data=dat,iter = 1000, chains = 4)


rt_site <- stanc(file="models/b_paretocounts_varyingintercepts_regression_yearsite.stan")
sm_site <- stan_model(stanc_ret = rt_site, verbose=FALSE)      # compile

fit_pareto_neon_varint_regression_siteyear  <- sampling(sm_site, data=dat,iter = 1000, chains = 4)


# saveRDS(fit_pareto_neon_varint_regression_siteyear, file = "models/fit_pareto_neon_varint_regression_siteyear.rds")
fit_pareto_neon_varint_regression_siteyear  <- readRDS(file = "models/fit_pareto_neon_varint_regression_siteyear.rds")

# Compare to sizeSpectra --------------------------------------------------
posts_reg <- fit_pareto_neon_varint_regression %>% as_draws_df() %>% 
  select(!starts_with("b_exp")) %>% 
  expand_grid(mat_s = unique(dat$mat_s)) %>% 
  ungroup() %>% 
  mutate(b_exp = bexp_overall + beta*mat_s)

posts_reg %>% 
  group_by(mat_s) %>% 
  median_qi(b_exp) %>% 
  ggplot(aes(x = mat_s, y = b_exp)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2) +
  geom_point(data = macro_fish_test, aes(y = b, x = mat_s)) + 
  ylim(-3, 0)



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
  geom_pointrange(aes(color = model), position = position_dodge(width = 0.1))


confint(lm(b ~ mat_site, data = macro_fish_test))
fit_pareto_neon_varint_regression %>% as_draws_df() %>% 
  median_qi(Intercept, beta_1)
