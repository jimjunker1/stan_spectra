library(tidyverse)
library(brms)
library(rstan)
library(tidybayes)
# library(sizeSpectra)
library(janitor)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# 1) load data ----------------------------------------
macro_fish_mat <- readRDS("data/macro_fish_mat.rds")

# 2) fit model ------------------------------------------------------------
# prior predictive
N = 1000
prior_sim <- tibble(beta = rnorm(N, 0, 1),
                    sigma_year = abs(rnorm(N, 0, 0.1)),
                    mu_year = rnorm(N, -1.5, 0.5),
                    .draw = 1:N) %>% 
  mutate(b_exp_overall = rnorm(nrow(.), mu_year, sigma_year)) %>% 
  expand_grid(mat_s = seq(-2, 2, length.out = 20)) %>% 
  mutate(b_exp_sim = b_exp_overall + beta*mat_s)

prior_sim %>% 
  ggplot(aes(x = mat_s, y = b_exp_sim, group = .draw)) + 
  geom_line(alpha = 0.2)

#load fitted model
fit_pareto_neon_varint_regression  <- readRDS(file = "models/fits/fit_pareto_neon_varint_regression.rds")

# fit full model (takes ~ 2-4 hours)
# fit_pareto_neon_varint_regression  <- stan(file = "models/b_paretocounts_varyingintercepts_regression.stan", 
#                                            data=list(x = macro_fish_mat$x,
#                                                      counts = macro_fish_mat$counts,
#                                                      N = nrow(macro_fish_mat),
#                                                      n_groups = length(unique(macro_fish_mat$year)),
#                                                      year = as.integer(as.factor(macro_fish_mat$year)),
#                                                      xmax = macro_fish_mat$xmax,
#                                                      xmin = macro_fish_mat$xmin,
#                                                      mat_s = macro_fish_mat$mat_s,
#                                                      site_no = as.integer(as.factor(macro_fish_mat$site_id)),
#                                                      n_sites = length(unique(macro_fish_mat$site_id))),
#                                            iter = 1000, chains = 4)

# saveRDS(fit_pareto_neon_varint_regression, file = "models/fits/fit_pareto_neon_varint_regression.rds")


# 3) extract posts --------------------------------------------------
posts_reg <- fit_pareto_neon_varint_regression %>% as_draws_df() %>% 
  select(!starts_with("b_exp")) %>% 
  expand_grid(mat_s = unique(macro_fish_mat$mat_s)) %>% 
  ungroup() %>% 
  mutate(b_exp = bexp_overall + beta*mat_s)


# 4) Plot -----------------------------------------------------------------

posts_reg %>% 
  group_by(mat_s) %>% 
  median_qi(b_exp) %>% 
  ggplot(aes(x = mat_s, y = b_exp)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2) +
  geom_point(data = macro_fish_mat, aes(y = b, x = mat_s)) + 
  ylim(-3, 0)

posts_reg_predgroup <- fit_pareto_neon_varint_regression %>% as_draws_df() %>% 
  clean_names() %>% 
  expand_grid(mat_s = unique(macro_fish_mat$mat_s)) %>% 
  pivot_longer(cols = starts_with("b_exp")) %>% 
  filter(name != "b_exp_overall") %>% 
  mutate(year = parse_number(name)) %>% 
  mutate(b_exp = value + beta*mat_s)

posts_reg_predgroup_summary <- posts_reg_predgroup %>% 
  group_by(mat_s, year) %>% 
  median_qi(b_exp) 
  
posts_reg_predgroup_summary %>% 
  ggplot(aes(x = mat_s, y = b_exp)) + 
  geom_line(aes(group = year)) + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper, group = year), alpha = 0.2) +
  geom_point(data = macro_fish_mat, aes(y = b, x = mat_s)) + 
  ylim(-3, 0)
