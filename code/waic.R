library(loo)
library(VGAM)
library(rstan)
library(brms)
library(tidyverse)
library(tidybayes)
library(hciR)

# waic on simulated data

# the model already calculates the log likelihood using generated quantities
# We'll compare stan's automation to doing it by hand after the fact with just R
fit_model = stan_model("models/old/b_paretocounts_singlesample.stan")

x = 0.002
lower = 0.001
upper = 1000
shape = 0.5

N = 1000
x = rtruncpareto(N, lower, upper, shape)
counts = rep(1, N)
lower = rep(lower, N)
upper = rep(upper, N)

stan_dat <- list(x = x,
                   N = N,
                   counts = counts,
                   xmax = upper,
                   xmin = lower)

fit <- sampling(object = fit_model,
                  data = stan_dat,
                  iter = 1000,
                  chains = 2,
                  open_progress = F,
                  verbose = F)
fit

loglik = extract_log_lik(fit)

loo::waic(loglik)


# calculate by hand
posts =as_draws_df(fit) 

# log likelihood
loglik_fit = posts %>% 
  select(lambda, .draw) %>% 
  expand_grid(stan_dat %>% as_tibble() %>% mutate(id = 1:nrow(.))) %>% 
  mutate(loglik = counts*(log((lambda+1) / ( xmax^(lambda+1) - xmin^(lambda+1))) + lambda*log(x)))

# convert to matrix
log_lik_matrix = loglik_fit %>% ungroup %>% 
  # filter(.draw < 3) %>% 
  select(loglik, .draw, id) %>% 
  # filter(.draw <= 100) %>% 
  pivot_wider(names_from = id, values_from = loglik) %>%
  select(-.draw) %>% 
  as.matrix() %>% unname()

# waic with loo
loo::waic(log_lik_matrix)


# use real data --------------------------------------------
fit_temponly = readRDS("models/fit_temponly.rds")
fit_interaction = readRDS("models/fit_interaction.rds")

dat = readRDS(file = "data/macro_fish_mat_siteminmax.rds") 


posts = as_draws_df(fit_temponly) %>% 
  pivot_longer(cols = contains('alpha_raw_site'),
               names_to = "site_group", values_to = "site_offset") %>% 
  pivot_longer(cols = contains("alpha_raw_year"),
               names_to = "year_group", values_to = "year_offset") %>% 
  mutate(site_id_int = parse_number(site_group),
         year_id = parse_number(year_group)) %>%
  filter(.draw <= 500) %>% 
  right_join(dat %>% ungroup %>% distinct(mat_s, site_id_int, year_id)) %>% 
  mutate(lambda = a + beta_mat*mat_s + sigma_year*year_offset + sigma_site*site_offset)

loglik_fit = posts %>% 
  filter(.draw <= 200) %>% 
  right_join(dat %>% distinct(year_id, site_id_int, id, x, xmax, xmin, counts, .keep_all = T)) %>% 
  mutate(loglik = counts*(log((lambda+1) / ( xmax^(lambda+1) - xmin^(lambda+1))) + lambda*log(x)))

log_lik_matrix_temp = loglik_fit %>% ungroup %>% 
  filter(.draw <= 200) %>%
  select(loglik, .draw, id) %>% 
  pivot_wider(names_from = id, values_from = loglik) %>%
  select(-.draw) %>% 
  as.matrix() %>% unname()

# interaction
posts = as_draws_df(fit_temponly) %>% 
  pivot_longer(cols = contains('alpha_raw_site'),
               names_to = "site_group", values_to = "site_offset") %>% 
  pivot_longer(cols = contains("alpha_raw_year"),
               names_to = "year_group", values_to = "year_offset") %>% 
  mutate(site_id_int = parse_number(site_group),
         year_id = parse_number(year_group)) %>%
  filter(.draw <= 500) %>% 
  right_join(dat %>% ungroup %>% distinct(mat_s, log_gpp_s, site_id_int, year_id)) %>% 
  mutate(lambda = a + beta_mat*mat_s + 
           # beta_gpp*log_gpp_s + 
           # beta_gpp_mat*log_gpp_s*mat_s + 
           sigma_year*year_offset + sigma_site*site_offset)

loglik_fit = posts %>% 
  # filter(.draw <= 200) %>% 
  right_join(dat %>% distinct(year_id, site_id_int, id, x, xmax, xmin, counts, .keep_all = T)) %>% 
  mutate(loglik = counts*(log((lambda+1) / ( xmax^(lambda+1) - xmin^(lambda+1))) + lambda*log(x)))

log_lik_matrix_interaction = loglik_fit %>% ungroup %>% 
  # filter(.draw <= 800) %>%
  select(loglik, .draw, id) %>% 
  pivot_wider(names_from = id, values_from = loglik) %>%
  select(-.draw) %>% 
  as.matrix() %>% unname()


waic_temp = loo::waic(log_lik_matrix_temp)
waic_interaction = loo::waic(log_lik_matrix_interaction)


waic_temp$estimates[[3]] - waic_interaction$estimates[[3]]

waic_tbl = tibble(waic = c(waic_temp$estimates[[3]], waic_interaction$estimates[[3]]),
       se = c(waic_temp$estimates[[6]], waic_interaction$estimates[[6]]),
       model = c("temp", "interaction"))

waic_tbl %>% 
  ggplot(aes(x = model, y = waic, ymin = waic - se, ymax = waic + se)) + 
  geom_pointrange()

# 
