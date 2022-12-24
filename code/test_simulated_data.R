library(tidyverse)
library(rstan)
library(tidybayes)
library(brms)
library(ggridges)

# simulate data -----------------------------------------------------------

set.seed(2345)
n_sim = 200
xmax = 10000
xmin = 0.001
b = -2.05
u = runif(n_sim, min = 0, max = 1)

temp = c(5, 10, 15, 20, 25, 30)
beta_temp = -0.2
b = -1.5 + beta_temp*temp
site = as.integer(1:5)
true_temp = seq(5, 30, length.out = 5)
years = as.integer(1:3)
sigma_site = rexp(1, 9)
sigma_year = rexp(1, 9)
a = -1.5

site_rand = tibble(sites, sigma_site) %>% 
  mutate(site_rand = rnorm(nrow(.), 0, sigma_site),
         true_temp = true_temp)

year_rand = tibble(years, sigma_year) %>% 
  mutate(year_rand = rnorm(nrow(.), 0, sigma_year))

individual_sims = site_rand %>% 
  expand_grid(year_rand) %>% 
  mutate(xmax = xmax,
         xmin = xmin,
         intercept = a,
         mat_s = scale(true_temp),
         beta_temp = beta_temp,
         b = intercept + site_rand + year_rand + beta_temp*mat_s) %>% 
  expand_grid(1:n_sim) %>%
  mutate(u = runif(nrow(.), min = 0, max = 1)) %>% 
  mutate(x = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1)),
         x = as.numeric(x),
         mat_s = as.numeric(mat_s), 
         no_m2 = 1)

count_sims = individual_sims %>% mutate(x = round(x,3)) %>% 
  group_by(x, sites, years, mat_s, xmin, xmax) %>% 
  count(x, name = "counts")


# compile models ----------------------------------------------------------
stan_spectra_mod_counts = stan_model("models/stan_spectra_mod.stan")
stan_spectra_mod_nocounts = stan_model("models/stan_spectra_mod_nocounts.stan")



# fit models --------------------------------------------------------------
# individuals
stan_data_individual = list(N = nrow(individual_sims),
                            mat_s = individual_sims$mat_s,
                            year = as.integer(as.factor(individual_sims$years)),
                            n_years = length(unique(individual_sims$years)),
                            n_sites = length(unique(individual_sims$sites)),
                            site = individual_sims$sites,
                            # counts = individual_sims$no_m2,
                            x = individual_sims$x,
                            xmin = individual_sims$xmin,
                            xmax = individual_sims$xmax)

fit_ind = sampling(object = stan_spectra_mod_nocounts, 
                   data = stan_data_individual,
                   iter = 1000, chains = 1, cores = 4)

saveRDS(fit_ind, file = "models/fit_ind.rds")

# counts

stan_data_counts = list(N = nrow(count_sims),
                            mat_s = count_sims$mat_s,
                            year = as.integer(as.factor(count_sims$years)),
                            n_years = length(unique(count_sims$years)),
                            n_sites = length(unique(count_sims$sites)),
                            site = count_sims$sites,
                            counts = count_sims$counts,
                            x = count_sims$x,
                            xmin = count_sims$xmin,
                            xmax = count_sims$xmax)

fit_counts = sampling(object = stan_spectra_mod_counts, 
                   data = stan_data_counts,
                   iter = 1000, chains = 1, cores = 4)

saveRDS(fit_counts, file = "models/fit_counts.rds")


# extract posteriors ------------------------------------------------------
post_ind = as_draws_df(fit_ind) %>% as_tibble()
post_counts = as_draws_df(fit_counts) %>% as_tibble()


# compare to known parameters ---------------------------------------------

known_params = tibble(beta = beta_temp,
                      a = a,
                      sigma_year = sigma_year,
                      sigma_site = sigma_site) %>% 
  pivot_longer(everything(), values_to = "true_values") %>% 
  expand_grid(likelihood = c("paretocounts", "truncated_pareto"))

# individuals
posts_and_true_ind = post_ind %>% 
  select(beta, a, contains("sigma"), .draw) %>% 
  pivot_longer(cols = -.draw) %>% 
  mutate(likelihood = "truncated_pareto")

# counts
posts_and_true_counts = post_counts %>% 
  select(beta, a, contains("sigma"), .draw) %>% 
  pivot_longer(cols = -.draw) %>% 
  mutate(likelihood = "paretocounts")

# combine and plot

ind_and_counts = bind_rows(posts_and_true_ind, 
          posts_and_true_counts) 


recover_parameters_plot = ind_and_counts %>% 
  ggplot(aes(x = name, y = value, fill = likelihood)) +
  geom_violin() + 
  # facet_wrap(~name, scales = "free") + 
  geom_point(data = known_params, aes(y = true_values), size = 1, 
             position = position_dodge(width = 0.9)) + 
  theme_default() + 
  scale_fill_grey(start  = 0.5, end = 1) +
  labs(y = "Values",
       x = "Parameter",
       fill = "Likelihood")

library(ggview)
ggview(recover_parameters_plot, units = "in", width = 7, height = 5)

ggsave(recover_parameters_plot, file = "plots/recover_parameters_plot.jpg", width = 7, height = 5,
       dpi = 500)



