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
beta_temp = -0.12
beta_gpp = 0.1
beta_mat_gpp = 0.8
b = -1.5 + beta_temp*temp
site = as.integer(1:5)
true_temp = seq(5, 30, length.out = 5)
true_gpp = rnbinom(5, mu = 4000, size = 1)
years = as.integer(1:3)
sigma_site = rexp(1, 9)
sigma_year = rexp(1, 9)
a = -1.5
sites = letters[1:5]

site_rand = tibble(sites, sigma_site) %>% 
  mutate(site_rand = rnorm(nrow(.), 0, sigma_site),
         true_temp = true_temp,
         true_gpp = true_gpp)

year_rand = tibble(years, sigma_year) %>% 
  mutate(year_rand = rnorm(nrow(.), 0, sigma_year))

individual_sims = site_rand %>% 
  expand_grid(year_rand) %>% 
  mutate(xmax = xmax,
         xmin = xmin,
         intercept = a,
         mat_s = scale(true_temp),
         gpp_s = scale(true_gpp),
         beta_temp = beta_temp,
         beta_gpp = beta_gpp,
         beta_mat_gpp = beta_mat_gpp,
         b = intercept + site_rand + year_rand + 
           beta_temp*mat_s + beta_gpp*gpp_s + beta_mat_gpp*mat_s*gpp_s) %>% 
  expand_grid(1:n_sim) %>%
  mutate(u = runif(nrow(.), min = 0, max = 1)) %>% 
  mutate(x = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1)),
         x = as.numeric(x),
         mat_s = as.numeric(mat_s), 
         gpp_s = as.numeric(gpp_s),
         no_m2 = 1,
         sites = as.integer(as.factor(sites)))

count_sims = individual_sims %>% mutate(x = round(x,3)) %>% 
  group_by(x, sites, years, mat_s, xmin, xmax, gpp_s) %>% 
  count(x, name = "counts")


# compile models ----------------------------------------------------------
stan_spectra_mod_counts = stan_model("models/stan_spectra_mod.stan")
stan_spectra_mod_nocounts = stan_model("models/stan_spectra_mod_nocounts.stan")
stan_spectra_mod_gpp_x_temp = stan_model("models/stan_spectra_mod_gpp_x_temp.stan")



# fit models --------------------------------------------------------------

stan_data_interaction = list(N = nrow(count_sims),
                            mat_s = count_sims$mat_s,
                            gpp_s = count_sims$gpp_s,
                            year = as.integer(as.factor(count_sims$years)),
                            n_years = length(unique(count_sims$years)),
                            n_sites = length(unique(count_sims$sites)),
                            site = count_sims$sites,
                            counts = count_sims$counts,
                            x = count_sims$x,
                            xmin = count_sims$xmin,
                            xmax = count_sims$xmax)

fit_interaction = sampling(object = stan_spectra_mod_gpp_x_temp, 
                   data = stan_data_interaction,
                   iter = 200, chains = 1, cores = 4)

saveRDS(fit_interaction, file = "models/fit_interaction.rds")


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



