library(rstan)
library(tidyverse)
library(janitor)

rstan_options(autowrite = TRUE)
rstan_options(threads_per_chain = 1)

# load data

macro_fish_mat_siteminmax = readRDS(file = "data/macro_fish_mat_siteminmax.rds") 

# compile model
stan_spectra_mod_gpp_x_temp = stan_model("models/stan_spectra_mod_gpp_x_temp.stan")
# stan_spectra_mod_gpp_x_temp_randsiteonly = stan_model("models/stan_spectra_mod_gpp_x_temp_randsiteonly.stan")


# make data and fit model ---------------------------------------------------------

dat = macro_fish_mat_siteminmax 

stan_data_interaction = list(N = nrow(dat),
                             mat_s = dat$mat_s,
                             gpp_s = dat$gpp_s,
                             year = as.integer(as.factor(dat$year)),
                             n_years = length(unique(dat$year)),
                             n_sites = length(unique(dat$site_id_int)),
                             site = dat$site_id_int,
                             counts = dat$counts,
                             x = dat$x,
                             xmin = dat$xmin,
                             xmax = dat$xmax)

fit_interaction = sampling(object = stan_spectra_mod_gpp_x_temp, 
                           data = stan_data_interaction,
                           iter = 2000, chains = 2, cores = 4)

saveRDS(fit_interaction, file = "models/fit_interaction.rds")
