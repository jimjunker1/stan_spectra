library(rstan)

rstan_options(autowrite = TRUE)
rstan_options(threads_per_chain = 1)
options(mc.cores = parallel::detectCores())
# load data

# compile model
# stan_spectra_mod_gpp_x_temp = stan_model("models/stan_spectra_mod_gpp_x_temp.stan")
stan_model = stan_model("models/stan_spectra_mod_temponly.stan")


dat = readRDS(file = "data/macro_fish_mat_siteminmax.rds")  

stan_data_interaction = list(N = nrow(dat),
                             mat_s = dat$mat_s,
                             # gpp_s = dat$log_gpp_s,
                             year = as.integer(as.factor(dat$year)),
                             n_years = length(unique(dat$year)),
                             n_samples = length(unique(dat$sample_id_int)),
                             n_sites = length(unique(dat$site_id_int)),
                             sample = dat$sample_id_int,
                             site = dat$site_id_int,
                             counts = dat$counts,
                             x = dat$x,
                             xmin = dat$xmin,
                             xmax = dat$xmax)

fit = sampling(object = stan_model,
               data = stan_data_interaction,
               iter = 2000, chains = 4, cores = 4)

saveRDS(fit, file = "models/fit_temperature.rds")

