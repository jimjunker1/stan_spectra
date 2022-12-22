library(rstan)

count_sims = readRDS("data/count_sims.rds")
stan_spectra_mod_gpp_x_temp = stan_model("models/stan_spectra_mod_gpp_x_temp.stan")

stan_data_interaction = list(N = nrow(count_sims),
                             mat_s = count_sims$mat_s,
                             gpp_s = count_sims$gpp_s,
                             year = as.integer(as.factor(count_sims$year)),
                             n_years = length(unique(count_sims$year)),
                             n_sites = length(unique(count_sims$site_id_int)),
                             site = count_sims$site_id_int,
                             counts = count_sims$counts,
                             x = count_sims$x,
                             xmin = count_sims$xmin,
                             xmax = count_sims$xmax)

fit_interaction = sampling(object = stan_spectra_mod_gpp_x_temp, 
                           data = stan_data_interaction,
                           iter = 200, chains = 2, cores = 4)

saveRDS(fit_interaction, file = "models/fit_interaction.rds")