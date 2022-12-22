library(rstan)
library(tidyverse)

rstan_options(autowrite = TRUE)
rstan_options(threads_per_chain = 1)

# load data
gpp = readRDS("data/gpp_means.rds") %>% clean_names() %>% 
  rename(gpp = mean, gpp_sd = sd) %>% 
  mutate(log_gpp = log(gpp),
         gpp_s = scale(gpp, center = T, scale = T),
         gpp_s = as.numeric(gpp_s))

macro_fish_mat_siteminmax = readRDS(file = "data/macro_fish_mat_siteminmax.rds") %>% # this is really site max and global min
  left_join(gpp)

# compile models
stan_spectra_mod_counts = stan_model("models/stan_spectra_mod.stan")
stan_spectra_mod_nocounts = stan_model("models/stan_spectra_mod_nocounts.stan")

# fit counts ----------------------------------------

# make stan data

data_counts = macro_fish_mat_siteminmax 

stan_data = list(N = nrow(data_counts),
                 mat_s = data_counts$mat_s,
                 year = as.integer(as.factor(data_counts$year)),
                 n_years = length(unique(data_counts$year)),
                 n_sites = length(unique(data_counts$site_id_int)),
                 site = data_counts$site_id_int,
                 # counts = data_counts$no_m2,
                 x = data_counts$dw,
                 xmin = data_counts$xmin,
                 xmax = data_counts$xmax)

# fit model 
mod_spectra_counts = sampling(object = stan_spectra_mod_counts, 
                                data = stan_data,
                                iter = 2000, chains = 2, cores = 4)




# fit no counts ----------------------------------------

# make stan data

data_nocounts = macro_fish_mat_siteminmax %>% 
  sample_n(1000)

stan_data = list(N = nrow(data_nocounts),
                 mat_s = data_nocounts$mat_s,
                 year = as.integer(as.factor(data_nocounts$year)),
                 n_years = length(unique(data_nocounts$year)),
                 n_sites = length(unique(data_nocounts$site_id_int)),
                 site = data_nocounts$site_id_int,
                 # counts = data_nocounts$no_m2,
                 x = data_nocounts$dw,
                 xmin = data_nocounts$xmin,
                 xmax = data_nocounts$xmax)

# fit model 
mod_spectra_nocounts = sampling(object = stan_spectra_mod_nocounts, 
                              data = stan_data,
                              iter = 50, chains = 1, cores = 4)







# fit interaction ---------------------------------------------------------

count_sims = macro_fish_mat_siteminmax %>% 
  filter(!is.na(gpp_s)) %>% 
  mutate(site_id_int = as.integer(as.factor(site_id))) %>% 
  group_by(ID) %>% 
  sample_frac(1)

count_sims = readRDS("data/count_sims.rds")

stan_spectra_mod_gpp_x_temp = stan_model("models/stan_spectra_mod_gpp_x_temp.stan")
stan_spectra_mod_gpp_x_temp_randyearonly = stan_model("models/stan_spectra_mod_gpp_x_temp_randyearonly.stan")
stan_spectra_mod_gpp_x_temp_norand = stan_model("models/stan_spectra_mod_gpp_x_temp_norand.stan")

stan_data_interaction = list(N = nrow(count_sims),
                             mat_s = count_sims$mat_s,
                             gpp_s = count_sims$gpp_s,
                             year = as.integer(as.factor(count_sims$year)),
                             n_years = length(unique(count_sims$year)),
                             # n_sites = length(unique(count_sims$site_id_int)),
                             # site = count_sims$site_id_int,
                             counts = count_sims$counts,
                             x = count_sims$x,
                             xmin = count_sims$xmin,
                             xmax = count_sims$xmax)

fit_interaction = sampling(object = stan_spectra_mod_gpp_x_temp_norand, 
                           data = stan_data_interaction,
                           iter = 200, chains = 2, cores = 4)

saveRDS(fit_interaction, file = "models/fit_interaction.rds")



count_sims %>% ungroup() %>% distinct(mat_s, gpp_s, mat_site, gpp, site_id) %>% 
  ggplot(aes(x = mat_s, y = gpp_s)) +
  # geom_smooth() +
  geom_point()
