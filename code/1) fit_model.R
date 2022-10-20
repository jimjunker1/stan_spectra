library(rstan)

# 1) load data ----------------------------------------
macro_fish_mat <- readRDS("data/macro_fish_mat.rds") %>% 
  mutate(site_id_int = as.integer(as_factor(site_id)))

# 2) make stan data
stan_data = list(N = nrow(macro_fish_mat),
                 mat_s = macro_fish_mat$mat_s,
                 year = as.integer(as.factor(macro_fish_mat$year)),
                 n_years = length(unique(macro_fish_mat$year)),
                 n_sites = length(unique(macro_fish_mat$site_id_int)),
                 site = macro_fish_mat$site_id_int,
                 counts = macro_fish_mat$no_m2,
                 x = macro_fish_mat$dw,
                 xmin = macro_fish_mat$xmin,
                 xmax = macro_fish_mat$xmax)

# 3) fit model ------------------------------------------------------------
mod_spectra = stan(file = "models/stan_spectra_mod.stan", 
                   data = stan_data,
                   iter = 1000, chains = 2, cores = 4)


saveRDS(mod_spectra, file = "models/mod_spectra.rds")

mod_spectra

