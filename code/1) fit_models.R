library(rstan)

rstan_options(autowrite = TRUE)
rstan_options(threads_per_chain = 1)

# load data
macro_fish_mat_siteminmax = readRDS(file = "data/macro_fish_mat_siteminmax.rds") # this is really site max and global min
macro_fish_mat_globalminmax  = readRDS(file = "data/macro_fish_mat_globalminmax.rds") 
macro_fish_mat = readRDS(file = "data/macro_fish_mat.rds")

# fit with global xmins and local xmaxs ----------------------------------------

# make stan data
stan_data = list(N = nrow(macro_fish_mat_siteminmax),
                 mat_s = macro_fish_mat_siteminmax$mat_s,
                 year = as.integer(as.factor(macro_fish_mat_siteminmax$year)),
                 n_years = length(unique(macro_fish_mat_siteminmax$year)),
                 n_sites = length(unique(macro_fish_mat_siteminmax$site_id_int)),
                 site = macro_fish_mat_siteminmax$site_id_int,
                 counts = macro_fish_mat_siteminmax$no_m2,
                 x = macro_fish_mat_siteminmax$dw,
                 xmin = macro_fish_mat_siteminmax$xmin,
                 xmax = macro_fish_mat_siteminmax$xmax)

# fit model 
mod_spectra_siteminmax = stan(file = "models/stan_spectra_mod.stan", 
                              data = stan_data,
                              iter = 1000, chains = 2, cores = 4)


saveRDS(mod_spectra_siteminmax, file = "models/sandbox/mod_spectra_siteminmax.rds")

mod_spectra_siteminmax

# fit with global xmins/xmaxs -------------------------------------------


# make stan data
stan_data = list(N = nrow(macro_fish_mat_globalminmax),
                 mat_s = macro_fish_mat_globalminmax$mat_s,
                 year = as.integer(as.factor(macro_fish_mat_globalminmax$year)),
                 n_years = length(unique(macro_fish_mat_globalminmax$year)),
                 n_sites = length(unique(macro_fish_mat_globalminmax$site_id_int)),
                 site = macro_fish_mat_globalminmax$site_id_int,
                 counts = macro_fish_mat_globalminmax$no_m2,
                 x = macro_fish_mat_globalminmax$dw,
                 xmin = macro_fish_mat_globalminmax$xmin,
                 xmax = macro_fish_mat_globalminmax$xmax)

# fit model 
mod_spectra_globalminmax = stan(file = "models/stan_spectra_mod.stan", 
                   data = stan_data,
                   iter = 1000, chains = 2, cores = 4)


saveRDS(mod_spectra_globalminmax, file = "models/sandbox/mod_spectra_globalminmax.rds")

mod_spectra_globalminmax





# fit original model with sample level xmins/xmaxs ------------------------


# make stan data
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

# fit model 
mod_spectra = stan(file = "models/stan_spectra_mod.stan", 
                   data = stan_data,
                   iter = 1000, chains = 2, cores = 4)


saveRDS(mod_spectra, file = "models/mod_spectra.rds")

mod_spectra





