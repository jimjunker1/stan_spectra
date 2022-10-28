
site_ids = macro_fish_dw %>% 
  ungroup() %>% 
  distinct(site_id) %>% 
  mutate(site_id_int = as.integer(as.factor(site_id)))


macro_fish_thin = macro_fish_dw %>% 
  group_by(dw, mat_s, year, site_id, xmin, xmax, ID) %>% 
  summarize(no_m2 = sum(no_m2)) %>% 
  left_join(site_ids)

saveRDS(macro_fish_thin, file = "data/macro_fish_thin.rds")


mod_spectra = stan(file = "models/stan_spectra_mod_test.stan", 
                   data = list(N = nrow(macro_fish_thin),
                               mat_s = macro_fish_thin$mat_s,
                               year = as.integer(as.factor(macro_fish_thin$year)),
                               n_years = length(unique(macro_fish_thin$year)),
                               n_sites = length(unique(macro_fish_thin$site_id_int)),
                               site = macro_fish_thin$site_id_int,
                               counts = macro_fish_thin$no_m2,
                               x = macro_fish_thin$dw,
                               xmin = macro_fish_thin$xmin,
                               xmax = macro_fish_thin$xmax),
                   iter = 100, chains = 1, cores = 4)


mod_spectra

site_ids = macro_fish_dw %>% 
  ungroup() %>% 
  distinct(site_id) %>% 
  mutate(site_id_int = as.integer(as.factor(site_id)))




hist(macro_fish_thin$dw, xlim = range(0, 200000))
hist(sim_data_counts$x_sim, xlim = range(0, 200000))

sim_data_counts %>% 
  ggplot(aes(x = n, y = x_sim)) + 
  geom_point() + 
  scale_x_log10()

macro_fish_thin %>% 
  ggplot(aes(x = no_m2, y = dw)) + 
  geom_point() + 
  scale_y_log10() +
  scale_x_log10() +
  NULL
