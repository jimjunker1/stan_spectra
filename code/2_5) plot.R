library(sizeSpectra)
library(tidyverse)
library(viridis)

id_summaries %>% 
  ggplot(aes(x = mat_s, y = lambda)) + 
  geom_pointrange(aes(ymin = .lower, ymax = .upper), 
                  position = position_jitter(width = 0.04)) +
  geom_line(data = temp_x_summaries, aes(color = log_gpp_s, group = log_gpp_s)) + 
  geom_ribbon(data = temp_x_summaries, aes(ymin = .lower, ymax = .upper, fill = log_gpp_s, group = log_gpp_s),
              alpha = 0.2) + 
  # scale_color_viridis() + 
  theme_default()


gpp_x_summaries %>% 
  ggplot(aes(x = log_gpp_s, y = lambda)) + 
  # geom_pointrange(aes(ymin = .lower, ymax = .upper), 
  #                 position = position_jitter(width = 0.04)) +
  geom_line(data = gpp_x_summaries, aes(color = mat_s, group = mat_s)) + 
  geom_ribbon(data = gpp_x_summaries, aes(ymin = .lower, ymax = .upper, fill = mat_s, group = mat_s),
              alpha = 0.2) + 
  # scale_color_viridis() + 
  theme_default() + 
  facet_wrap(~mat_s)


as_draws_df(mod_spectra) %>% 
  select(beta_mat) %>% 
  filter(beta_mat > 0)

# There is a 94% probability that gpp is positively related to lambda.



# refit without varying intercepts
# precompile code
stan_dat = stan_model("models/old/b_paretocounts_singlesample.stan")

sim_data_list = macro_fish_mat_siteminmax %>% 
  group_by(site_id_int) %>% 
  group_split()

single_fit_posts = list()

for(i in 1:length(sim_data_list)){
  group = unique(sim_data_list[[i]]$site_id_int)
  
  stan_dat = list(x = sim_data_list[[i]]$x,
                  N = nrow(sim_data_list[[i]]),
                  counts = sim_data_list[[i]]$counts,
                  # mat_s = sim_data[[i]]$predictor,
                  # n_years = length(unique(sim_data[[i]]$id)),
                  # year = as.integer(as.factor(sim_data[[i]]$id)),
                  xmax = sim_data_list[[i]]$xmax,
                  xmin = sim_data_list[[i]]$xmin)
  
  single_fits = sampling(object = fit_model,
                         data = stan_dat,
                         chains = 1, 
                         iter = 1000)
  
  single_fit_posts[[i]] = as_draws_df(single_fits) %>% 
    mutate(group = group)
}


single_vs_all = bind_rows(single_fit_posts) %>% 
  group_by(group) %>% 
  median_qi(lambda) %>% 
  mutate(model = "single",
         ID = group) %>% 
  left_join(id_summaries %>% ungroup %>% distinct(mat_s, log_gpp_s,site_id, ID)) %>% 
  bind_rows(id_summaries %>% mutate(model = "full"))


single_vs_all %>% 
  ggplot(aes(x = mat_s, y = lambda, ymin = .lower, ymax = .upper, color = model)) +
  geom_pointrange(position = position_dodge(width = 0.1))


# biomass -----------------------------------------------------------------
biomass = readRDS(file = "data/biomass.rds")

biomass_posts = readRDS(file = "posteriors/biomass_posts.rds")

biomass_summary = biomass_posts %>% 
  group_by(mat_s, log_gpp_s) %>% 
  median_qi(.epred)


biomass_summary %>% 
  filter(log_gpp_s > min(log_gpp_s & log_gpp_s < max(log_gpp_s))) %>% 
  ggplot(aes(x = mat_s, y = .epred)) + 
  geom_line() + 
  geom_point(data = biomass, aes(y = sample_biomass)) + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2) + 
  scale_y_log10()


