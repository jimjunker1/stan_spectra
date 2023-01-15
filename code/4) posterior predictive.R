library(rstan)
library(tidyverse)
library(janitor)
library(tidybayes)
library(brms)

interaction_posts_wrangled = readRDS(file = "posteriors/interaction_posts_wrangled.rds")
temp_posts_wrangled = readRDS(file = "posteriors/temp_posts_wrangled.rds")
dat = readRDS(file = "data/macro_fish_mat_siteminmax.rds") 

sim_postpreds = function(model, sample_n = 600000){
  interaction_posts_wrangled %>% 
    filter(.draw <= 100) %>% 
    right_join(dat) %>% 
    mutate(u = runif(nrow(.), 0, 1)) %>% # uniform draw
    mutate(x = (u*xmax^(lambda+1) +  (1-u) * xmin^(lambda+1) ) ^ (1/(lambda+1))) %>% 
    mutate(data = "y_rep") %>% 
    bind_rows(dat %>% 
                mutate(data = "y_raw") %>% 
                mutate(.draw = 0)) %>% 
    rename(sim = x) %>% 
    group_by(.draw) %>% 
    sample_n(sample_n, weight = counts, replace = T) %>% 
    select(sim, .draw, data, site_id, year, sample_id_int, xmin, xmax, counts)
  
}

id = as.integer(runif(1, 1, 118))

sim_postpreds(interaction_posts_wrangled) %>%
  filter(sample_id_int == id) %>% 
  ggplot(aes(x = .draw, y = sim*counts, color = data, group = .draw)) + 
  geom_violin() +
  scale_y_log10() +
  facet_wrap(sample_id_int ~ site_id) +
  NULL

simposts = sim_postpreds(temp_posts_wrangled)

simposts %>% 
  mutate(data = as.factor(data)) %>% 
  mutate(data = fct_relevel(data, "y_rep")) %>% 
  ggplot(aes(x = sim*counts, color = data, group = .draw, alpha = data)) +
  geom_density() +
  scale_x_log10() +
  facet_wrap(~sample_id_int, scales = "free_y") +
  scale_alpha_manual(values = c(0.01, 0.01)) +
  ggthemes::scale_color_colorblind() +
  theme_default() +
  NULL

median = dat %>% 
  # group_by(sample_id_int) %>%
  summarize(median = mean(x*counts))


simposts_median = simposts %>% 
  group_by(.draw) %>% 
  summarize(median = mean(sim*counts))


simposts_median %>% 
  # filter(sample_id_int <= 30) %>%
  ggplot(aes(x = median)) + 
  geom_histogram() + 
  # facet_wrap(~sample_id_int, scales = "free_x") +
  geom_vline(data = median, aes(xintercept = median)) +
  scale_x_log10() +
  NULL




# main effects only -------------------------------------------------------

fit_interaction = readRDS(file = "models/fit_interaction.rds")

posts = as_draws_df(fit_interaction) %>% 
  expand_grid(dat %>% ungroup %>% distinct(mat_s, log_gpp_s)) %>% 
  mutate(lambda = a + beta_mat*mat_s + beta_gpp*log_gpp_s + beta_gpp_mat*log_gpp_s*mat_s) 

post_preds = posts %>% 
  filter(.draw <= 10) %>% 
  left_join(dat, by = c("log_gpp_s", "mat_s")) %>% 
  mutate(u = runif(nrow(.), 0, 1)) %>% # uniform draw
  mutate(x = (u*xmax^(lambda+1) +  (1-u) * xmin^(lambda+1) ) ^ (1/(lambda+1))) %>% 
  mutate(data = "y_rep") %>% 
  bind_rows(dat %>% 
              mutate(data = "y_raw") %>% 
              mutate(.draw = 0)) %>% 
  rename(sim = x) %>% 
  group_by(.draw) %>% 
  sample_n(1000, weight = counts, replace = T) %>% 
  select(sim, .draw, data, site_id, year, sample_id_int)

post_preds %>% 
  ggplot(aes(x = sim, color = data, group = .draw)) + 
  geom_density() +
  scale_x_log10()





