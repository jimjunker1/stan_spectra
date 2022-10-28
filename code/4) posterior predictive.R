library(rstan)
library(tidyverse)
library(janitor)
library(tidybayes)
library(brms)

macro_fish_thin <- readRDS("data/macro_fish_thin.rds") %>% ungroup 
mod_spectra <- readRDS("models/mod_spectra.rds")
bayes_sim_tibble = readRDS("data/bayes_sim_tibble.rds")

# get posts
summary_posts = as_draws_df(mod_spectra) %>% 
  summarize(b = mean(a),
            sd = sd(a))

# draw b's from posterior
b = rnorm(10, mean = summary_posts_ibts$b[1],
          sd = summary_posts_ibts$sd[1])


# simulate the full data set ----------------------------------------------

y_reps = macro_fish_thin %>% 
  mutate(u = runif(nrow(.), min = 0, max = 1)) %>% # uniform draw
  expand_grid(b = b) %>%                           # add b draws
  mutate(y_rep = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1)), # simulate
         sim = as.integer(as.factor(b))) %>% # add identifier
  select(y_rep, sim, no_m2) %>% 
  mutate(data = "y_rep") %>% 
  rename(y = y_rep) %>%
  bind_rows(macro_fish_thin %>% mutate(sim = 0,
                                       data = "y_raw") %>% 
              select(data, sim, dw, no_m2) %>% 
              rename(y = dw)) 

y_reps %>% 
  ggplot(aes(x = y*no_m2, fill = data, group = sim)) + 
  geom_histogram(bins = 50) +
  facet_wrap(~sim) +
  scale_x_log10() +
  NULL


macro_fish_thin %>% 
  mutate(u = runif(nrow(.), min = 0, max = 1)) %>% # uniform draw
  expand_grid(b = b) %>%                           # add b draws
  mutate(y_rep = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1)), # simulate
         sim = as.integer(as.factor(b))) %>% 
  filter(sim <= 1) %>% 
  ggplot(aes(x = y_rep*no_m2, y = dw*no_m2)) +
  geom_point(size = 0.1) + 
  scale_x_log10() +
  scale_y_log10() + 
  geom_abline()


macro_fish_thin %>% 
  mutate(u = runif(nrow(.), min = 0, max = 1)) %>% # uniform draw
  expand_grid(b = b) %>%                           # add b draws
  mutate(y_rep = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1)), # simulate
         sim = as.integer(as.factor(b)),
         obs_exp = y_rep - dw) %>% 
  group_by(sim) %>% 
  mutate(greater_than = case_when(obs_exp > 0 ~ 1, 
                                  TRUE ~ 0)) %>% 
  summarize(p_value = sum(greater_than)/nrow(.))


# simulate by site ----------------------------------------------
macro_fish_thin %>% 
  mutate(year = as.integer(as.factor(macro_fish_thin$year))) %>% 
  distinct(site_id_int, year, xmin, xmax, no_m2) %>% 
  rename(counts = no_m2)

