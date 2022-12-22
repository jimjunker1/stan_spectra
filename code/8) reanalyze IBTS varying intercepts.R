library(sizeSpectra)
library(rstan)
library(janitor)
library(tidyverse)
library(brms)
library(ggthemes)
library(patchwork)
library(tidybayes)
rstan_options("auto_write" = TRUE)


# reanalyze IBTS data from Edwards et al. 2020 ----------------------------

# analyze regression with varying intercepts ---------------------
data("IBTS_data")

stan_varint_only = stan_model("models/stan_varint_only.stan")

ibts_data = IBTS_data %>% 
  clean_names() %>% 
  mutate(group = as.integer(as.factor(ibts_data$year)),
         x = body_mass,
         xmin = min(x),
         xmax = max(x))

stan_dat = list(N = nrow(ibts_data),
                x = ibts_data$x,
                xmax = ibts_data$xmax,
                xmin = ibts_data$xmin,
                group = as.integer(as.factor(ibts_data$group)),
                n_group = length(unique(ibts_data$group)),
                # n_sites = length(unique(ibts_data$site_id_int)),
                # site = ibts_data$site_id_int,
                counts = ibts_data$number,
                x = ibts_data$x,
                xmin = ibts_data$xmin,
                xmax = ibts_data$xmax)

ibts_varintonly_mod = sampling(object = stan_varint_only,
                               data = stan_dat,
                               chains = 2,
                               iter = 1000,
                               cores = 2)

# priors on 12/22 run were normal(-1.5, 0.5) for alpha and exp(8) for sigma
saveRDS(ibts_varintonly_mod, file = "models/ibts_varintonly_mod.rds")

# posts
ibts_varintonly_posts = as_draws_df(ibts_varintonly_mod) %>% 
  select(-beta_mat) %>% # this was mistakenly declared as a variable, but was not in the model equation. It shouldn't have affected the varying interctps
  pivot_longer(cols = contains("alpha_raw")) %>% 
  mutate(group = parse_number(name),
         lambda = a + sigma_group*value,
         method = "Bayesian Hierarchical") %>% 
  left_join(ibts_data %>% clean_names() %>% ungroup %>% distinct(year, group)) 

ibts_varintonly_post_summary = ibts_varintonly_posts %>% 
  group_by(method, year) %>% 
  median_qi(lambda) 

edwards_ibts_years = bind_rows(edwards_mle_estimates, ibts_varintonly_post_summary) %>% 
  left_join(ibts_data %>% group_by(year) %>% tally())

edwards_ibts_years %>% 
  mutate(method = as.factor(method),
         method = fct_relevel(method, "MLE")) %>% 
  ggplot(aes(x = year, y = lambda, shape = method, color = method)) +
  geom_pointrange(aes(ymin = .lower, ymax = .upper),
                  position = position_dodge(width = 0.4)) +
  scale_color_grey(start = 0.3, end = 0.2) +
  scale_shape_manual(values = c(1, 16)) +
  geom_hline(aes(yintercept = mean(ibts_varintonly_posts$a))) + 
  theme_default() + 
  labs(x = "Year",
       y = "\u03bb",
       color = "Method",
       shape = "Method")


