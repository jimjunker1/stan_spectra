library(tidyverse)
library(rstan)
library(brms)
library(tidybayes)

n_sim = 1000
xmax = 10000
xmin = 0.001
b = -1.5
u = runif(n_sim, min = 0, max = 1)

sims = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))

sim_data <- tibble(x = round(sims,3)) %>% 
  count(x, name = "counts") %>% 
  mutate(xmin = xmin,
         xmax = xmax)

fit_b_pareto_sim <- stan(file = "models/b_paretocounts_singlesample.stan",
                         data = list(x = sim_data$x,
                                     N = nrow(sim_data),
                                     counts = sim_data$counts,
                                     xmax = sim_data$xmax,
                                     xmin = sim_data$xmin),
                         iter = 1000,
                         chains = 2) 

saveRDS(fit_b_pareto_sim, file = "models/fit_b_pareto_sim.rds")

fit_b_pareto_sim = readRDS(file = "models/fit_b_pareto_sim.rds")

summary_posts = as_draws_df(fit_b_pareto_sim) %>% 
  summarize(b = mean(b_exp),
            sd = sd(b_exp))
  
b_sim = rnorm(10, mean = summary_posts$b[1],
          sd = summary_posts$sd[1])


y_rep_data = tibble(b = b_sim) %>% 
  expand_grid(u = runif(10000, min = 0, max = 1)) %>% 
  mutate(sim = as.integer(as.factor(b)),
         xmin = xmin,
         xmax = xmax,
         x = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1)),
         data = "y_rep") %>% 
  mutate(x = round(x, 3)) %>% 
  group_by(sim, data) %>% 
  count(x, name = "counts") %>% 
  bind_rows(sim_data %>% mutate(sim = 0, data = "y_raw") %>% select(sim, x, counts, data))

library(ggthemes)
post_pred_sim = y_rep_data %>% 
  ggplot(aes(x = x*counts)) + 
  geom_histogram(aes(fill = data), bins = 100) + 
  facet_wrap(~sim, scales = "free_y") +
  scale_x_log10() +
  theme_default() +
  scale_fill_colorblind() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) +
  labs(fill = "",
       x = "Individual Dry Mass per square meter",
       title = "Posterior Predicive Check with simulated data") +
  NULL

saveRDS(post_pred_sim, file = "plots/post_pred_sim.rds")


y_rep_data <- sim_data %>% 
  expand_grid(b = b_sim) %>% 
  mutate(u = runif(nrow(.), min = 0, max = 1)) %>% 
  mutate(x = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1)),
         data = "y_rep",
         sim = as.integer(as.factor(b))) %>% 
  select(x, sim, counts, xmin, xmax, data) %>% 
  bind_rows(sim_data %>% mutate(data = "y",
                                sim = 0))

y_rep_data  %>% 
  ggplot(aes(x = sim, color = data, y = x*counts)) +
  geom_jitter(width = 0.1, height = 0) + 
  scale_y_log10() + 
  geom_violin(aes(group = sim))

y_rep_data  %>% 
  ggplot(aes(x = x*counts, fill = data, group = sim)) +
  geom_histogram(bins = 50) + 
  facet_wrap(~sim) +
  scale_x_log10() +
  NULL

# repeat with sizeSpectra data --------------------------------------------

ibts = IBTS_data %>% 
  filter(Year == min(Year)) %>% 
  mutate(x = bodyMass,
         xmax = max(bodyMass),
         xmin = min(bodyMass),
         counts = Number)


fit_b_pareto_ibts <- stan(file = "models/b_paretocounts_singlesample.stan",
                         data = list(x = ibts$x,
                                     N = nrow(ibts),
                                     counts = ibts$counts,
                                     xmax = ibts$xmax,
                                     xmin = ibts$xmin),
                         iter = 1000,
                         chains = 2) 

summary_posts_ibts = as_draws_df(fit_b_pareto_ibts) %>% 
  summarize(b = mean(b_exp),
            sd = sd(b_exp))

b = rnorm(10, mean = summary_posts_ibts$b[1],
          sd = summary_posts_ibts$sd[1])

y_rep_data <- ibts %>% 
  expand_grid(b = b) %>% 
  mutate(u = runif(nrow(.), min = 0, max = 1)) %>% 
  mutate(x = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1)),
         data = "y_rep",
         sim = as.integer(as.factor(b))) %>% 
  select(x, sim, Number, xmin, xmax, data) %>% 
  bind_rows(ibts %>% mutate(data = "y",
                                sim = 0) %>% 
              select(x, sim, Number, xmin, xmax, data)) %>% 
  mutate(counts = as.integer(Number/min(Number))) %>% 
  filter(counts <= 1000) %>% 
  uncount(counts)

y_rep_data  %>% 
  ggplot(aes(x = sim, color = data, y = x*Number)) +
  # geom_jitter(width = 0.1, height = 0) + 
  scale_y_log10() + 
  geom_violin(aes(group = sim))

y_rep_data  %>% 
  ggplot(aes(x = x*Number, color = data, group = sim)) +
  geom_density() + 
  # scale_x_log10() +
  NULL
