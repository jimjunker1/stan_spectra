library(tidyverse)
library(rstan)
library(brms)
library(tidybayes)

# How many individuals? ----------------------------------------------------------------
# make empty tibble 
sample_size_sims = tibble(iter = numeric(),
                    n = numeric(),
                    b = numeric(),
                    sd = numeric())
# precompile code
fit_model = stan_model("models/old/b_paretocounts_singlesample.stan")

# iterate over 1:10 replicates for each of 11 sample sizes
for (i in 1:10) {
  for (j in 2^seq(1, 11)) {
    for(k in seq(-2, -1.2, length.out = 3)) {
    n_sim = 3000
    xmax = 215000
    xmin = 0.003
    b = k
    u = runif(n_sim, min = 0, max = 1)
    
    sims = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))
    
    
    sim_data <- tibble(x = round(sims,3)) %>% 
      sample_n(j) %>% 
      count(x, name = "counts") %>% 
      mutate(xmin = xmin,
             xmax = xmax)
    
    stan_dat <- list(x = sim_data$x,
                     N = nrow(sim_data),
                     counts = sim_data$counts,
                     xmax = sim_data$xmax,
                     xmin = sim_data$xmin)
    
    fit <- sampling(object = fit_model,
                    data = stan_dat,
                    iter = 1000,
                    chains = 2,
                    open_progress = F,
                    verbose = F)
    
    b_exp = as.data.frame(rstan::extract(fit, pars = "b_exp")) 
    
    sample_size_sims <- sample_size_sims %>%
      add_row(iter = i,
              n = j,
              b = mean(b_exp$b_exp),
              sd = sd(b_exp$b_exp))
  }
  }
}

saveRDS(sample_size_sims, "models/sample_size_simulations/sample_size_sims.rds")

sample_size_sims = readRDS("models/sample_size_simulations/sample_size_sims.rds")

b_lines = tibble(b = c(-2, -1.6, -1.2),
                 b_known = c(-2, -1.6, -1.2),
                 sd = NA) %>% 
  pivot_longer(cols = -b_known)

sample_size_sims %>% 
  mutate(b_known = rep(rep(c(-2, -1.6, -1.2),11),10)) %>% 
  pivot_longer(cols = c(b, sd)) %>% 
  ggplot(aes(x = n, y = value)) + 
  facet_wrap(b_known~name, scales = "free", ncol = 2) +
  geom_point() +
  labs(x = "Number of individuals measured (before combining to densities)") + 
  scale_x_log10() +
  geom_hline(data = b_lines, aes(yintercept = value))


sample_size_plot = sample_size_sims %>% 
  mutate(b_known = rep(rep(c(-2, -1.6, -1.2),11),10)) %>% 
  pivot_longer(cols = c(b, sd)) %>% 
  filter(name == "b") %>% 
  ggplot(aes(x = n, y = value)) + 
  facet_wrap(~b_known, scales = "free", ncol = 3) +
  geom_point(size = 0.1) +
  labs(x = "Number of individual body sizes in a sample",
       y = "lambda") + 
  scale_x_log10() +
  geom_hline(data = b_lines, aes(yintercept = value)) + 
  theme_default()

library(ggview)
ggview(sample_size_plot, width = 6, height = 2, units = "in")
ggsave(sample_size_plot, width = 6, height = 2, units = "in", 
       file = "plots/sample_size_plot.jpg")
saveRDS(sample_size_plot, file = "plots/sample_size_plot.rds")

sample_size_sims %>% 
  mutate(b_known = rep(rep(c(-2, -1.6, -1.2),11),10)) %>% 
  # pivot_longer(cols = c(b, sd)) %>% 
  ggplot(aes(x = n, y = b)) + 
  facet_wrap(~b_known, scales = "free", ncol = 1) +
  geom_pointrange(aes(ymin = b - sd, ymax = b + sd),position = position_jitter(width = 0.06, height = 0),
                  size = 0.2) +
  labs(x = "Number of individuals measured (before combining to densities)") + 
  scale_x_log10() +
  geom_hline(data = b_lines, aes(yintercept = value))





# How many size/density groups needed? ----------------------------------------------------------------
# make empty tibble 
sample_sizedensity_sims = tibble(iter = numeric(),
                          n = numeric(),
                          b = numeric(),
                          sd = numeric())
# precompile code
fit_model = stan_model("models/old/b_paretocounts_singlesample.stan")

# iterate over 1:10 replicates for each of 11 sample sizes
for (i in 1:10) {
  for (j in seq(2,60, length.out = 8)) {
    for(k in seq(-2.2, -1.2, length.out = 3)) {
      n_sim = 3000
      xmax = 10000
      xmin = 0.001
      b = k
      u = runif(n_sim, min = 0, max = 1)
      
      sims = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))
      
      
      sim_data <- tibble(x = round(sims,3)) %>% 
        # sample_n(j) %>%
        count(x, name = "counts") %>% 
        mutate(xmin = xmin,
               xmax = xmax) %>% 
        sample_n(j)
      
      stan_dat <- list(x = sim_data$x,
                       N = nrow(sim_data),
                       counts = sim_data$counts,
                       xmax = sim_data$xmax,
                       xmin = sim_data$xmin)
      
      fit <- sampling(object = fit_model,
                      data = stan_dat,
                      iter = 1000,
                      chains = 2,
                      open_progress = F,
                      verbose = F)
      
      b_exp = as.data.frame(extract(fit, pars = "b_exp")) 
      
      sample_sizedensity_sims <- sample_sizedensity_sims %>%
        add_row(iter = i,
                n = j,
                b = mean(b_exp$b_exp),
                sd = sd(b_exp$b_exp))
    }
  }
}

saveRDS(sample_sizedensity_sims, "models/sample_size_simulations/sample_sizedensity_sims.rds")

b_lines = tibble(b = c(-2.2, -1.7, -1.2),
                 b_known = c(-2.2, -1.7, -1.2),
                 sd = NA) %>% 
  pivot_longer(cols = -b_known)

sample_sizedensity_sims %>% 
  mutate(b_known = rep(rep(c(-2.2, -1.7, -1.2),8),10)) %>% 
  pivot_longer(cols = c(b, sd)) %>% 
  ggplot(aes(x = n, y = value)) + 
  facet_wrap(b_known~name, scales = "free", ncol = 2) +
  geom_point() +
  labs(x = "Number of individuals measured (before combining to densities)") + 
  scale_x_log10() +
  geom_hline(data = b_lines, aes(yintercept = value))


sample_sizedensity_sims %>% 
  mutate(b_known = rep(rep(c(-2.2, -1.7, -1.2),8),10)) %>% 
  # pivot_longer(cols = c(b, sd)) %>% 
  ggplot(aes(x = n, y = b)) + 
  facet_wrap(~b_known, scales = "free", ncol = 1) +
  geom_pointrange(aes(ymin = b - sd, ymax = b + sd),position = position_jitter(width = 0.006, height = 0),
                  size = 0.2) +
  labs(x = "Number of individuals measured (before combining to densities)") + 
  scale_x_log10() +
  geom_hline(data = b_lines, aes(yintercept = value))


# How many size/density groups for NEON data? -----------------------------

# make stan data
neon_data = readRDS(file = "data/macro_fish_mat_siteminmax.rds")

# load fitted model
mod_spectra_siteminmax = readRDS(file = "models/sandbox/mod_spectra_siteminmax.rds")

# make empty tibble 

colnames = names(as.data.frame(extract(mod_spectra_siteminmax)))

sample_sizedensity_sims_neon = data.frame(matrix(ncol = 34)) %>% 
  mutate_all(as.numeric)

colnames(sample_sizedensity_sims_neon) = colnames

sample_sizedensity_sims_neon = sample_sizedensity_sims_neon %>% mutate(iter = as.numeric(NA),
                                                                       n = as.numeric(NA),
                                                                       statistic = NA)

# precompile code
fit_neon_model = stan_model("models/stan_spectra_mod.stan")

# iterate over 1:10 replicates for each of 11 sample sizes
for (i in 1:2) {
  for (j in seq(0.1, 0.5, length.out = 3)) {
      
      data = neon_data %>% 
        group_by(site_id, year) %>% 
        sample_frac(j)
      
      stan_dat = list(N = nrow(data),
                       mat_s = data$mat_s,
                       year = as.integer(as.factor(data$year)),
                       n_years = length(unique(data$year)),
                       n_sites = length(unique(data$site_id_int)),
                       site = data$site_id_int,
                       counts = data$no_m2,
                       x = data$dw,
                       xmin = data$xmin,
                       xmax = data$xmax)
      
      fit <- sampling(object = fit_neon_model,
                      data = stan_dat,
                      iter = 1000,
                      chains = 1,
                      open_progress = F,
                      verbose = F)
      
      mean = as_tibble(as.data.frame(extract(fit))) %>% 
        lapply(mean) %>% as_tibble() %>% 
        mutate(iter = i, n = j, 
               statistic = "mean")
      
      sd = as_tibble(as.data.frame(extract(fit))) %>% 
        lapply(sd) %>% as_tibble() %>% 
        mutate(iter = i, n = j,
               statistic = "sd")
      
      sample_sizedensity_sims_neon <- bind_rows(sample_sizedensity_sims_neon, mean, sd)
      }
}

mean_full = as_tibble(as.data.frame(extract(mod_spectra_siteminmax))) %>% 
  lapply(mean) %>% as_tibble() %>% 
  mutate(iter = 1, n = 1, 
         statistic = "mean")

sd_full = as_tibble(as.data.frame(extract(mod_spectra_siteminmax))) %>% 
  lapply(sd) %>% as_tibble() %>% 
  mutate(iter = 1, n = 1,
         statistic = "sd")


sample_sizedensity_sims_neon = bind_rows(sample_sizedensity_sims_neon, mean_full, sd_full)
saveRDS(sample_sizedensity_sims_neon, "models/sample_size_simulations/sample_sizedensity_sims_neon.rds")


sample_sizedensity_sims_neon %>% 
  select(-lp__) %>% 
  janitor::remove_empty(which = "rows") %>% 
  pivot_longer(cols = c(-statistic, -iter, -n)) %>% 
  pivot_wider(names_from = statistic, values_from = value) %>% 
  ggplot(aes(x = n, y = sd)) + 
  geom_jitter(height = 0, width = 0.02)


a = sample_sizedensity_sims_neon %>% 
  filter(n != 1) %>% 
  select(-lp__) %>% 
  janitor::remove_empty(which = "rows") %>% 
  pivot_longer(cols = c(-statistic, -iter, -n)) %>% 
  pivot_wider(names_from = statistic, values_from = value) %>% 
  rename(mean_partial = mean,
         sd_partial = sd,
         n_partial = n)

b = sample_sizedensity_sims_neon %>% 
  filter(n == 1) %>% 
  select(-lp__) %>% 
  janitor::remove_empty(which = "rows") %>% 
  pivot_longer(cols = c(-statistic, -iter, -n)) %>% 
  pivot_wider(names_from = statistic, values_from = value) %>% 
  rename(mean_full = mean,
         sd_full = sd,
         n_full = n)

compare_params = left_join(a, b)

compare_params %>% 
  ggplot(aes(y = mean_partial, x = mean_full)) + 
  geom_point() + 
  facet_wrap(~n_partial, scales = "free")

compare_params %>% 
  ggplot(aes(y = sd_partial, x = sd_full)) + 
  geom_point() + 
  facet_wrap(~n_partial, scales = "free")



