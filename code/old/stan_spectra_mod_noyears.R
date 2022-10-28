library(tidyverse)
library(brms)
library(rstan)
library(tidybayes)
library(janitor)
library(sizeSpectra)
library(ggridges)

# 1) load data ----------------------------------------
macro_fish_dw <- readRDS("data/macro_fish_dw.rds") 

# 2) fit model ------------------------------------------------------------
# prior predictive
N = 1000
prior_sim <- tibble(beta = rnorm(N, 0, 0.1),
                    sigma_year = abs(rnorm(N, 0, 0.1)),
                    sigma_site = abs(rnorm(N, 0, 0.1)),
                    alpha_year_raw = rnorm(N, 0, 5),
                    alpha_site_raw = rnorm(N, 0, 5),
                    a = rnorm(N, -1.5, 0.2),
                    .draw = 1:N) %>% 
  mutate(intercept = a + sigma_year*alpha_year_raw + sigma_site*alpha_site_raw) %>% 
  expand_grid(mat_s = seq(-2, 2, length.out = 20)) %>% 
  mutate(y_pred = intercept + beta*mat_s)

prior_sim %>% 
  ggplot(aes(x = mat_s, y = y_pred, group = .draw)) + 
  geom_line(alpha = 0.2)


# fit full model
macro_fish_thin = macro_fish_dw %>% 
  group_by(ID) %>% 
  mutate(xmin = min(dw),
         xmax = max(dw)) %>% 
  distinct() %>% 
  ungroup() %>% 
  # sample_n(size = 18000) %>%
  # slice(1:50) %>%
  filter(site_id %in% c("ARIK", "BIGC", "BLDE", "CUPE")) %>%
  {.}

site_id = macro_fish_dw %>% ungroup() %>% distinct(site_id, mat_s)
site_id_site = tibble(site_id = site_id$site_id, 
                      site = as.integer(as.factor(unique(macro_fish_dw$site_id)))) %>% 
  left_join(site_id)


mod_spectra = stan(file = "models/b_testnc_multiple_clusters_noncentered.stan", 
                data = list(N = nrow(macro_fish_thin),
                            mat_s = macro_fish_thin$mat_s,
                            year = as.integer(as.factor(macro_fish_thin$year)),
                            n_years = length(unique(macro_fish_thin$year)),
                            n_sites = length(unique(macro_fish_thin$site_id)),
                            site = as.integer(as.factor(macro_fish_thin$site_id)),
                            counts = macro_fish_thin$no_m2,
                            x = macro_fish_thin$dw,
                            xmin = macro_fish_thin$xmin,
                            xmax = macro_fish_thin$xmax),
                iter = 500, chains = 1, cores = 4)


# saveRDS(mod_spectra, file = "models/mod_spectra.rds")


bayes_intercepts <- mod_spectra %>% 
  as_draws_df() %>% 
  pivot_longer(cols = contains("raw_site")) %>% 
  mutate(site = as.integer(parse_number(name)),
         model = "bayes") %>% 
  # left_join(site_id_site) %>%
  mutate(offset = sigma_site*value) %>%
  mutate(a_site = a + beta*mat_s) %>%
  group_by(name, site_id, mat_s) %>% 
  median_qi(a_site) 


bayes_intercepts %>% 
  ggplot(aes(x = mat_s, y = a_site)) + 
  geom_pointrange(aes(ymin = .lower, ymax = .upper),color = "red") +
  ylim(-3, 0)

bayes_draws <- mod_spectra %>% 
  as_draws_df() %>% 
  pivot_longer(cols = contains("raw_site")) %>% 
  mutate(offset = sigma_site*value) %>%
  mutate(a_site = a + offset) %>% 
  mutate(site = as.integer(parse_number(name)),
         model = "bayes") %>% 
  left_join(site_id_site)

bayes_draws %>% 
  ggplot(aes(x = site_id, y = a_site)) +
  geom_violin() 



# Plot --------------------------------------------------------------------

### Simulate regression lines

# 1) To make regression lines, simulate sequence between xmin and xmax for each sample,
# then join raw data identifiers and posteriors estimates of b (a_site, upper, lower, etc)...

# function to generate log sequence (other wise all of the numbers are too large for x to plot correctly)
# https://stackoverflow.com/questions/23901907/create-a-log-sequence-across-multiple-orders-of-magnitude
# logarithmic spaced sequence

lseq <- function(from=xmin, to=xmax, length.out=7) {
  exp(seq(log(from), log(to), length.out = length.out))
}

sim_grid <- macro_fish_thin %>%
  distinct(group, xmin, xmax) %>% 
  pivot_longer(cols = c(xmin, xmax)) %>%
  arrange(group, value) %>% 
  group_by(group) %>% 
  complete(value = lseq(min(value), max(value), length.out = 800)) %>% 
  select(-name) %>% 
  rename(x = value) %>% 
  left_join(macro_fish_thin %>% distinct(group, ID, site_id, year_month, xmin, xmax)) %>%  
  left_join(bayes_intercepts %>% 
              select(site_id, a_site, .lower, .upper)) 


# 2) ...then simulate prob x >= x for each x. y_plb estimates come from line 155 here: https://github.com/andrew-edwards/fitting-size-spectra/blob/master/code/PLBfunctions.r


line_sim <- sim_grid %>% 
  mutate(y_plb_med = 1 - (x^(a_site+1) - xmin^(a_site+1))/(xmax^(a_site+1) - xmin^(a_site+1)), # simulate prob x>=x
         y_plb_lower = 1 - (x^(.lower+1) - xmin^(.lower+1))/(xmax^(.lower+1) - xmin^(.lower+1)),
         y_plb_upper = 1 - (x^(.upper+1) - xmin^(.upper+1))/(xmax^(.upper+1) - xmin^(.upper+1))) %>%
  filter(y_plb_med > 0) %>% 
  filter(y_plb_lower > 0) %>% 
  filter(y_plb_upper > 0) %>% 
  arrange(group, x) %>%
  mutate(x = round(x, 5)) %>% 
  distinct(group, x, .keep_all = T)

### Simulate raw data (i.e., data that "would" have been collected after accounting for no_m2)

# 3) To simulate raw data, get counts of organisms and cumulative counts...
dat_bayes_counts = macro_fish_thin %>% 
  mutate(group = paste(site_id, year_month, sep = "_")) %>% 
  select(dw, no_m2, group) %>% ungroup() %>% 
  group_by(dw, group) %>% 
  summarize(Count = sum(no_m2)) %>% 
  arrange(group, desc(dw)) %>% 
  group_by(group) %>% 
  mutate(cumSum = cumsum(Count),
         cumProp = cumSum / sum(Count),
         length = ceiling(sum(Count))) 

# 4) then generate sequence of values to simulate over
dat_bayes_sim <- dat_bayes_counts %>% 
  dplyr::group_by(group, length) %>% 
  dplyr::summarize(min_cumProp = min(cumProp)) %>% 
  dplyr::group_by(group) %>% 
  dplyr::do(dplyr::tibble(cumPropsim = seq(.$min_cumProp, 1, length = .$length/10))) # dividing by something reduces file size by limiting iterations, but check for accuracy

# 5) then simulate cumulative proportion data to plot against MLE estimates by group
# make lists first
dat_bayes_simlist <- dat_bayes_sim %>% dplyr::group_by(group) %>% dplyr::group_split() 
dat_bayes_countslist <- dat_bayes_counts %>% dplyr::group_by(group) %>% dplyr::group_split() 
bayes_sim = list() # empty list to population

# simulate data with for loop
for(i in 1:length(dat_bayes_simlist)){
  bayes_sim[[i]] = dat_bayes_simlist[[i]] %>% dplyr::as_tibble() %>% 
    dplyr::mutate(dw = dat_bayes_countslist[[i]][findInterval(dat_bayes_simlist[[i]]$cumPropsim,
                                                                     dat_bayes_countslist[[i]]$cumProp), ]$dw)
}

# 6) Create data frame with "raw" data to plot
bayes_sim_tibble <- dplyr::bind_rows(bayes_sim) # dots to plot...very large file


# 7) Make plot

line_sim %>%  
  ggplot(aes(x = x, y = y_plb_med)) + 
  geom_line() +
  geom_ribbon(aes(ymin = y_plb_lower, ymax = y_plb_upper), alpha = 0.2) +
  # scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~group, scales = "free") +
  geom_point(data = bayes_sim_tibble , aes(x = dw, y = cumPropsim)) +
  NULL




# Posterior Predictive ----------------------------------------------------

# simulate y_pred from the posterior
sim_macro <- macro_fish_thin %>% 
  left_join(bayes_intercepts %>% 
              select(site_id, a_site, .lower, .upper)) %>%
  expand_grid(sim = 1:10) %>% 
  mutate( u = runif(nrow(.)),
          b = a_site) %>% 
  mutate(y_pred = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))) %>% # simulate y_pred (via Edwards github for rPLB) - confirmed in code/rplb_by_hand.R
  # mutate(y_pred = no_m2*y_pred) %>% 
  group_by(ID) %>% 
  mutate(rank = rank(y_pred)) %>% 
  mutate(group = paste(site_id, year_month, sep = "_"))

# combine with y_raw
y_pred = sim_macro %>% 
  group_by(ID, y_pred, site_id, sim) %>%
  mutate(y_pred = round(y_pred, 8)) %>%
  select(ID,site_id, y_pred, sim) %>%
  # count() %>% 
  mutate(model = "y_pred") %>% 
  rename(dw = y_pred) %>% 
  # bind_rows(macro_fish_thin %>%
  #             select(ID,site_id, dw, no_m2) %>%
  #             mutate(model = "y_raw",
  #                    sim = 0,
  #                    dw = dw)) %>%
  bind_rows(bayes_sim_tibble %>% 
              mutate(model = "y_bayes_sim",
                     sim = -1,
                     site_id = str_sub(group, 1, 4)) %>% 
              mutate(ID = as.integer(as.factor(group)))) %>% 
  right_join(macro_fish_thin %>% distinct(site_id, ID))

# plot
y_pred %>% 
  filter(dw != 0) %>% 
  # mutate(ID = as.factor(ID),
  #        sim = as.factor(sim)) %>% 
  ggplot(aes(y = dw, fill = model, x = sim)) + 
  geom_boxplot(aes(group = interaction(model, sim)),
               outlier.shape = NA) + 
  # geom_jitter(width = 0.1, height = 0, size = 0.1) +
  facet_wrap(~site_id) +
  scale_y_log10() +
  NULL


