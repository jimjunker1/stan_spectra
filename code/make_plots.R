library(rstan)
library(tidyverse)
library(janitor)

macro_fish_thin <- readRDS("data/macro_fish_thin.rds") 
mod_spectra <- readRDS("models/mod_spectra.rds")
mle_mat <- read_csv("C:/Users/Jeff.Wesner/OneDrive - The University of South Dakota/USD/Github Projects/neon_size_spectra/data/derived_data/mle_mat.csv") %>% 
  select(-ID) %>% 
  left_join(macro_fish_thin %>% distinct(site_id, mat_s))
prior_pred_b_vs_mats = readRDS("plots/prior_pred_b_vs_mats.rds")


# Plot regressions
sim_regressions <- mod_spectra %>% 
  as_draws_df() %>% 
  pivot_longer(cols = contains("raw_site")) %>% 
  mutate(site = as.integer(parse_number(name)),
         model = "bayes") %>% 
  select(name, value, site, sigma_site, a, beta, .draw) %>% 
  left_join(macro_fish_thin %>% ungroup %>% distinct(site_id, site_id_int, mat_s) %>% rename(site = site_id_int)) %>%
  mutate(offset = sigma_site*value) %>%
  mutate(a_site = a + beta*mat_s) %>% 
  group_by(name, site_id, mat_s) %>% 
  median_qi(a_site) 

sim_regressions %>% 
  ggplot(aes(x = mat_s, y = a_site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2) +
  theme_ggdist() + 
  geom_point(data = mle_mat %>% distinct(b, mat_s), aes(y = b)) + 
  labs(y = "b exponent",
       x = "Standardized Temperature",
       caption = "Dots are the old b estimates from the SizeSpectra package. Regression lines are from the new Stan model")


# prior vs posterior

# make data
post_prior_sims = mod_spectra %>% 
  as_draws_df() %>% 
  pivot_longer(cols = contains("raw_site")) %>% 
  mutate(site = as.integer(parse_number(name)),
         model = "bayes") %>% 
  select(name, value, site, sigma_site, a, beta, .draw) %>% 
  left_join(macro_fish_thin %>% ungroup %>% distinct(site_id, site_id_int, mat_s) %>% rename(site = site_id_int)) %>%
  mutate(offset = sigma_site*value) %>%
  mutate(a_site = a + beta*mat_s) %>% 
  select(.draw, a_site, mat_s) %>% 
  mutate(model = "Posterior") %>% 
  rename(b_only = a_site,
         sims = .draw) %>% 
  bind_rows(prior_pred_b_vs_mats$data %>% mutate(model = "Prior") %>% 
              select(model, b_only, sims, mat_s))

# plot
post_v_prior = post_prior_sims %>% 
  mutate(model = fct_relevel(model, "Prior")) %>% 
  filter(sims <= 300) %>% 
  ggplot(aes(x = mat_s, y = b_only)) + 
  geom_line(aes(group = sims), alpha = 0.2) + 
  facet_wrap(~model) + 
  labs(y = "b exponent",
       x = "Standardized mean annual temperature (deg C)") + 
  theme_ggdist()


saveRDS(post_v_prior, file = "plots/post_v_prior.rds")


# Plot site means
sim_site_means <- mod_spectra %>% 
  as_draws_df() %>% 
  pivot_longer(cols = contains("raw_site")) %>% 
  mutate(offset = sigma_site*value) %>%
  mutate(a_site = a + offset) %>% 
  mutate(site = as.integer(parse_number(name)),
         model = "bayes") %>% 
  left_join(macro_fish_thin %>% ungroup %>% distinct(site_id, site_id_int, mat_s) %>% rename(site = site_id_int))

site_means = sim_site_means %>% 
  ggplot(aes(x = reorder(site_id, a_site), y = a_site)) +
  geom_violin(aes(group = mat_s)) + 
  geom_point(data = mle_mat %>% ungroup %>% distinct(site_id, b), aes(x = reorder(site_id, b), y = b),
             shape = 21, size = 1) + 
  coord_flip() + 
  theme_ggdist() + 
  labs(y = "b exponent",
       x = "NEON site",
       caption = "Dots are the old b estimates from the SizeSpectra package. Posteriors are from the new Stan model")

saveRDS(site_means, file = "plots/site_means.rds")

# plot varying intercept sd's
var_int_sds = mod_spectra %>% 
  as_draws_df() %>% 
  pivot_longer(cols = contains(c("raw_year","raw_site"))) %>% 
  mutate(year_site = str_sub(name, 11, 14)) %>% 
  group_by(.draw, year_site, sigma_year, sigma_site) %>% 
  summarize(mean_value = mean(value)) %>% 
  pivot_wider(names_from = year_site, values_from = mean_value) %>% 
  mutate(year_sd = sigma_year*year,
         site_sd = sigma_site*site) %>% 
  pivot_longer(cols = c(year_sd, site_sd))


var_int_sds %>% 
  ggplot(aes(x = value, y = name)) + 
  stat_halfeye()




# Plot x versus prob x<=X -------------------------------------------------

# 1) To make regression lines, simulate sequence between xmin and xmax for each sample,
# then join raw data identifiers and posteriors estimates of b (a_site, upper, lower, etc)...

# function to generate log sequence (other wise all of the numbers are too large for x to plot correctly)
# https://stackoverflow.com/questions/23901907/create-a-log-sequence-across-multiple-orders-of-magnitude
# logarithmic spaced sequence

lseq <- function(from=xmin, to=xmax, length.out=7) {
  exp(seq(log(from), log(to), length.out = length.out))
}

sim_grid <- macro_fish_thin %>%
  mutate(group = ID) %>% 
  ungroup() %>% 
  distinct(group, xmin, xmax) %>% 
  pivot_longer(cols = c(xmin, xmax)) %>%
  arrange(group, value) %>% 
  group_by(group) %>% 
  complete(value = lseq(min(value), max(value), length.out = 800)) %>% 
  select(-name) %>% 
  rename(x = value) %>% 
  left_join(macro_fish_thin %>% mutate(group = ID) %>% distinct(group, ID, site_id, year, xmin, xmax)) %>%  
  left_join(sim_regressions %>% 
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
  mutate(group = ID) %>% 
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
saveRDS(bayes_sim_tibble, file = "data/bayes_sim_tibble.rds")

# 7) Make plot

line_sim %>%
  filter(site_id == "ARIK") %>% 
  ggplot(aes(x = x, y = y_plb_med)) + 
  geom_line() +
  geom_ribbon(aes(ymin = y_plb_lower, ymax = y_plb_upper), alpha = 0.2) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~group, scales = "free") +
  geom_point(data = bayes_sim_tibble %>%  
               filter(group == 1:9)  , aes(x = dw, y = cumPropsim)) +
  NULL

