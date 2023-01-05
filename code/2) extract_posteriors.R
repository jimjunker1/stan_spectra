library(tidyverse)
library(brms)
library(rstan)
library(tidybayes)
library(janitor)
library(sizeSpectra)
library(ggridges)


# lambda models -----------------------------------------------------------

# model
mod_spectra = readRDS(file = "models/fit_interaction.rds")

# data
macro_fish_mat_siteminmax = readRDS(file = "data/macro_fish_mat_siteminmax.rds") 

# posts for ids
id_posts <- mod_spectra %>% 
  as_draws_df() %>% 
  pivot_longer(cols = contains("raw_site")) %>% 
  mutate(ID = as.integer(parse_number(name)),
         model = "bayes") %>% 
  select(name, value, ID, sigma_site, a, contains("beta"), .draw) %>% 
  left_join(macro_fish_mat_siteminmax %>% 
              ungroup %>% 
              distinct(ID, site_id, site_id_int, mat_s, log_gpp_s)) %>%
  mutate(offset = sigma_site*value) %>%
  mutate(lambda = a + beta_mat*mat_s + beta_gpp*log_gpp_s + beta_gpp_mat*log_gpp_s*mat_s + offset) 

id_summaries = id_posts %>% 
  group_by(ID, site_id, mat_s, log_gpp_s) %>% 
  median_qi(lambda) 

saveRDS(id_summaries, file = "posteriors/id_summaries.rds")

# posts for regression
# temperature on the x axis with 3 gpp levels
gpp_conds = macro_fish_mat_siteminmax %>% ungroup %>% 
  distinct(log_gpp_s) %>% 
  summarize(log_gpp_s = quantile(log_gpp_s, c(0.1, 0.5, 0.9)))

mat_conds = seq(min(macro_fish_mat_siteminmax$mat_s), max(macro_fish_mat_siteminmax$mat_s), 
                length.out = 15)

temp_x_regression = as_draws_df(mod_spectra) %>% select(contains("beta"), a) %>% 
  expand_grid(mat_s = mat_conds) %>% 
  expand_grid(gpp_conds) %>% 
  mutate(lambda =  a + beta_mat*mat_s + beta_gpp*log_gpp_s + beta_gpp_mat*log_gpp_s*mat_s)

temp_x_summaries = temp_x_regression %>% 
  group_by(mat_s, log_gpp_s) %>% 
  median_qi(lambda)

saveRDS(temp_x_summaries, file = "posteriors/temp_x_summaries.rds")


# gpp on the x axis with 3 temp levels

mat_conds = macro_fish_mat_siteminmax %>% ungroup %>% 
  distinct(mat_s) %>% 
  summarize(mat_s = quantile(mat_s, c(0.1, 0.5, 0.9)))

gpp_conds = seq(min(macro_fish_mat_siteminmax$log_gpp_s), max(macro_fish_mat_siteminmax$log_gpp_s), 
                length.out = 15)

gpp_x_regression = as_draws_df(mod_spectra) %>% select(contains("beta"), a) %>% 
  expand_grid(log_gpp_s = gpp_conds) %>% 
  expand_grid(mat_conds) %>% 
  mutate(lambda =  a + beta_mat*mat_s + beta_gpp*log_gpp_s + beta_gpp_mat*log_gpp_s*mat_s)

gpp_x_summaries = gpp_x_regression %>% 
  group_by(mat_s, log_gpp_s) %>% 
  median_qi(lambda)


saveRDS(gpp_x_summaries, file = "posteriors/gpp_x_summaries.rds")


# biomass models ----------------------------------------------------------

gpp_conds = macro_fish_mat_siteminmax %>% ungroup %>% 
  distinct(log_gpp_s) %>% 
  summarize(log_gpp_s = quantile(log_gpp_s, c(0.1, 0.5, 0.9)))

mat_conds = tibble(mat_s = seq(min(macro_fish_mat_siteminmax$mat_s), max(macro_fish_mat_siteminmax$mat_s), 
                length.out = 15))

mod_biomass = readRDS(file = "models/mod_biomass.rds")

library(tidybayes)
biomass_posts = mat_conds %>% 
  expand_grid(gpp_conds) %>% 
  add_epred_draws(mod_biomass, re_formula = NA)

saveRDS(biomass_posts, file = "posteriors/biomass_posts.rds")



# other -------------------------------------------------------------------



sim_regressions %>% 
  ggplot(aes(x = mat_s, y = a_site)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2) +
  theme_ggdist() + 
  geom_point(data = mle_mat %>% distinct(b, mat_s), aes(y = b)) + 
  labs(y = "b exponent",
       x = "Standardized Temperature",
       caption = "Dots are the old b estimates from the SizeSpectra package. Regression lines are from the new Stan model")



sim_site_means <- mod_spectra %>% 
  as_draws_df() %>% 
  pivot_longer(cols = contains("raw_site")) %>% 
  mutate(offset = sigma_site*value) %>%
  mutate(a_site = a + offset) %>% 
  mutate(site = as.integer(parse_number(name)),
         model = "bayes") %>% 
  left_join(macro_fish_mat %>% ungroup %>% distinct(site_id, site_id_int, mat_s) %>% rename(site = site_id_int))

sim_site_means %>% 
  ggplot(aes(x = reorder(site_id, a_site), y = a_site)) +
  geom_violin(aes(group = mat_s)) + 
  # geom_point(data = mle_mat %>% ungroup %>% distinct(site_id, b), aes(x = reorder(site_id, b), y = b),
  #            shape = 21, size = 1) + 
  coord_flip() + 
  theme_ggdist() + 
  labs(y = "b exponent",
       x = "NEON site",
       caption = "Dots are the old b estimates from the SizeSpectra package. Posteriors are from the new Stan model")


mod_regression = mod_spectra_siteminmax %>% 
  as_draws_df() %>% 
  select(a, beta, .draw) %>% 
  expand_grid(mat_s = seq(-3, 3), length.out = 100) %>% 
  mutate(b = a + beta*mat_s) %>% 
  group_by(mat_s) %>% 
  median_qi(b)


mod_regression %>% 
  ggplot(aes(x = mat_s, y = b)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2) +
  ylim(-1.6, -1.3)


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

sim_grid <- macro_fish_mat %>%
  mutate(group = ID) %>% 
  ungroup() %>% 
  distinct(group, xmin, xmax) %>% 
  pivot_longer(cols = c(xmin, xmax)) %>%
  arrange(group, value) %>% 
  group_by(group) %>% 
  complete(value = lseq(min(value), max(value), length.out = 800)) %>% 
  select(-name) %>% 
  rename(x = value) %>% 
  left_join(macro_fish_mat %>% mutate(group = ID) %>% distinct(group, ID, site_id, year, xmin, xmax)) %>%  
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
dat_bayes_counts = macro_fish_mat %>% 
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




# Posterior Predictive ----------------------------------------------------

# 1) get posterior mean and sd for chosen parameters
post_mean_sd_parameter = as_draws_df(mod_spectra) %>% as_tibble() %>% clean_names() %>% 
  select(a) %>%  # change for different parameters/groups
  summarize(b = mean(a),
            sd = sd(a)) %>% 
  expand_grid(site_id = macro_fish_mat %>% ungroup %>% distinct(site_id) %>% pull())

# 2) simulate y_pred
sim_ypred <- macro_fish_mat %>% 
  left_join(post_mean_sd_parameter) %>%       # add posterior mean/sd/upper/lower
  expand_grid(sim = 1:10) %>%                 # number of data sets to simulate
  mutate(u = runif(nrow(.))) %>%              # uniform sample for simulation
  mutate(y_pred = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))) %>% # simulate y_pred (via Edwards github for rPLB) - confirmed in code/rplb_by_hand.R
  group_by(ID) %>% 
  mutate(rank = rank(y_pred)) %>%             # rank within groups
  mutate(group = paste(site_id, year, sep = "_")) %>%  # clean up to combine with raw data
  mutate(y_pred = round(y_pred, 5)) %>%
  select(ID,site_id, y_pred, sim) %>%
  group_by(ID, y_pred, site_id, sim) %>%
  count(y_pred) %>%
  mutate(model = "y_pred") %>% 
  rename(dw = y_pred,
         no_m2 = n)

# 3) get counts of organisms and cumulative counts...
sim_bayes_counts = sim_ypred %>% 
  mutate(group = paste(ID, sim, sep = "_")) %>% 
  select(dw, no_m2, group) %>% ungroup() %>% 
  group_by(dw, group) %>% 
  summarize(Count = sum(no_m2)) %>% 
  arrange(group, desc(dw)) %>% 
  group_by(group) %>% 
  mutate(cumSum = cumsum(Count),
         cumProp = cumSum / sum(Count),
         length = ceiling(sum(Count))) 

# 4) then generate sequence of values to simulate over
sim_bayes_sim <- sim_bayes_counts %>% 
  dplyr::group_by(group, length) %>% 
  dplyr::summarize(min_cumProp = min(cumProp)) %>% 
  dplyr::group_by(group) %>% 
  dplyr::do(dplyr::tibble(cumPropsim = seq(.$min_cumProp, 1, length = .$length/10))) # dividing by something reduces file size by limiting iterations, but check for accuracy

# 5) then simulate cumulative proportion data to plot against MLE estimates by group
# make lists first
sim_bayes_simlist <- sim_bayes_sim %>% dplyr::group_by(group) %>% dplyr::group_split() 
sim_bayes_countslist <- sim_bayes_counts %>% dplyr::group_by(group) %>% dplyr::group_split() 
bayes_sim_ypred = list() # empty list to population

# simulate data with for loop
for(i in 1:length(sim_bayes_simlist)){
  bayes_sim_ypred[[i]] = sim_bayes_simlist[[i]] %>% dplyr::as_tibble() %>% 
    dplyr::mutate(dw = sim_bayes_countslist[[i]][findInterval(sim_bayes_simlist[[i]]$cumPropsim,
                                                              sim_bayes_countslist[[i]]$cumProp), ]$dw)
}

# 6) Create data frame with simulated data to plot
bayes_sim_ypred_tibble <- dplyr::bind_rows(bayes_sim_ypred) %>% # dots to plot...very large file
  mutate(model = "y_pred")



# 7) combine with raw sims
y_pred_raw = bind_rows(bayes_sim_ypred_tibble %>% separate(group, c("group","sim", sep = "_")) %>% mutate(group = as.numeric(group),
                                                                                                          sim = as.numeric(sim)), 
                       bayes_sim_tibble %>% mutate(model = "y_raw", sim = -1))

# sim_raw = bayes_sim_tibble %>%
#   mutate(model = "y_bayes_sim",
#          sim = -1,
#          site_id = str_sub(group, 1, 4)) %>%
#   mutate(ID = as.integer(as.factor(group))) %>%
#   select(-site_id) %>%
#   right_join(macro_fish_mat %>% ungroup %>% distinct(site_id, ID))



# combine y_raw and y_pred
# y_pred_raw = bind_rows(sim_ypred, sim_raw)


# plot
y_pred_raw %>% 
  # filter(dw != 0) %>% 
  # mutate(ID = as.factor(ID),
  #        sim = as.factor(sim)) %>% 
  ggplot(aes(y = dw, fill = model, x = sim)) + 
  geom_boxplot(aes(group = interaction(model, sim)),
               outlier.shape = NA) + 
  # geom_jitter(width = 0.1, height = 0, size = 0.1) +
  # facet_wrap(~site_id) +
  scale_y_log10() +
  NULL


