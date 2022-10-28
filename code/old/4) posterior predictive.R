library(rstan)
library(tidyverse)
library(janitor)

macro_fish_thin <- readRDS("data/macro_fish_thin.rds") 
mod_spectra <- readRDS("models/mod_spectra.rds")
bayes_sim_tibble = readRDS("data/bayes_sim_tibble.rds")


# posterior predictive using example from here:https://stackoverflow.com/questions/66081593/how-to-generate-data-from-a-distribution-whose-cdf-is-not-in-closed-form
# this produces the empirical cdf

# 1) get pdf

#
#   Data generation for b not equal to -1
#
#    xmax is the largest value in the population
#    xmin is the smallest value in the population
#

post_pred = macro_fish_thin %>% 
  ungroup() %>% 
  mutate(b = rnorm(1, -1.3, 0.001),
         u = runif(nrow(.), min = 0, max = 1),
         y_pred = ( u*( xmax^(b+1) - xmin^(b+1) ) + xmin^(b+1) )^(1/(b+1))) %>% 
  select(dw, no_m2, y_pred) %>% 
  arrange(no_m2) %>% 
  mutate(counts = as.integer(no_m2/min(no_m2))) 



post_pred %>% 
  mutate(y_pred = y_pred/no_m2) %>% 
  pivot_longer(cols = c( -counts,-no_m2)) %>%
  ggplot(aes(x = name, y = value)) +
  geom_boxplot(aes(group = name)) +
  scale_y_log10() 



post_pred %>% 
  # pivot_longer(cols = -no_m2) %>% 
  ggplot(aes(x = dw, y = y_pred/no_m2)) +
  geom_point() + 
  scale_x_log10() + 
  scale_y_log10()


n = 1
xmin = 0.1
xmax = 10000
b = -1.5
counts = 5


u =runif(n, min = 0, max =1)

#  Step 2

#  if b not equal to -1

x = ( u*( xmax^(b+1) - xmin^(b+1) ) + xmin^(b+1) )^(1/(b+1))

rep(x, counts)


#  if b equal to -1

x = exp( u * ( log(xmax) - log(xmin) )+ log(xmin) )




# Posterior Predictive ----------------------------------------------------

# 1) get posterior mean and sd for chosen parameters
post_mean_sd_parameter = as_draws_df(mod_spectra) %>% as_tibble() %>% clean_names() %>% 
  select(contains) %>%  # change for different parameters/groups
  summarize(b = mean(a),
            sd = sd(a)) %>% 
  expand_grid(site_id = macro_fish_thin %>% ungroup %>% distinct(site_id) %>% pull())

# 2) simulate y_pred
sim_ypred <- macro_fish_thin %>% 
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

# # 3) get counts of organisms and cumulative counts...
# sim_bayes_counts = sim_ypred %>% 
#   mutate(group = paste(ID, sim, sep = "_")) %>% 
#   select(dw, no_m2, group) %>% ungroup() %>% 
#   group_by(dw, group) %>% 
#   summarize(Count = sum(no_m2)) %>% 
#   arrange(group, desc(dw)) %>% 
#   group_by(group) %>% 
#   mutate(cumSum = cumsum(Count),
#          cumProp = cumSum / sum(Count),
#          length = ceiling(sum(Count))) 
# 
# # 4) then generate sequence of values to simulate over
# sim_bayes_sim <- sim_bayes_counts %>% 
#   dplyr::group_by(group, length) %>% 
#   dplyr::summarize(min_cumProp = min(cumProp)) %>% 
#   dplyr::group_by(group) %>% 
#   dplyr::do(dplyr::tibble(cumPropsim = seq(.$min_cumProp, 1, length = .$length/10))) # dividing by something reduces file size by limiting iterations, but check for accuracy
# 
# # 5) then simulate cumulative proportion data to plot against MLE estimates by group
# # make lists first
# sim_bayes_simlist <- sim_bayes_sim %>% dplyr::group_by(group) %>% dplyr::group_split() 
# sim_bayes_countslist <- sim_bayes_counts %>% dplyr::group_by(group) %>% dplyr::group_split() 
# bayes_sim_ypred = list() # empty list to population
# 
# # simulate data with for loop
# for(i in 1:length(sim_bayes_simlist)){
#   bayes_sim_ypred[[i]] = sim_bayes_simlist[[i]] %>% dplyr::as_tibble() %>% 
#     dplyr::mutate(dw = sim_bayes_countslist[[i]][findInterval(sim_bayes_simlist[[i]]$cumPropsim,
#                                                               sim_bayes_countslist[[i]]$cumProp), ]$dw)
# }
# 
# # 6) Create data frame with simulated data to plot
# bayes_sim_ypred_tibble <- dplyr::bind_rows(bayes_sim_ypred) %>% # dots to plot...very large file
#   mutate(model = "y_pred") %>% 
#   separate(group, c("group","sim", sep = "_")) %>% 
#   mutate(group = as.numeric(group),
#          sim = as.numeric(sim))


# 7) combine with raw sims
y_pred_raw = bind_rows(sim_ypred,
                       macro_fish_thin %>% mutate(model = "y_raw",
                                                  sim = 0) %>% 
                         mutate(dw = dw*no_m2))



# plot
y_pred_raw %>% 
  ggplot(aes(y = dw, fill = model, x = sim)) + 
  geom_violin(aes(group = interaction(model, sim)),
              outlier.shape = NA) + 
  scale_y_log10() +
  NULL


y_pred_raw %>% 
  filter(group == 1:10) %>% 
  ggplot(aes(y = dw, fill = model, x = sim)) + 
  geom_violin(aes(group = interaction(model, sim))) + 
  scale_y_log10() +
  facet_wrap(~group, scales = "free_y") +
  NULL



