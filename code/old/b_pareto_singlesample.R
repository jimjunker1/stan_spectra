library(tidyverse)
library(brms)
library(rstan)
library(tidybayes)
library(janitor)
library(sizeSpectra)
library(ggridges)

N = 1000
b = -1.4
xmin = 0.001
xmax = 10000
sim_nocounts = tibble(x = rPLB(n = 1000, b = b, xmin = xmin, xmax = xmax))

sim_data = tibble(x = sim_nocounts$x,
                  xmin = xmin,
                  xmax = xmax)


test_nocounts = stan(file = "models/b_pareto_singlesample.stan",
                     data = list(N = nrow(sim_data),
                                 x = sim_data$x,
                                 xmin = sim_data$xmin,
                                 xmax = sim_data$xmax),
                     iter = 500, chains = 1)

nocounts_post <- test_nocounts %>% 
  as_draws_df()

nocounts_posts_summary <- nocounts_post %>% 
  median_qi(b_exp)

nocounts_post_predict <- nocounts_post %>% select(-lp__) %>% 
  sample_n(10) %>% 
  expand_grid(y_rep = 1:1000) %>% 
  mutate(x = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1)),
         sim = as.integer(as.factor(.draw))) %>% 
  mutate(model = "ypred")


post_raw_compare <- nocounts_post_predict %>% 
  bind_rows(sim_data %>% mutate(model = "y_raw",
                                sim = 0)) 

post_raw_compare %>% 
  ggplot(aes(x = sim, y = x + 0.1, color = model)) + 
  geom_point(aes(group = sim), position = position_jitter(width = 0.2, height  =0)) + 
  scale_y_log10()
### Simulate regression lines

# 1) To make regression lines, simulate sequence between xmin and xmax for each sample,
# then join raw data identifiers and posteriors estimates of b (a_site, upper, lower, etc)...

# function to generate log sequence (other wise all of the numbers are too large for x to plot correctly)
# https://stackoverflow.com/questions/23901907/create-a-log-sequence-across-multiple-orders-of-magnitude
# logarithmic spaced sequence

lseq <- function(from=xmin, to=xmax, length.out=7) {
  exp(seq(log(from), log(to), length.out = length.out))
}

sim_grid_single <- sim_data %>%
  distinct(xmin, xmax) %>% 
  pivot_longer(cols = c(xmin, xmax)) %>%
  arrange(value) %>% 
  # group_by(group) %>% 
  # complete(value = seq(min(value), max(value), length.out = 100000)) %>% 
  complete(value = lseq(min(value), max(value), length.out = 1000)) %>% 
  select(-name) %>% 
  rename(y_raw = value) %>% 
  mutate(xmin = unique(sim_data$xmin), 
         xmax = unique(sim_data$xmax),
         a_site = nocounts_posts_summary$b_exp,
         .lower = nocounts_posts_summary$.lower,
         .upper = nocounts_posts_summary$.upper) 


# 2) ...then simulate prob x >= x for each x. y_plb estimates come from line 155 here: https://github.com/andrew-edwards/fitting-size-spectra/blob/master/code/PLBfunctions.r


line_sim_single <- sim_grid_single %>% 
  mutate(x = y_raw) %>% 
  mutate(y_plb_med = 1 - (x^(a_site+1) - xmin^(a_site+1))/(xmax^(a_site+1) - xmin^(a_site+1)), # simulate prob x>=x
         y_plb_lower = 1 - (x^(.lower+1) - xmin^(.lower+1))/(xmax^(.lower+1) - xmin^(.lower+1)),
         y_plb_upper = 1 - (x^(.upper+1) - xmin^(.upper+1))/(xmax^(.upper+1) - xmin^(.upper+1))) %>%
  filter(y_plb_med > 0) %>% 
  filter(y_plb_lower > 0) %>% 
  filter(y_plb_upper > 0) %>% 
  arrange(x) %>%
  mutate(x = round(x, 5)) %>% 
  distinct(x, .keep_all = T)

### Simulate raw data (i.e., data that "would" have been collected after accounting for no_m2)

# 3) To simulate raw data, get counts of organisms and cumulative counts...
dat_bayes_counts_single = sim_data %>% 
  # mutate(group = paste(site_id, year_month, sep = "_")) %>% 
  select(x) %>% ungroup() %>% 
  group_by(x) %>%
  mutate(no_m2 = 1) %>% 
  summarize(Count = sum(no_m2)) %>% 
  arrange(desc(x)) %>% 
  # group_by(group) %>% 
  mutate(cumSum = cumsum(Count),
         cumProp = cumSum / sum(Count),
         length = ceiling(sum(Count))) 

# 4) then generate sequence of values to simulate over
dat_bayes_sim_single <- dat_bayes_counts_single %>% 
  dplyr::group_by(length) %>% 
  dplyr::summarize(min_cumProp = min(cumProp)) %>% 
  # dplyr::group_by(group) %>% 
  dplyr::do(dplyr::tibble(cumPropsim = seq(.$min_cumProp, 1, length = .$length/10))) # dividing by something reduces file size by limiting iterations, but check for accuracy

# 5) then simulate cumulative proportion data to plot against MLE estimates by group
# make lists first
dat_bayes_simlist <- list(dat_bayes_sim_single)
dat_bayes_countslist <- list(dat_bayes_counts_single)
bayes_sim = list() # empty list to population

# simulate data with for loop
for(i in 1:length(dat_bayes_simlist)){
  bayes_sim[[i]] = dat_bayes_simlist[[i]] %>% dplyr::as_tibble() %>% 
    dplyr::mutate(dw = dat_bayes_countslist[[i]][findInterval(dat_bayes_simlist[[i]]$cumPropsim,
                                                              dat_bayes_countslist[[i]]$cumProp), ]$x)
}

# 6) Create data frame with "raw" data to plot
bayes_sim_tibble <- dplyr::bind_rows(bayes_sim) # dots to plot...very large file


# 7) Make plot

line_sim_single %>%  
  ggplot(aes(x = x, y = y_plb_med)) + 
  geom_line() +
  geom_ribbon(aes(ymin = y_plb_lower, ymax = y_plb_upper), alpha = 0.2) +
  scale_x_log10() +
  scale_y_log10() +
  # facet_wrap(~group, scales = "free") +
  geom_point(data = bayes_sim_tibble , aes(x = dw, y = cumPropsim)) +
  NULL
