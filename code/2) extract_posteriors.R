library(tidyverse)
library(brms)
library(rstan)


# lambda models -----------------------------------------------------------

# models
fit_interaction = readRDS(file = "models/fit_interaction_sites.rds")
fit_temponly = readRDS(file = "models/fit_temponly.rds")
fit_interaction_sites = readRDS(file = "models/fit_interaction_sites.rds")

# data
dat = readRDS(file = "data/macro_fish_mat_siteminmax.rds") 

# extract posts
posts_interaction = as_draws_df(fit_interaction)
posts_temponly = as_draws_df(fit_temponly)

posts_interaction %>% select(contains("sigma")) %>% 
  pivot_longer(cols = everything()) %>% 
  group_by(name) %>% 
  median_qi(value)

posts_pivoted = posts_interaction %>% 
  pivot_longer(cols = contains("alpha_raw_sample"),
               names_to = "sample_id", values_to = "sample_offset") %>% 
  pivot_longer(cols = contains("alpha_raw_site"),
               names_to = "site_id_int", values_to = "site_offset") %>% 
  pivot_longer(cols = contains("alpha_raw_year"),
               names_to = "year_id", values_to = "year_offset")

posts_data = posts_pivoted %>% 
  filter(.draw <= 500) %>% 
  mutate(sample_id_int = as.integer(parse_number(sample_id)),
         site_id_int = as.integer(parse_number(site_id_int)),
         year_id = as.integer(parse_number(year_id)))  %>% 
  right_join(dat, by = c("sample_id_int", "site_id_int", "year_id"))


line_posts = posts_interaction %>% 
  filter(.draw <= 1000) %>% 
  expand_grid(dat %>% distinct(log_gpp_s)) %>% 
  expand_grid(mat_s = dat %>% distinct(mat_s) %>% median_qi(mat_s) %>% 
                select(-.width, -.point, -.interval) %>%
                pivot_longer(cols = everything()) %>% pull()) %>% 
  mutate(lambda = a + beta_mat*mat_s + beta_gpp*log_gpp_s + beta_gpp_mat*log_gpp_s*mat_s) %>% 
  group_by(log_gpp_s, mat_s) %>% 
  median_qi(lambda)

sample_posts = posts_data %>% 
  mutate(lambda = a + beta_mat*mat_s + beta_gpp*log_gpp_s + beta_gpp_mat*log_gpp_s*mat_s +
           year_offset*sigma_year + sample_offset*sigma_sample + site_offset*sigma_site) %>% 
  group_by(site_id, year, sample_id, log_gpp_s) %>%
  median_qi(lambda)


sample_posts %>% 
  ggplot(aes(x = log_gpp_s, y = lambda)) + 
  geom_pointrange(aes(ymin = .lower, ymax = .upper), 
                  position = position_jitter(width = 0.1), size = 0.1) + 
  geom_line(data = line_posts, aes(group = as.factor(log_gpp_s))) + 
  geom_ribbon(data = line_posts, aes(ymin = .lower, ymax = .upper), alpha = 0.2) +
  facet_wrap(~mat_s)

fit_interaction

# posts for ids

wrangle_posts = function(model, name = NA, data = NA, group1 = "alpha_raw_site", group2 = "alpha_raw_year"){
  model %>% 
    as_draws_df() %>% 
    pivot_longer(cols = contains(group1),
                 names_to = "sample_id_int", values_to = "sample_offset") %>% 
    pivot_longer(cols = contains(group2),
                 names_to = "year_id", values_to = "year_offset") %>% 
    mutate(sample_id_int = as.integer(parse_number(sample_id_int)),
           year_id = as.integer(as.factor(parse_number(year_id))),
           model = name) %>% 
    select(sample_id_int, year_id, sample_offset, year_offset, sigma_site, sigma_year, a, contains("beta"), .draw, model) %>% 
    right_join(data %>% 
                 ungroup %>% 
                 distinct(site_id, year_id, year, sample_id_int, mat_s, log_gpp_s))  
}


interaction_posts_wrangled = wrangle_posts(model = fit_interaction, data = dat, name = "interaction") %>% 
  mutate(sample_offset = sigma_site*sample_offset,
         year_offset = sigma_year*year_offset) %>%
  mutate(lambda = a + beta_mat*mat_s + beta_gpp*log_gpp_s + beta_gpp_mat*log_gpp_s*mat_s + 
           sample_offset +
           year_offset) %>% 
  ungroup

temp_posts_wrangled = wrangle_posts(model = fit_temponly, data = dat, name = "temponly") %>% 
  mutate(sample_offset = sigma_site*sample_offset,
         year_offset = sigma_year*year_offset) %>%
  mutate(lambda = a + beta_mat*mat_s + 
           sample_offset +
           year_offset) %>% 
  ungroup

interactionsite_posts_wrangled = as_draws_df(fit_interaction_sites) %>% 
  pivot_longer(cols = contains("raw_year"),
               names_to = "year_id", values_to = "year_offset") %>% 
  pivot_longer(cols = contains("raw_sample"),
               names_to = "sample_id_int", values_to = "sample_offset") %>% 
  mutate(sample_id_int = as.integer(parse_number(sample_id_int)),
         year_id = as.integer(parse_number(year_id))) %>% 
  right_join(dat %>% ungroup %>% distinct(year_id, sample_id_int)) %>% 
  pivot_longer(cols = contains("raw_site"),
               names_to = "site_id_int", values_to = "site_offset")  %>% 
  mutate(site_id_int = as.integer(parse_number(site_id_int))) %>% 
  right_join(dat %>% ungroup %>% distinct(year_id, sample_id_int, site_id_int), by = c("year_id", "sample_id_int", "site_id_int")) 
  


saveRDS(interaction_posts_wrangled, file = "posteriors/interaction_posts_wrangled.rds")
saveRDS(temp_posts_wrangled, file = "posteriors/temp_posts_wrangled.rds")
saveRDS(interactionsite_posts_wrangled, file = "posteriors/interactionsite_posts_wrangled.rds")

id_summaries = id_posts %>% 
  group_by(ID, site_id, mat_s, log_gpp_s) %>% 
  median_qi(lambda) 

saveRDS(id_summaries, file = "posteriors/id_summaries.rds")

# posts for regression
# temperature on the x axis with 3 gpp levels
gpp_conds = dat %>% ungroup %>% 
  distinct(log_gpp_s) %>% 
  summarize(log_gpp_s = quantile(log_gpp_s, c(0.1, 0.5, 0.9)))

mat_conds = seq(min(dat$mat_s), max(dat$mat_s), 
                length.out = 15)

temp_x_regression = as_draws_df(fit_interaction) %>% select(contains("beta"), a) %>% 
  expand_grid(mat_s = mat_conds) %>% 
  expand_grid(gpp_conds) %>% 
  mutate(lambda =  a + beta_mat*mat_s + beta_gpp*log_gpp_s + beta_gpp_mat*log_gpp_s*mat_s)

temp_x_summaries = temp_x_regression %>% 
  group_by(mat_s, log_gpp_s) %>% 
  median_qi(lambda)

saveRDS(temp_x_summaries, file = "posteriors/temp_x_summaries.rds")


# gpp on the x axis with 3 temp levels

mat_conds = dat %>% ungroup %>% 
  distinct(mat_s) %>% 
  summarize(mat_s = quantile(mat_s, c(0.1, 0.5, 0.9)))

gpp_conds = seq(min(dat$log_gpp_s), max(dat$log_gpp_s), 
                length.out = 15)

gpp_x_regression = as_draws_df(fit_interaction) %>% select(contains("beta"), a) %>% 
  expand_grid(log_gpp_s = gpp_conds) %>% 
  expand_grid(mat_conds) %>% 
  mutate(lambda =  a + beta_mat*mat_s + beta_gpp*log_gpp_s + beta_gpp_mat*log_gpp_s*mat_s)

gpp_x_summaries = gpp_x_regression %>% 
  group_by(mat_s, log_gpp_s) %>% 
  median_qi(lambda)


saveRDS(gpp_x_summaries, file = "posteriors/gpp_x_summaries.rds")


# biomass models ----------------------------------------------------------

gpp_conds = dat %>% ungroup %>% 
  distinct(log_gpp_s) %>% 
  summarize(log_gpp_s = quantile(log_gpp_s, c(0.1, 0.5, 0.9)))

mat_conds = tibble(mat_s = seq(min(dat$mat_s), max(dat$mat_s), 
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



sim_site_means <- fit_interaction %>% 
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



