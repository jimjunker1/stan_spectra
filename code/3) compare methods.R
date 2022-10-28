library(tidyverse)
library(sizeSpectra)
library(brms)
library(tidybayes)
source("code/functions.R")

# load data
macro_fish_mat = readRDS(file = "data/macro_fish_mat.rds")
macro_fish_mat_globalminmax  = readRDS(file = "data/macro_fish_mat_globalminmax.rds")
macro_fish_mat_siteminmax = readRDS(file = "data/macro_fish_mat_siteminmax.rds") # this is really site max and global min

# load models
mod_spectra = readRDS("models/mod_spectra.rds")
mod_spectra_globalminmax = readRDS("models/sandbox/mod_spectra_globalminmax.rds")
mod_spectra_siteminmax = readRDS("models/sandbox/mod_spectra_siteminmax.rds")

# note: the mle code below produces the same estimates as the code in sizeSpectra
# this was confirmed by JSW on 10/21/2022.

# sample specific xmin/xmax -----------------------------------------------

sim_site_means <- mod_spectra %>% 
  as_draws_df() %>% 
  pivot_longer(cols = contains("raw_site")) %>% 
  mutate(offset = sigma_site*value) %>%
  mutate(a_site = a + offset) %>% 
  mutate(site = as.integer(parse_number(name)),
         model = "bayes") %>% 
  left_join(macro_fish_mat %>% ungroup %>% distinct(site_id, site_id_int, mat_s) %>% rename(site = site_id_int))


site_means_bayes = sim_site_means %>% 
  group_by(site, site_id, mat_s) %>% 
  median_qi(a_site)

# maximum likelihood

min_negLL = calc_exponent_maxlik(macro_fish_mat)

min_negLL_sites = min_negLL %>% 
  left_join(macro_fish_mat %>% distinct(ID, site_id, date)) %>% 
  left_join(site_means_bayes %>% select(site_id, a_site))
  
# compare bayes and MLE
sim_site_means %>% 
  ggplot(aes(y = reorder(site_id, -a_site), x = a_site)) +
  geom_violin() +
  geom_point(data = min_negLL_sites, aes(x = b_exp), position = position_jitter(width = 0.01, height = 0),
             shape = 21) +
  labs(title = "Bayes model (violins) versus MLE (dots)",
       subtitle = "Sample xmin and sample xmax")



# global xmin/xmax --------------------------------------------------------


sim_site_means_globalminmax <- mod_spectra_globalminmax %>% 
  as_draws_df() %>% 
  pivot_longer(cols = contains("raw_site")) %>% 
  mutate(offset = sigma_site*value) %>%
  mutate(a_site = a + offset) %>% 
  mutate(site = as.integer(parse_number(name)),
         model = "bayes") %>% 
  left_join(macro_fish_mat %>% ungroup %>% distinct(site_id, site_id_int, mat_s) %>% rename(site = site_id_int))


site_means_bayes_globalminmax = sim_site_means_globalminmax %>% 
  group_by(site, site_id, mat_s) %>% 
  median_qi(a_site)

min_negLL = calc_exponent_maxlik(macro_fish_mat_globalminmax)

min_negLL_globalminmax = min_negLL %>% 
  left_join(macro_fish_mat_globalminmax %>% distinct(ID, site_id, date)) %>% 
  left_join(site_means_bayes %>% select(site_id, a_site))

min_negLL_sites_globalminmax = min_negLL_globalminmax %>% 
  left_join(macro_fish_mat_globalminmax %>% distinct(ID, site_id, date)) %>% 
  left_join(site_means_bayes_globalminmax %>% select(site_id, a_site))

# compare bayes and MLE
sim_site_means_globalminmax %>% 
  ggplot(aes(y = reorder(site_id, -a_site), x = a_site)) +
  geom_violin() +
  geom_point(data = min_negLL_sites_globalminmax, aes(x = b_exp), position = position_jitter(width = 0.01, height = 0),
             shape = 21) +
  labs(title = "Bayes model (violins) versus MLE (dots)",
       subtitle = "Global xmin and global xmax")


# global xmin and site xmax -----------------------------------------------


sim_site_means_siteminmax <- mod_spectra_siteminmax %>% 
  as_draws_df() %>% 
  pivot_longer(cols = contains("raw_site")) %>% 
  mutate(offset = sigma_site*value) %>%
  mutate(a_site = a + offset) %>% 
  mutate(site = as.integer(parse_number(name)),
         model = "bayes") %>% 
  left_join(macro_fish_mat %>% ungroup %>% distinct(site_id, site_id_int, mat_s) %>% rename(site = site_id_int))


site_means_bayes_siteminmax = sim_site_means_siteminmax %>% 
  group_by(site, site_id, mat_s) %>% 
  median_qi(a_site)


min_negLL_siteminmax = calc_exponent_maxlik(macro_fish_mat_siteminmax)

min_negLL_sites_siteminmax = min_negLL_siteminmax %>% 
  left_join(macro_fish_mat_siteminmax %>% distinct(ID, site_id, date)) %>% 
  left_join(site_means_bayes_siteminmax %>% select(site_id, a_site))  # keep this global to compare with global model

# compare bayes and MLE
sim_site_means_siteminmax %>% 
  ggplot(aes(y = reorder(site_id, -a_site), x = a_site)) +
  geom_violin() +
  geom_point(data = min_negLL_sites_siteminmax, aes(x = b_exp), position = position_jitter(width = 0.01, height = 0),
             shape = 21) +
  labs(title = "Bayes model (violins) versus MLE (dots)",
       subtitle = "Global xmin and site xmax")





# compare xmin/xmax choices -------------------------------------------------------------
library(ggridges)

# site level predictions
all_posts = bind_rows(sim_site_means %>% mutate(model = "sample xmin/xmax"),
          sim_site_means_siteminmax %>% mutate(model = "site xmin/xmax"),
          sim_site_means_globalminmax %>% mutate(model = "global xmin/sample xmax")) 
  
all_posts %>% 
ggplot(aes(x = a_site, y = reorder(site_id, -a_site), fill = model)) + 
  geom_density_ridges(alpha = 0.5)

# regression slopes

site_slope = mod_spectra_siteminmax %>% 
  as_draws_df() %>% 
  select(a, beta, .draw) %>% 
  mutate(model = "site xmin/xmax")

global_slope = mod_spectra_globalminmax %>% 
  as_draws_df() %>% 
  select(a, beta, .draw) %>% 
  mutate(model = "global xmin/sample xmax")

sample_slope = mod_spectra %>% 
  as_draws_df() %>% 
  select(a, beta, .draw) %>% 
  mutate(model = "sample xmin/xmax")

all_slopes = bind_rows(site_slope, global_slope, sample_slope)
  
all_slopes %>% 
  ggplot(aes(x = model, y = beta)) + 
  stat_halfeye() +
  geom_hline(yintercept = 0)



# regressions
site_reg = mod_spectra_siteminmax %>% 
  as_draws_df() %>% 
  select(a, beta, .draw) %>% 
  expand_grid(mat_s = seq(-3, 3, length.out = 25)) %>% 
  mutate(b = a + beta*mat_s) %>% 
  group_by(mat_s) %>% 
  median_qi(b) %>% 
  mutate(model = "site xmin/xmax")

global_reg = mod_spectra_globalminmax %>% 
  as_draws_df() %>% 
  select(a, beta, .draw) %>% 
  expand_grid(mat_s = seq(-3, 3, length.out = 25)) %>% 
  mutate(b = a + beta*mat_s) %>% 
  group_by(mat_s) %>% 
  median_qi(b) %>% 
  mutate(model = "global xmin/sample xmax")

sample_reg = mod_spectra %>% 
  as_draws_df() %>% 
  select(a, beta, .draw) %>% 
  expand_grid(mat_s = seq(-3, 3, length.out = 25)) %>% 
  mutate(b = a + beta*mat_s) %>% 
  group_by(mat_s) %>% 
  median_qi(b) %>% 
  mutate(model = "sample xmin/xmax")

all_reg = bind_rows(sample_reg, global_reg, site_reg)


all_reg %>% 
  ggplot(aes(x = mat_s, y = b, fill = model, color = model)) + 
  geom_line() + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2) +
  ylim(-2, -1) +
  facet_wrap(~model, nrow = 1)



# raw sizeSpectra code for MLE ------------------------------------------------------
id = 13

dat = macro_fish_mat 

test = dat %>% filter(ID == id) %>% 
  rename(bodyMass = dw, 
         Number = no_m2)


valCounts = test %>% 
  select(bodyMass, Number) %>% 
  ungroup

MLE.valCounts = valCounts %>% 
  select(bodyMass, Number) %>% 
  group_by(bodyMass) %>% 
  summarize(Count = sum(Number)) %>% 
  arrange(desc(bodyMass))

x = test$bodyMass
c = test$Number
sum.log.x = sum(log(x))
xmin = min(x)
xmax = max(x)
sumCounts = sum(MLE.valCounts$Count)
MLE.K = dim(MLE.valCounts)[1]
MLE.sumCntLogMass = sum(MLE.valCounts$Count * log(MLE.valCounts$bodyMass))
PL.bMLE.counts = 1/( log(xmin) - MLE.sumCntLogMass/sumCounts) - 1

negLL.PLB.counts = function(b, x, c, K=length(c), xmin=min(x), xmax=max(x),
                            sumclogx = sum(c * log(x)))
{
  if(xmin <= 0 | xmin >= xmax | length(x) != K | length(c) != K |
     c[1] == 0 | c[K] == 0 | min(c) < 0)
    stop("Parameters out of bounds in negLL.PLB.counts")
  n = sum(c)
  if(b != -1)
  { neglogLL = -n * log( ( b + 1) / (xmax^(b + 1) - xmin^(b + 1)) ) -
    b * sumclogx
  } else
  { neglogLL = n * log( log(xmax) - log(xmin) ) + sumclogx
  }
  return(neglogLL)
}

PLB.minLL =  nlm(negLL.PLB.counts, p=PL.bMLE.counts,
                 x=MLE.valCounts$bodyMass, c=MLE.valCounts$Count, K=MLE.K,
                 xmin=xmin, xmax=xmax, sumclogx=MLE.sumCntLogMass)

# compare (these should be the same)
PLB.minLL$estimate   # sizeSpectra

min_negLL %>% filter(ID == id) %>% select(b_exp) # code by hand

