library(tidyverse)
library(sizeSpectra)

mod_spectra = readRDS("models/mod_spectra.rds")
macro_fish_mat <- readRDS("data/macro_fish_mat.rds") %>% 
  mutate(site_id_int = as.integer(as_factor(site_id))) %>%
  mutate(IDname = paste0(site_id, date),
         ID = as.integer(as.factor(IDname))) %>% 
  dplyr::mutate(Year = ID,
                Number = no_m2,
                bodyMass = dw) # different names that work with sizeSpectra package

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

# compute maximum likelihood estimates
# note: the code below produces the same estimates as the code in sizeSpectra
# this was confirmed by JSW on 10/21/2022.
b_exp = seq(-2.5, -0.5, length.out = 100)

negLL = macro_fish_mat %>% 
  expand_grid(b_exp = b_exp) %>% 
  mutate(negLL = -no_m2*(log((b_exp+1) / ( xmax^(b_exp+1) - xmin^(b_exp+1))) + b_exp*log(dw))) %>% 
  select(ID, negLL, b_exp) %>% 
    group_by(ID, b_exp) %>% 
    summarize(sumnegLL = sum(negLL))

min_negLL = negLL %>% 
  group_by(ID) %>% 
  filter(sumnegLL == min(sumnegLL))

min_negLL_sites = min_negLL %>% 
  left_join(macro_fish_mat %>% distinct(ID, site_id, date)) %>% 
  left_join(site_means_bayes %>% select(site_id, a_site))
  
# compare bayes and MLE
sim_site_means %>% 
  ggplot(aes(y = reorder(site_id, -a_site), x = a_site)) +
  geom_violin() +
  geom_point(data = min_negLL_sites, aes(x = b_exp), position = position_jitter(width = 0.01, height = 0),
             shape = 21) +
  labs(title = "Bayes model (violins) versus MLE (dots)")





# sizeSpectra code for MLE ------------------------------------------------------


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

