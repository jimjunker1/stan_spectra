library(tidyverse)
library(sizeSpectra)

# Compute log-likelihood by hand. Conform that it recaptures the assigned b
b = -1.5
xmin = 0.001
xmax = 1000
N = 2000

# simulate data from known parameters
sims2 <- tibble(x = rPLB(n = N, b = b, xmin = xmin, xmax = xmax)) %>% 
  mutate(x = round(x, 3)) %>% 
  count(x, name = "counts") %>% 
  mutate(xmin = xmin, 
         xmax = xmax) %>% 
  expand_grid(b_exp = seq(-8, 4, length.out = 30)) %>% 
  group_by(b_exp) %>% 
  # compute log likelihood via Edwards 2020 eqn2, modified for counts via Supplemental of that papers
  mutate(log_lik = counts*(log((b_exp+1) / ( xmax^(b_exp+1) - xmin^(b_exp+1))) + b_exp*log(x))) %>% 
  group_by(b_exp) %>% 
  summarize(negsum_log_lik = -sum(log_lik)) # compute negative log-likelihood

# check
sims2 %>% 
  ggplot(aes(x = b_exp, y = negsum_log_lik)) + 
  geom_point() + 
  geom_vline(xintercept = b)


