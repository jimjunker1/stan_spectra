library(tidyverse)
library(sizeSpectra)

n_sim = 100000
xmax = 10000
xmin = 0.001
b = -1.5
u = runif(n_sim)

sim_hand = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))
sim_plb = rPLB(n = n_sim, b = b, xmin =xmin, xmax =xmax)

tibble(sim_hand = sim_hand,
       sim_plb = sim_plb) %>% 
  pivot_longer(cols = everything()) %>% 
  ggplot(aes(x = value, fill = name)) + 
  geom_density(alpha = 0.5) + 
  scale_x_log10()





