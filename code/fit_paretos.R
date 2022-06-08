library(tidyverse)
library(brms)
library(rstan)
library(tidybayes)
library(sizeSpectra)


# sim data
N_ind = 2000
xmin = 0.01
xmax = 1000
b = -2.4
x_ind = rPLB(n = N_ind, b = b, xmin = xmin, xmax = xmax) # simulate individual masses from bounded power law


# Fit truncated Pareto with counts (same as eightMethods.count MLE, but as Bayes) ----------------------------------------

# simulate counts by rounding x's and counting the number of individuals in each rounded "bin"
fake <- tibble(x = round(x_ind,4)) %>% count(x, name = "c") # assume

# make vectors for Stan
x = fake$x # 
c = fake$c
N = length(c)

# fit stan model
fit_counts <- stan(file = "models/stan_pareto_truncated_counts.stan",
                   data = list(N = N, 
                               x = x,
                               c = c,
                               xmin = xmin,
                               xmax = xmax))

# Fit truncated Pareto with individuals (same as eightMethodsMEE MLE, but as Bayes) ----------------------------------------

fit_individuals = stan(file = "models/stan_pareto_truncated_individualsMEE.stan",
                       data = list(N = N_ind, 
                                   x = x_ind,
                                   xmin = xmin,
                                   xmax = xmax),
                       chains = 4,
                       iter = 2000)






# Compare to eightMethodsX -------------------------------------------------

# eightMethods.count
eight.data = tibble(Number = c, bodyMass = x, Year = 1980)
eight_result_counts <- eightMethods.count(data = eight.data)

# compare
stan_result_counts = fit_counts %>% as_draws_df() %>% median_qi(b) %>% 
  mutate(Method = "stan_counts",
         confMin = .lower, confMax = .upper) %>% 
  select(Method, b, confMin, confMax)

eight_result_counts %>% as_tibble() %>% bind_rows(stan_result_counts)  


# eightMethodsMEE
eight_resultMEE <- eightMethodsMEE(x = x_ind)

# compare
stan_result_individuals <- fit_individuals %>% as_draws_df() %>% median_qi(b) %>% 
  mutate(Method = "stan_individuals")

eight_result_individuals = tibble(b = eight_resultMEE$hMLE.list$b,
                                  .lower = eight_resultMEE$hMLE.list$confVals[1],
                                  .upper = eight_resultMEE$hMLE.list$confVals[2],
                                  Method = "eightMethodsMEE")

bind_rows(stan_result_individuals, eight_result_individuals)
