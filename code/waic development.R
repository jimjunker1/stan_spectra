library(brms)
library(tidyverse)
library(janitor)


iris

test = brm(Sepal.Length ~ Sepal.Width + (1|Species),
           data = iris,
           family = gaussian(),
           iter = 1000, 
           chain = 1)


waic(test, re_formula = NULL)


# by hand ------------------
posts = as_draws_df(test)

post_mu = test$data %>% 
  mutate(setosa01 = case_when(grepl("setosa", Species) ~ 1, TRUE ~ 0),
         versicolor01 = case_when(grepl("versicolor", Species) ~ 1, TRUE ~ 0),
         virginica01 = case_when(grepl("virginica", Species) ~ 1, TRUE ~ 0)) %>% 
  rownames_to_column() %>% 
  expand_grid(posts) %>% 
  janitor::clean_names() %>% 
  glimpse() %>% 
  mutate(mu = b_intercept + b_sepal_width*sepal_width + 
           r_species_setosa_intercept*setosa01+ 
           r_species_versicolor_intercept*versicolor01+ 
           r_species_virginica_intercept*virginica01)

post_logprob = post_mu %>% 
  mutate(logprob = dnorm(sepal_length, mu, sigma, log = T))

post_lppd = post_logprob %>% 
  group_by(rowname) %>% 
  mutate(n_samples = max(draw)) %>% 
  summarize(lppd = rethinking::log_sum_exp(logprob) - log(n_samples),
            pwaic = var(logprob)) %>% 
  distinct(lppd, rowname, pwaic)

-2*(sum(post_lppd$lppd) - sum(post_lppd$pwaic))




# by hand without data matrix of 0,1's ------------------------------------

posts = as_draws_df(test)

post_mu = test$data %>% 
  # mutate(setosa01 = case_when(grepl("setosa", Species) ~ 1, TRUE ~ 0),
  #        versicolor01 = case_when(grepl("versicolor", Species) ~ 1, TRUE ~ 0),
  #        virginica01 = case_when(grepl("virginica", Species) ~ 1, TRUE ~ 0)) %>% 
  rownames_to_column() %>% 
  expand_grid(posts) %>% 
  janitor::clean_names() %>% 
  glimpse() %>% 
  pivot_longer(cols = contains("r_species")) %>% 
  mutate(flag = str_detect(name, as.character(species))) %>% 
  filter(flag == T) %>% 
  mutate(mu = b_intercept + b_sepal_width*sepal_width +
           value)

post_logprob = post_mu %>% 
  mutate(logprob = dnorm(sepal_length, mu, sigma, log = T))

post_lppd = post_logprob %>% 
  group_by(rowname) %>% 
  mutate(n_samples = max(draw)) %>% 
  summarize(lppd = rethinking::log_sum_exp(logprob) - log(n_samples),
            pwaic = var(logprob)) %>% 
  distinct(lppd, rowname, pwaic)

-2*(sum(post_lppd$lppd) - sum(post_lppd$pwaic))





post_lppd = loglik_fit %>% 
  group_by(rowname) %>% 
  mutate(n_samples = max(.draw)) %>% 
  summarize(lppd = rethinking::log_sum_exp(loglik) - log(n_samples),
            pwaic = var(loglik)) %>% 
  distinct(lppd, rowname, pwaic)

-2*(sum(post_lppd$lppd) - sum(post_lppd$pwaic))







