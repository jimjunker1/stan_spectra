library(tidyverse)
library(tidybayes)
library(brms)



# Compute log-likelihood by hand. Conform that it recaptures the a --------
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





# Simulate regression varying intercepts paretocounts -----------------------------------------------------

xmin = 0.01
xmax = 1000

sim_params = tibble(beta = rnorm(5 -0, 1),
                    intercept = rnorm(5, -1.5, 0.1)) %>% 
  mutate(year = 1:nrow(.)) %>% 
  expand_grid(site = 1:5) %>% 
  mutate(site_no = as.integer(as.factor(paste0(year,"_", site)))) %>% 
  expand_grid(matsims = -2:2) %>% 
  mutate(b = intercept + rnorm(nrow(.),0, 0.1) + beta*matsims,
         index = 1:nrow(.))

sim_params %>% 
  ggplot(aes(x = matsims, y = b)) + 
  geom_smooth(aes(group = year), method = "lm", se = F) + 
  geom_point(aes(color = site))

sim_data = sim_params %>% 
  expand_grid(row = 1:100) %>%
  mutate(xmin = xmin,
         xmax = xmax,
         u = runif(nrow(.))) %>% 
  mutate(x_sim = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))) %>% # rPLB from Edwards 
  arrange(index, x_sim) %>% 
  group_by(index) %>% 
  mutate(rank = rank(x_sim))

sim_data_counts = sim_data %>% 
  group_by(xmin, xmax,index, matsims, year, site_no, site) %>%
  mutate(x_sim = round(x_sim, 4)) %>% 
  count(x_sim) 


# test stan model
mod_test = stan(file = "models/b_testnc_multiple_clusters_noncentered.stan", 
                data = list(N = nrow(sim_data_counts),
                            x = sim_data_counts$x_sim,
                            xmin = sim_data_counts$xmin,
                            xmax = sim_data_counts$xmax,
                            counts = sim_data_counts$n,
                            mat_s = sim_data_counts$matsims,
                            n_years = length(unique(sim_data_counts$year)),
                            n_sites = length(unique(sim_data_counts$site)),
                            site = as.integer(as.factor(sim_data_counts$site)),
                            year = as.integer(as.factor(sim_data_counts$year))),
                iter = 100, chains = 1)


posts_reg <- mod_test %>% 
  as_draws_df() %>% 
  expand_grid(mat_s = seq(-2, 2, length.out = 10)) %>% 
  group_by(.draw) %>% 
  mutate(sigma_sim_year = rnorm(1, 0, sigma_year),
         sigma_sim_site = rnorm(1, 0, sigma_site)) %>% 
  mutate(b_pred = a + mat_s*beta + sigma_sim_site + sigma_sim_year )

mod_test %>% 
  as_draws_df() %>% 
  pivot_longer(cols = contains("raw_site")) %>%
  mutate(bsite = a + value*sigma_site) %>% 
  group_by(name) %>% 
  median_qi(bsite)

posts_reg %>% 
  ggplot(aes(x = mat_s, y = b_pred)) + 
  geom_line(aes(group = .draw), alpha = 0.1) + 
  geom_point(data = sim_params, aes(x = matsims, y = b))


# Simulate non-centered regression varying intercepts paretocounts -----------------------------------------------------

xmin = 0.01
xmax = 1000
n_draws = 20

sim_params = tibble(beta = rnorm(n_draws, 0, 0.5),
                    intercept = rnorm(n_draws, -1.5, 0.5)) %>% 
  mutate(year = 1:nrow(.)) %>% 
  expand_grid(site = 1:5) %>% 
  mutate(site_no = as.integer(as.factor(paste0(year,"_", site)))) %>% 
  expand_grid(matsims = -2:2) %>% 
  mutate(b = intercept + beta*matsims,
         index = 1:nrow(.)) %>% 
  group_by(site) %>% 
  mutate(sigma_site = rexp(1, 1),
         sigma_year = rexp(1, 1),
         alpha_raw_site = rnorm(1, 0, 1),
         alpha_raw_year = rnorm(1, 0, 1)) %>% 
  mutate(b_year = b + sigma_year*alpha_raw_year,
         b_site = b + sigma_site*alpha_raw_site,
         b_year_site = b + sigma_site*alpha_raw_site + sigma_year*alpha_raw_year)

sim_params %>% 
  pivot_longer(cols = c(b, b_year, b_site, b_year_site)) %>% 
  ggplot(aes(x = matsims, y = value)) + 
  facet_wrap(~name, ncol = 4) +
  geom_smooth(aes(group = year), method = "lm", se = F) + 
  geom_point() + 
  ylim(-8, 8)

sim_data = sim_params %>%  
  expand_grid(row = 1:10) %>%
  mutate(xmin = xmin,
         xmax = xmax,
         u = runif(nrow(.))) %>% 
  mutate(x_sim = (u*xmax^(b+1) +  (1-u) * xmin^(b+1) ) ^ (1/(b+1))) %>% # rPLB from Edwards 
  arrange(index, x_sim) %>% 
  group_by(index) %>% 
  mutate(rank = rank(x_sim))

sim_data_counts = sim_data %>% 
  group_by(xmin, xmax,index, matsims, year, site_no, site) %>%
  mutate(x_sim = round(x_sim, 2)) %>% 
  count(x_sim) 

write_csv(sim_data_counts, file = "data/sim_data_counts.csv")

# test stan model
mod_test_varint = stan(file = "models/b_testnc_multiple_clusters_noncentered.stan", 
                data = list(N = nrow(sim_data_counts),
                            x = sim_data_counts$x_sim,
                            xmin = sim_data_counts$xmin,
                            xmax = sim_data_counts$xmax,
                            counts = sim_data_counts$n,
                            mat_s = sim_data_counts$matsims,
                            n_years = length(unique(sim_data_counts$year)),
                            n_sites = length(unique(sim_data_counts$site)),
                            site = as.integer(as.factor(sim_data_counts$site)),
                            year = as.integer(as.factor(sim_data_counts$year))),
                iter = 100, chains = 1,
                control=list(adapt_delta=0.95))


posts_reg <- mod_test %>% 
  as_draws_df() %>% 
  expand_grid(mat_s = seq(-2, 2, length.out = 10)) %>% 
  group_by(.draw) %>% 
  mutate(sigma_sim_year = rnorm(1, 0, sigma_year),
         sigma_sim_site = rnorm(1, 0, sigma_site)) %>% 
  mutate(b_pred = a + mat_s*beta + sigma_sim_site + sigma_sim_year )

mod_test %>% 
  as_draws_df() %>% 
  pivot_longer(cols = contains("raw_site")) %>%
  mutate(bsite = a + value*sigma_site) %>% 
  group_by(name) %>% 
  median_qi(bsite)

posts_reg %>% 
  ggplot(aes(x = mat_s, y = b_pred)) + 
  geom_line(aes(group = .draw), alpha = 0.1) + 
  geom_point(data = sim_params, aes(x = matsims, y = b))

