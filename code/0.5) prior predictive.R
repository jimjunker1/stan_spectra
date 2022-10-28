library(tidyverse)

# prior predictive
N = 1000
prior_sim <- tibble(beta = rnorm(N, 0, 0.1),
                    sigma_year = abs(rnorm(N, 0, 0.1)),
                    sigma_site = abs(rnorm(N, 0, 0.1)),
                    alpha_year_raw = rnorm(N, 0, 5),
                    alpha_site_raw = rnorm(N, 0, 5),
                    a = rnorm(N, -1.5, 0.2),
                    .draw = 1:N) %>% 
  mutate(intercept = a + sigma_year*alpha_year_raw + sigma_site*alpha_site_raw) %>% 
  expand_grid(mat_s = seq(-2, 2, length.out = 20)) %>% 
  mutate(y_pred = intercept + beta*mat_s)

prior_sim %>% 
  ggplot(aes(x = mat_s, y = y_pred, group = .draw)) + 
  geom_line(alpha = 0.2) + 
  labs(y = "y_pred (of the b exponent)",
       x = "Stream mean annual temperature (standardized)")


