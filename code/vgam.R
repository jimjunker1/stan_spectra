library(loo)

log_lik1 <- extract_log_lik(fit_temponly, parameter_name = "lp__")
waic1 <- loo(fit_temponly)
loo::waic(log_lik)


log_lik <- extract(fit_temponly)$lp__

?dnorm()

dnorm(2, 1, 0.2, log = TRUE)


x = 0.002
lower = 0.001
upper = 10000
shape = 1.4


library(VGAM)
dtruncpareto(x, lower, upper, shape, log = FALSE)

N = 1000
x = rtruncpareto(N, lower, upper, shape)
counts = rep(1, N)

fit_model = stan_model("models/old/b_paretocounts_singlesample.stan")

stan_dat <- list(x = x,
                   N = N,
                   counts = ,
                   xmax = upper,
                   xmin = lower)

fit <- sampling(object = fit_model,
                  data = stan_dat,
                  iter = 1000,
                  chains = 2,
                  open_progress = F,
                  verbose = F)
  

