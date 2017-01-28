library(rstan)
library(bayesplot)

N <- 500
x1 <- rnorm(N, 2, 4)
x2 <- rnorm(N, 0, 1)
y <- 2 + 0.5*x1 - 1.5*x2 + rnorm(N, 0, 3)
sd_het <- runif(N, 1,1000)
z <- 2 + 0.5*x1 - 1.5*x2 + rnorm(N, 0, sd_het)
w <- rep(NA, N)
w <- 2 + 0.5*x1 - 1.5*x2 + rnorm(N, 0, c(rep(3,250),rep(1,250),rep(10,250),rep(2,250)))
sd_cor <- mvtnorm::rmvnorm(1, rep(0,500), diag(500))
sd_cor <- t(sd_cor)

dat <- list(y = y, X = cbind(rep(1,500), x1, x2))
dat$K <- ncol(dat$X)
dat$N <- nrow(dat$X)

fit_hom <- stan("ls_hom.stan", data = dat, chains = 4, iter = 1000, cores = 4)
print(fit_hom, digits = 2)

dat$y <- z
fit_het <- stan("ls_het.stan", data = dat, chains = 4, iter = 2000, cores = 4)
print(fit_het, digits = 2)
fit_clu <- lm(w ~ x1 + x2)

post_het <- extract(fit_het, inc_warmup = TRUE, permuted = FALSE, pars = c("beta", "sigma"))
color_scheme_set("mix-blue-pink")
mcmc_trace(post_het, n_warmup = 1000, regex_pars = c("beta"), facet_args = list(nrow = 2, labeller = label_parsed))

plot(fit_hom$residuals)
plot(fit_het$residuals)
plot(fit_clu$residuals)

qqnorm(fit_hom$residuals)
qqnorm(fit_het$residuals)
qqnorm(fit_clu$residuals)