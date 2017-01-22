setwd("~/Desktop/git/bayesian-econometrics/spatial")

library(rstan)
library(spdep)
library(mvtnorm)
library(maptools)
library(bayesplot)

chi.poly <- readShapePoly("foreclosures.shp")
chi.poly <- chi.poly[1:5,]                                  # use a smaller data set
chidata <- chi.poly@data                                     # extract data from the spatial polygon
list.queen<-poly2nb(chi.poly, queen=TRUE)                    # if queen = TRUE not specified then rook method used for neighbors
W_list <- nb2listw(list.queen, style="W", zero.policy=TRUE)  # style='W' -> row standardized
W <- listw2mat(W_list)                                       # creates the weight matrix used in stan

plot(readShapePoly("foreclosures.shp"))
plot(W_list, coordinates(chi.poly), points = FALSE, add = TRUE, col = "red")

X <- cbind(rep(1,nrow(chidata)), log(chidata$est_fcs_rt), log(chidata$totpop))
# X <- cbind(rep(1,nrow(chidata)), chidata$est_fcs_rt, chidata$bls_unemp)
y <- log(chidata$violent)

N <- nrow(X)
K <- ncol(X)

dat <- list("N" = N, "K" = K, "X" = X, "y" = y, "W" = W)

### SAR STAN

sar <- stan_model("sar.stan")
start <- proc.time()
fit_sar <- sampling(sar, data = dat, chains = 4, iter = 500, cores = 4)
end <- proc.time()
end["elapsed"] - start["elapsed"]
print(fit_sar, digits = 2, pars = c("lambda", "beta", "sigma"))
opt_sar <- optimizing(sar, data = dat, iter = 10000, verbose = TRUE)
print(opt_sar, digits = 2)

post_sar <- extract(fit_sar, inc_warmup = TRUE, permuted = FALSE, pars = c("lambda", "beta", "sigma"))
color_scheme_set("mix-blue-pink")
mcmc_trace(post_sar, n_warmup = 50, facet_args = list(nrow = 2, labeller = label_parsed))

### compare with spdep

# SAR
sar.chi<-lagsarlm(log(chidata$violent) ~ log(chidata$est_fcs_rt) + log(chidata$totpop), data=chi.poly@data, W_list)
summary(sar.chi)
coef(sar.chi)
# SAR using 2SLS
sar2sls.chi<-stsls(violent~est_fcs_rt+bls_unemp, data=chi.poly@data, W_list)
summary(sar2sls.chi)

### SEM STAN

sem <- stan_model("sem.stan")

start <- proc.time()
fit_sem <- sampling(sem, data = dat, chains = 4, iter = 2000, cores = 4)
end <- proc.time()
end["elapsed"] - start["elapsed"]
print(fit_sem, digits = 2, pars = c("beta","lambda", "sigma"))

opt_sem <- optimizing(sem, data = dat, iter = 2000, verbose = TRUE)
print(opt_sem)

post_sem <- extract(fit_sem, inc_warmup = TRUE, permuted = FALSE, pars = c("beta","lambda", "sigma"))
color_scheme_set("mix-blue-pink")
mcmc_trace(post_sem, n_warmup = 500, facet_args = list(nrow = 2, labeller = label_parsed))

### compare with spdep

errorsalm.chi<-errorsarlm(I(log(violent))~I(log(est_fcs_rt))+I(log(totpop)), data=chi.poly@data, W_list)
summary(errorsalm.chi)
coef(errorsalm.chi)


### Simulated Data SAR

N <- 50
W_bin <- matrix(rep(0, N * N), nrow = N)
W_bin[lower.tri(W_bin)] <- rbinom(choose(N,2), 1, 0.5)
W_bin <- W_bin + t(W_bin)
W <- apply(W_bin, 2, function(x){x/rowSums(W_bin)})
I <- diag(N)
lambda <- 0.5
sigma <- 0.3
Sigma <- solve((I - lambda * W) %*% (I - lambda * t(W))) * sigma
X <- cbind(rep(1,N),rnorm(N, 0, 1), rnorm(N, 3, 1))
beta <- c(3, 2.5, -1.5)
mu <- solve(I - lambda * W) %*% X %*% beta
y <- c(mvtnorm::rmvnorm(1, mu, Sigma))

dat <- list("N" = nrow(X), "K" = ncol(X), "y" = y, "X" = X, "W" = W)

start <- proc.time()
fit_sar <- stan(file = "sar.stan", data = dat, chains = 4, iter = 2000, cores = 4)
end <- proc.time()
cat("model run time:", round((end["elapsed"] - start["elapsed"])/60,2), "minutes")
print(fit_sar, digits = 2)
traceplot(fit_sar, inc_warmup = TRUE)

check_sar <- spdep::lagsarlm(y ~ x1 + x2, data = data.frame(cbind("y" = dat$y, "x1" = dat$X[,2], "x2" = dat$X[,3])),
                      listw = mat2listw(W, style = "W"))
coef(check_sar)
