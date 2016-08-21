#########################################################
### Bayesian Vector Auto Regression (VAR) Using RStan ###
#########################################################

library("rstan")
library("vars")

### using the vars package ###
data(Canada)
coef(VAR(Canada[,c(2,4)], p = 1, type = "none"))
VAR(Canada[,c(2,4)], p = 1)

### using rstan ###
canada_data <- list()
canada_data$Y <- Canada[,c(2,4)]
canada_data$N <- dim(Canada[,c(2,4)])[1]
canada_data$K <- dim(Canada[,c(2,4)])[2]
stanc("/Users/imad/Desktop/git/bayesian-econometrics/var.stan")
stan_var <- stan(file = "/Users/imad/Desktop/git/bayesian-econometrics/var.stan", data = canada_data, cores = 4, chains = 4, iter = 2000)
print(stan_var, digits = 3)
