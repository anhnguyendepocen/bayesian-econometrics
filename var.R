#########################################################
### Bayesian Vector Auto Regression (VAR) Using RStan ###
#########################################################

library("rstan")
library("vars")

### using the vars package ###
data(Canada)
VAR(Canada[,c(2,4)], p = 1, type = "none")
VAR(Canada[,c(2,4)], p = 1)

### using rstan ###
canada_data <- list()
canada_data$Y <- Canada[,c(1,2)]
canada_data$N <- dim(Canada[,c(1,2)])[1]
canada_data$K <- dim(Canada[,c(1,2)])[2]

canada_var <- VAR(Canada[,c(1,2)], p = 1)
predict_canada_var <- predict(canada_var)

stanc("/Users/imad/Desktop/git/bayesian-econometrics/var.stan")
stan_var <- stan(file = "/Users/imad/Desktop/git/bayesian-econometrics/var.stan", data = canada_data, cores = 4, chains = 4, iter = 2000)
traceplot(stan_var)

print(stan_var, pars = c("C","A","sigma","lp__"))
print(canada_var)

stan_var_post <- list()
stan_var_post$extract <- extract(stan_var)
stan_var_post$predict <- cbind(
  colMeans(stan_var_post$extract$Y_rep[,1,]),
  colMeans(stan_var_post$extract$Y_rep[,2,])
)

# plot generated quantities

par(mfrow=c(2,1))

plot(c(Canada[,1]), type = "l", lwd = 2, main = "Canadian Employment (1980-2000)", xlab = "Time (Quarterly)", ylab = "Employment", xaxt = "n")
axis(side = 1, at = seq(0, 83, 1), labels = NA)
axis(side = 1, at = seq(0, 80, 4), labels = seq(1980,2000), tck = 0, las = 3, cex.axis = 0.8)
lines(append(NA, stan_var_post$predict[,1]), col = "red", lwd = 2)

plot(c(Canada[,2]), type = "l", lwd = 2, main = "Canadian Productivity (1980-2000)", xlab = "Time (Quarterly)", ylab = "Productivity", xaxt = "n")
axis(side = 1, at = seq(0, 83, 1), labels = NA)
axis(side = 1, at = seq(0, 80, 4), labels = seq(1980,2000), tck = 0, las = 3, cex.axis = 0.8)
lines(append(NA, stan_var_post$predict[,2]), col = "cornflowerblue", lwd = 2)

dev.off()

var_model <- stan_model(file = "/Users/imad/Desktop/git/bayesian-econometrics/var.stan")
opt_var <- optimizing(var_model, data = canada_data, draws = 1000, algorithm = "Newton")
opt_var$par
coef(VAR(Canada[,c(2,4)], p = 1, type = "none"))
