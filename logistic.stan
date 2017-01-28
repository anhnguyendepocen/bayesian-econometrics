// logistic regeression
data {
  int<lower=1> N;               // observations
  int<lower=1> K;               // predictors
  int<lower=0, upper=1> y[N];   // outcome variable
  matrix[N,K] X;                // design matrix
}
parameters {
  vector[K] beta;     // coefficients
}
model {
  target += bernoulli_logit_lpmf(y | X * beta);
  // target += bernoulli_lpmf(y | inv_logit(X * beta)); // equivalent to above
  // priors
  target += normal_lpdf(beta | 0, 1);
}
