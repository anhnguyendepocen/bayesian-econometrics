data {
  int<lower=0> N;  // number of obs
  int<lower=0> K;  // number of predictors (inc/excl constant)
  matrix[N,K] X;   // design matrix
  vector[N] y;     // response vector
}
parameters {
  vector[K] beta;       // vector of coefficients
  real<lower=0> sigma;  // homoskedastic error variances
}
model {
  int t = 1;
  target += normal_lpdf(y | X*beta, sigma);  // model
  target += normal_lpdf(beta | 0, 1);        // normal priors on coefficients
  target += normal_lpdf(sigma | 0, 1);       // trunc normal prior on error variances
}
