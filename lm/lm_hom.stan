// Least Squares (homoskedastic)
data {
  int<lower=1> N;     // observations
  int<lower=1> K;     // predictors
  vector[N] y;        // outcome variable
  matrix[N,K] X;      // design matrix
}
parameters {
  vector[K] beta;       // coefficients on predictors
  real<lower=0> sigma;  // variance of y (constant)
}
model {
  target += normal_lpdf(y | X * beta, sigma);
  // priors
  target += normal_lpdf(beta | 0, 1);
  target += normal_lpdf(sigma | 0, 1);
}
