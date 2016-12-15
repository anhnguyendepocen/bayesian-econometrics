// Spatial Error Model (SEM)
data {
  int<lower=0> N;  // number of obs
  int<lower=0> K;  // number of predictors
  matrix[N,K] X;   // design matrix
  matrix[N,N] W;   // spatial weights matrix
  vector[N] y;     // response vector
}
transformed data {
  matrix[N,N] I;
  I = diag_matrix(rep_vector(1.0, N));
}
parameters {
  vector[K] beta;
  real<lower=0,upper=1> lambda;  // spatial correlation
  real<lower=0> sigma;
}
transformed parameters {
  cov_matrix[N] Sigma;
  Sigma = ((I - lambda * W) * (I - lambda * W')) / sigma;
}
model {
  target += multi_normal_prec_lpdf(y | X*beta, Sigma);
  target += beta_lpdf(lambda | 1, 2);  // prior on spatial correlation
}
