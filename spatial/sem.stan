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
model {
  matrix[N,N] Sigma;
  matrix[N,N] weight_stuff;
  weight_stuff = I - lambda * W;
  Sigma = tcrossprod(weight_stuff) / sigma;

  target += multi_normal_prec_lpdf(y | X*beta, Sigma);
  target += cauchy_lpdf(lambda | 0, 2);  // prior on spatial correlation
  target += cauchy_lpdf(beta | 0, 3);
  target += cauchy_lpdf(sigma | 0, 3);
}

