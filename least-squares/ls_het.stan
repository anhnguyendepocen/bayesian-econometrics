data {
  int<lower=0> N;
  int<lower=0> K;
  matrix[N,K] X;
  vector[N] y;
}
parameters {
  vector[K] beta;
  vector<lower=0>[N] sigma;
}
model {
  target += normal_lpdf(y | X*beta, sigma);
  target += normal_lpdf(beta | 0, 1);
  target += normal_lpdf(sigma | 0, 1);
}
