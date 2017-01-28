// least squares (heteroskedastic)
data {
  int<lower=1> N;
  int<lower=1> K;
  vector[N] y;
  matrix[N,K] X;
}
parameters {
  vector[K] beta;
  vector[N] sigma;
}
model {
  target += normal_lpdf(y | X * beta, sigma);
  target += normal_lpdf(beta | 0, 1);
  target += normal_lpdf(sigma | 0, 1);
}
