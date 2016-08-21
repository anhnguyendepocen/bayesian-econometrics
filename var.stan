# Vector Auto Regression (VAR)

data {
  int<lower=0> N;  // length of series
  int<lower=0> K;  // number of variables
  matrix[N,K] Y;   // time series
}
transformed data {
  matrix[N-1,K] Y_lead;
  matrix[N-1,K] Y_lag;
  for (n in 1:N-1) {
    Y_lead[n] = Y[n];
  }
  for (n in 2:N) {
    Y_lag[n-1] = Y[n];
  }
}
parameters {
  matrix[K,K] A;   // matrix of coefficients
  vector<lower=1>[K] sigma;
}
model {
  matrix[K,N-1] eta;
  eta = A*Y_lag';
  for (n in 1:N-1) {
    Y_lead'[,n] ~ normal(eta[,n], sigma);
  }
}
