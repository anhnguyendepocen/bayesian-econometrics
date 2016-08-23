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
    Y_lag[n] = Y[n];
  }
  for (n in 2:N) {
    Y_lead[n-1] = Y[n];
  }
  print("Y_lead = ", Y_lead);
  print("Y_lag = ", Y_lag);
}
parameters {
  vector[K] C;              // constant for each equation
  matrix[K,K] A;            // matrix of coefficients
  vector<lower=0>[K] sigma; // standard deviation for each equation
}
model {
  matrix[K,N-1] eta;
  for (k in 1:K) {
    eta[k,] = C[k] + (A*Y_lag')[k,];    // linear predictor
  }
  #eta = C+(A*Y_lag');
  for (n in 1:N-1) {
    // Y_lead'[,n] ~ normal(eta[,n], sigma);              // model
    target += normal_lpdf(Y_lead'[,n] | eta[,n], sigma);  // new syntax
  }
  // include priors
  for (k in 1:K){
    C[k] ~ normal(0,1);
    for (j in 1:K) {
      A[j,k] ~ normal(0,1);
    }
  }
}
generated quantities {
  matrix[K,N-1] Y_rep;
  matrix[K,N-1] eta_rep;
  for (k in 1:K) {
    eta_rep[k,] = C[k] + (A*Y_lag')[k,];
  }
  for (n in 1:N-1) {
    for (k in 1:K) {
      Y_rep[k,n] = normal_rng(eta_rep[k,n], sigma[k]);  // predictions
    }
  }
}
