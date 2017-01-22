// Spatial Autocorrelation Model (SEM)
// try demeaning predictors
// try non centered vs centered parameterization
data {
  int<lower=0> N;           // number of obs
  int<lower=0> K;           // number of predictors
  matrix[N,K] X;            // design matrix
  matrix<lower=0>[N,N] W;   // spatial weight matrix
  vector[N] y;              // response vector
}
transformed data {
  matrix[N,N] I;
  I = diag_matrix(rep_vector(1.0, N));
}
parameters {
  real<lower=0,upper=1> rho;  // note the support is (0,1) not (-1,1)
  vector[K] beta;             // parameter vector
  real<lower=0> sigma;        // variance
}
model {
  matrix[N,N] Sigma;
  matrix[N,N] weight_stuff;
  weight_stuff = I - rho * W;
  Sigma = tcrossprod(weight_stuff) /sigma;

  target += multi_normal_prec_lpdf(y | inverse(weight_stuff) * X * beta, Sigma);
  target += cauchy_lpdf(rho | 0, 2);       // prior on spatial autocorrelation
  target += cauchy_lpdf(beta | 0, 3);     // priors on predictor parameters
  target += cauchy_lpdf(sigma | 0, 3);    // prior on variation
  for (i in 1:N)
    target += cauchy_lpdf(Sigma[i,] | 0, 3);    // prior on covariance matrix
}

