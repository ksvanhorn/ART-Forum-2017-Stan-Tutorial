// Stan model code for linear model
data {
  int N;  
  int K;
  vector[N] y;
  matrix[N,K] X;
}
parameters {
  real beta0;
  vector[K] beta;
  real<lower=0> sigma;
}
model {
  y ~ normal(beta0 + X*beta, sigma);
}