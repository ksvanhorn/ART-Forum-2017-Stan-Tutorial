data {
  int<lower=0> N1;  //obserations in data set 1
  int<lower=0> N2;  //observations in data set 2
  int<lower=0> K1;  //number of variables observed in data set 1 only
  int<lower=0> K2;  //number of variables observed in data set 2 only
  int<lower=0> Kb;  //number of variables observed in both data sets
  vector[K1] y1[N1];  
  vector[K2] y2[N2]; 
  vector[Kb] yb[N1 + N2];
}
parameters {
  vector[K1 + K2 + Kb] mu;
  corr_matrix[K1 + K2 + Kb] Omega;
  vector<lower=0>[K1 + K2 + Kb] tau;
  vector[K2] y1mis[N1];
  vector[K1] y2mis[N2];
} 
transformed parameters{
  vector[K1 + K2 + Kb] y[N1 + N2];
  for (n in 1:N1) {
    for (k in 1:K1) y[n][k] = y1[n][k];
    for (k in 1:K2) y[n][K1 + k] = y1mis[n][k];
    for (k in 1:Kb) y[n][K1 + K2 + k] = yb[n][k];
  }
  for (n in 1:N2) {
    for (k in 1:K1) y[N1+n][k] = y2mis[n][k];
    for (k in 1:K2) y[N1+n][K1+k] = y2[n][k];
    for (k in 1:Kb) y[N1+n][K1+K2+k] = yb[N1+n][k];
  }
}
model {
  //priors
  mu ~ normal(0, 100);  
  tau ~ cauchy(0,2.5);
  Omega ~ lkj_corr(2);
  //likelihood
  y ~ multi_normal(mu, quad_form_diag(Omega, tau));
}
