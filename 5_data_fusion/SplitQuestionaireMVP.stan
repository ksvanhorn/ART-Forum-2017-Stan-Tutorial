functions {
  int mysum(int[,] a) {
    int s;
    s = 0;
    for (i in 1:size(a))
      s = s + sum(a[i]);
    return s;
  }
} 
data {
  int<lower=0> K1;     // number of vars only observed in data set 1
  int<lower=0> K2;     // number of vars only observed in data set 2
  int<lower=0> Kb;     // number of vars observed in both data sets
  int<lower=0> N1;     // number of observations in data set 1
  int<lower=0> N2;     // number of observations in data set 2     
  int<lower=0,upper=2> y[N1+N2, K1+K2+Kb];  // should contain zeros in missing positions
}
transformed data {
  int<lower=1, upper=N1+N2> n_pos[mysum(y)];
  int<lower=1, upper=K1+K2+Kb> k_pos[size(n_pos)];
  int<lower=1, upper=N1+N2> n_neg[(N1+N2)*(K1+K2+Kb) - K2*N1 - K1*N2 - mysum(y)];
  int<lower=1, upper=K1+K2+Kb> k_neg[size(n_neg)];
  int<lower=0> N_pos; 
  int<lower=0> N_neg; 
  N_pos = size(n_pos);
  N_neg = size(n_neg);
  {
    int i; 
    int j;
    i = 1; 
    j = 1; 
    for (n in 1:N1) {                      //positions in observed y1
      for (k in 1:K1) { 
        if (y[n,k] == 1) {
          n_pos[i] = n; 
          k_pos[i] = k; 
          i = i + 1;
        } else {
          n_neg[j] = n; 
          k_neg[j] = k;
          j = j + 1;
        }
      }
      for (k in (K1+K2+1):(K1+K2+Kb)) {
        if (y[n,k] == 1) {
          n_pos[i] = n; 
          k_pos[i] = k; 
          i = i + 1;
        } else {
          n_neg[j] = n; 
          k_neg[j] = k;
          j = j + 1;
        }
      }
    }
    for (n in (N1+1):(N1+N2)) {               //positions in observed y2
      for (k in (K1+1):(K1+K2+Kb)) { 
        if (y[n,k] == 1) {
          n_pos[i] = n; 
          k_pos[i] = k; 
          i = i + 1;
        } else {
          n_neg[j] = n; 
          k_neg[j] = k;
          j = j + 1;
        }
      }
    }
  }
}
parameters {
  vector[K1 + K2 + Kb] mu;
  corr_matrix[K1 + K2 + Kb] Omega;
  vector<lower=0>[N_pos] z_pos;
  vector<upper=0>[N_neg] z_neg;
  vector[K2] z1mis[N1];
  vector[K1] z2mis[N2];
} 
transformed parameters{
  vector[K1 + K2 + Kb] z[N1 + N2]; 
  vector[K2] y1mis[N1]; 
  vector[K1] y2mis[N2];
  for (i in 1:N_pos) 
    z[n_pos[i], k_pos[i]] = z_pos[i];
  for (i in 1:N_neg)
    z[n_neg[i], k_neg[i]] = z_neg[i];
  for (n in 1:N1) {
    for (k in 1:K2) {
      z[n, K1 + k] = z1mis[n, k];
      if (z1mis[n, k] > 0)
        y1mis[n, k] = 1; 
      if (z1mis[n, k] < 0)
        y1mis[n, k] = 0;
    }
  }
  for (n in 1:N2) {
    for (k in 1:K1) {
      z[N1 + n, k] = z2mis[n, k];
      if (z2mis[n, k] > 0)
        y2mis[n, k] = 1; 
      if (z2mis[n, k] < 0)
        y2mis[n, k] = 0;
    }
  }
}
model {
  mu ~ normal(0, 3);  
  Omega ~ lkj_corr(1);
  z ~ multi_normal(mu, Omega);
}

