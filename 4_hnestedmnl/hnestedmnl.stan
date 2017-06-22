# Copyright (C) 2015, The Modellers / Hall and Partners
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# version 2, as published by the Free Software Foundation.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
data {
  int<lower=1> R; // # of respondents
  int<lower=1> K; // # of product covariates; no intercept
  int<lower=1> G; // # of respondent covariates
  int<lower=1> S; // # of scenarios per respondent
  int<lower=2> C; // # of alternatives (choices) per scenario
  matrix[C, K] X[R, S]; // X[r,s] is covariate matrix of scenario s for respondent r.
  matrix[G, R] Z; // Z[,r] is vector of covariates describing respondent r
  int<lower=1,upper=C> Y1[R,S];  // forced-choice responses
  int<lower=0,upper=1> Y2[R,S];  // whether respondent would choose any of the options
}
parameters {
  real<lower=0, upper=1> lambda;
  vector<lower=0>[K+1] sigma_raw;
  matrix[K+1, G] Theta_raw;
  matrix[K+1, R] Epsilon;
}
transformed parameters {
  vector<lower=0>[K+1] sigma = 5 * sigma_raw;
  matrix[K+1, G] Theta = 10 * Theta_raw;
  matrix[K+1, R] Beta = Theta * Z + diag_pre_multiply(sigma, Epsilon);
}
model {
  lambda ~ uniform(0, 1);
  sigma_raw ~ normal(0, 1);
  to_vector(Theta_raw) ~ normal(0, 1);
  to_vector(Epsilon) ~ normal(0, 1);
  for (r in 1:R) {
    vector[K] b = Beta[2:(K+1), r];
    real alpha = Beta[1,r];
    for (s in 1:S) {
      vector[C] u = X[r,s] * b;
    	real u_buy = alpha + lambda * log_sum_exp(u);
    	Y1[r,s] ~ categorical_logit(u);
      Y2[r,s] ~ bernoulli_logit(u_buy);
    }
  }
}
