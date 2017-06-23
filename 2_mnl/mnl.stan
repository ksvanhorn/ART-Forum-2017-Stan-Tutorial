# Multinomial logit model
# Elea McDonnell Feit, eleafeit@gmail.com

data {
	int<lower=2> C; // # of alternatives (choices) in each scenario
	int<lower=1> N;	// # of observations
	int<lower=1> K; // # of covariates
	int<lower=1,upper=C> Y[N]; // observed choices
	matrix[C,K] X[N]; // matrix of attributes for each obs
}

parameters {
	vector[K] beta;
}

model {
  # priors
	// beta ~ normal(0,3); // often used in Gibbs sampling
	beta ~ cauchy(0, 2.5);  // better
	# model
	for (i in 1:N)
		Y[i] ~ categorical_logit(X[i]*beta);
}

