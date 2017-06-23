# Advanced techniques with lm in Stan:
#  - Testing with synthetic data 
#  - Comparing to other inferential engines in R
# Elea McDonnell Feit, eleafeit@gmail.com
# 14 June 2017

# ======== Test with synthetic data ===========
library(MASS)  # multivariate distributions for generating data

generate_lm_data <- function(N=100, # number of observations
                             mean.X=c(0, 0), Sigma.X=diag(1,2), # mean and covariance of ivs
                             beta.0=0, beta=c(1,-1), sigma.e=1) # parameters of linear model
{ 
  epsilon <- rnorm(n=N, mean=0, sd=sigma.e)
  X <- mvrnorm(n=N, mu=mean.X, Sigma=Sigma.X) 
  colnames(X) <- paste("X", 1:(ncol(X)), sep="")
  y <- beta.0 + X %*% beta + epsilon
  list(data=list(N=N, K=length(beta), y=as.vector(y), X=X), 
       true=list(beta.0=beta.0, beta=beta, sigma.e=sigma.e))  
}

set.seed(20130601)
d <- generate_lm_data(N=100)$data
str(d)

m.stan <- stan(file="lm.stan", data=d, iter=1000) # read model in from file
summary(m.stan)$summary
plot(m.stan)
get_stanmodel(m.stan)

# fancier visualizations of posteriors
library(bayesplot)
mcmc_combo(As.mcmc.list(m.stan))

# ========= Compare to  other inferential methods ==========

# utility function for combining y and X in a data frame
frameup <- function(data) {
  data.frame(y=data$y, data$X)
}

# compare to lm(), which estimates the model by OLS
m.lm <- lm(y ~ ., data=frameup(d))
summary(m.lm)

# compare to MCMCregress(), which estimates the model using Bayesian methods
library(MCMCpack)  # alternative Bayesian estimation for linear model
m.MCMCpack <- MCMCregress(y ~ ., data=frameup(d))
summary(m.MCMCpack)
mcmc_combo(as.mcmc.list(m.MCMCpack))

# pull out means and posterior std deviations for comparison
data.frame(stan=summary(m.stan)$summary[1:4,c(1,3)], 
           lm=rbind(coef(summary(m.lm))[1:3,1:2], c(NA, NA)), 
           MCMCpack=summary(m.MCMCpack)$statistics[1:4,1:2])

