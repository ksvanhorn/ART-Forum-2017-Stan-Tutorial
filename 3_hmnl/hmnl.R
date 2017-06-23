# Code for estimating multinomial logit models in Stan
# Elea McDonnell Feit, eleafeit@gmail.com
# Last updated: 22 June 2017

# Copyright 2017, Kevin Van Horn & Elea McDonnell Feit 
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
#
# You may obtain a copy of the License at
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# ========== Load libraries and set options ==========
library(rstan)
library(MASS)
library(shinystan)

# writes a compiled Stan program to the disk to avoid recompiling
rstan_options(auto_write=TRUE) 
# allows Stan chains to run in parallel on multiprocessor machines
options(mc.cores = parallel::detectCores())

# ========== Test Stan with synthetic data ============

# function to generate mnl data
generate_hmnl_data <- function(R=100, S=30, C=3, 
                               Theta=matrix(rep(1, 8), nrow=2), 
                               Sigma=diag(0.1, 4)){
  K <- ncol(Theta)
  G <- nrow(Theta)
  Y <- array(dim=c(R, S))
  X <- array(rnorm(R*S*C*K), dim=c(R, S, C, K)) # normal covariates
  Z <- array(dim=c(G, R))
  Z[1,] <- 1  # intercept
  if (G > 1) {
    Z[2:G,] <- rnorm(R*(G-1)) # normal covariates
  }
  Beta <- array(dim=c(K, R))
  for (r in 1:R) {
    Beta[,r] <- mvrnorm(n=1, mu=Z[,r]%*%Theta, Sigma=Sigma)
    for (s in 1:S)
      Y[r,s] <- sample(x=C, size=1, prob=exp(X[r,s,,]%*%Beta[,r])) # logit formula
   }
  list(R=R, S=S, C=C, K=K, G=G, Y=Y, X=X, Z=Z, 
       beta.true=beta, Theta.true=Theta, Sigma.true=Sigma)
}

d1 <- generate_hmnl_data()
str(d1)

test.stan <- stan(file="hmnl.stan", data=d1, iter=1000, chains=4) 

plot(test.stan, plotfun="trace", pars=("Theta"))

plot(test.stan, plotfun="trace", pars=c("tau"))

plot(test.stan, plotfun="trace", pars=("Omega"))

plot(test.stan, pars=c("Theta", "tau", "Omega"))

plot(test.stan, pars=c("Theta", "Sigma"))

summary(test.stan, pars=c("Theta"))$summary

summary(test.stan, pars=c("tau"))$summary

summary(test.stan, pars=c("Omega"))$summary

# ========= Read in Chocolate Data and Prep for Stan ==========
rm(list=ls()) # tidy up
choc.df <- read.csv("cbc_chocolate.csv")

# Coding the chocolate data (this ought to be a function)
choc.contrasts <- list(Brand = "contr.sum", Type = "contr.sum")
choc.coded <- model.matrix(~ Brand + Type, data = choc.df, contrasts = choc.contrasts)
choc.coded <- choc.coded[,2:ncol(choc.coded)] # remove intercept
# Fix the bad labels from contr.sum
choc.names <- c("BrandDove", "BrandGhirardelli", "BrandGodiva", "BrandHersheys", 
                "TypeDark", "TypeDarkNuts", "TypeMilk", "TypeMilkNuts")
colnames(choc.coded) <- choc.names
choc.df <- cbind(choc.df, choc.coded)

head(choc.df)

# Munge into Stan list
R <- length(unique(choc.df$Ind))
S <- length(unique(choc.df$Trial))
C <- max(choc.df$Alt)
K <- 9
Y <- array(dim=c(R, S))
X <- array(rnorm(R*S*C*K), dim=c(R, S, C, K)) 
Z <- array(1, dim=c(1, R)) # intercept only
for (r in 1:R) { # respondents
  for (s in 1:S){ # choice scenarios
    scenario <- choc.df[choc.df$Ind==unique(choc.df$Ind)[r] & 
                        choc.df$Trial==unique(choc.df$Trial)[s], ]
    X[r,s,,] <- data.matrix(scenario[,c(7, 9:16)]) # price and coded brand and type
    Y[r,s] <- scenario$Alt[as.logical(scenario$Chosen)]
  }
}

choc.standata <- list(C=C, K=K, R=R, S=S, G=1, Y=Y, X=X, Z=Z)
str(choc.standata)

rm(Y, X, Z, R, S, r, s, C, K, choc.contrasts, scenario, choc.coded)

# ========== Run Stan and Check Convergence =========
choc.stan <- stan(file="hmnl.stan", data=choc.standata)

plot(choc.stan, plotfun="trace", pars=("Theta"))

plot(choc.stan, plotfun="trace", pars=c("tau"))

plot(choc.stan, plotfun="trace", pars=c("Omega[1,2]"))

plot(choc.stan, plotfun="trace", pars=paste("Beta[", 1:9, ",1]", sep="")) # resp 1

summary(choc.stan)$summary[,c("Rhat", "n_eff")]

# Convergence checking gets tedious with a lot of parameters, so automate
check_fit <- function(fit) {
  summ <- summary(fit)$summary
  range_rhat <- range(summ[ , 'Rhat'])
  rhat_ok <- 0.99 <= range_rhat[1] && range_rhat[2] <= 1.1
  range_neff <- range(summ[ , 'n_eff'])
  neff_ok <- range_neff[1] >= 400
  sp <- rstan::get_sampler_params(fit, inc_warmup=FALSE)
  max_divergent <- max(sapply(sp, function(p){ sum(p[ , 'divergent__']) }))
  no_divergent <- max_divergent == 0
  
  list(ok = rhat_ok && neff_ok && no_divergent,
       range_rhat = range_rhat,
       range_neff = range_neff,
       max_divergent = max_divergent)
}

check_fit(choc.stan)

# ========== Summarize Posterior ==========
choc.names
plot(choc.stan, pars=c("Theta", "tau"))

plot(choc.stan, pars=c("Omega"))

plot(choc.stan, pars=paste("Beta[", 1:9, ",1]", sep="")) + ggtitle("Respondent 1: Likes Milk Chocolate")

plot(choc.stan, pars=paste("Beta[", 1:9, ",2]", sep="")) + ggtitle("Respondent 2: Likes Dark Chocolate")

launch_shinystan(choc.stan)

# ========== Simulate Shares from Model =========
# function for computing shares from beta.draws for an hmnl model
shares.hmnl.post <- function(beta.draws, X) # X is attribute matrix for scenario
{
  R <- dim(beta.draws)[3]  # respondents
  D <- dim(beta.draws)[1]  # draws
  shares <- array(NA, dim=c(nrow(X), R, D))
  for (d in 1:D) {
    for (r in 1:R) {
      beta <- beta.draws[d,,r] 
      V <- exp(X %*% beta)  
      shares[,r,d] <- V/sum(V)
    }
  }
  shares
}

choc.standata$X[1,1,,] # note option 2 is dark
shares <- shares.hmnl.post(extract(choc.stan, pars=c("Beta"))$Beta, 
                          choc.standata$X[1,1,,])
str(shares)

apply(shares, 1, quantile, probs=c(0.5, 0.025, 0.975))
apply(shares, 1:2, quantile, probs=c(0.5, 0.025, 0.975))

