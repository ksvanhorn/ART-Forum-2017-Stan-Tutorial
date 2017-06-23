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
library(shinystan)

# writes a compiled Stan program to the disk to avoid recompiling
rstan_options(auto_write=TRUE) 
# allows Stan chains to run in parallel on multiprocessor machines
options(mc.cores = parallel::detectCores())

# ========== Test Stan with synthetic data ============

# function to generate mnl data
generate_mnl_data <- function(N=1000, C=3, beta=c(1, -2)){
  K <- length(beta)
  Y <- rep(NA, N)
  X <- list(NULL) 
  for (i in 1:N) {
    X[[i]] <- matrix(rnorm(C*K), ncol=K) # normal covariates  
    Y[i] <- sample(x=C, size=1, prob=exp(X[[i]]%*%beta)) # logit formula
  }
  list(N=N, C=C, K=K, Y=Y, X=X)
}

d1 <- generate_mnl_data(N=1000, C=3, beta=c(1, -2))
str(d1)

test.stan <- stan(file="mnl.stan", data=d1, iter=1000, chains=4) 
summary(test.stan)$summary

summary(test.stan)$c_summary

plot(test.stan, plotfun="trace")

plot(test.stan)

plot(test.stan, plotfun="hist")

plot(test.stan, plotfun="dens")

# ========= Read in and inspect Chocolate data =========
choc.df <- read.csv("cbc_chocolate.csv")
head(choc.df)

summary(choc.df)

mosaicplot(~ Brand + Chosen, data=choc.df)

mosaicplot(~ Type + Chosen, data=choc.df)

mosaicplot(~ Price + Chosen, data=choc.df)

# ========== Prep the Chocolate data for Stan =========
# illustrating alternative coding schemes (not required)
choc.df[1:10,]
model.matrix(~ Brand, data = choc.df)[1:10,]
model.matrix(~ Brand, data = choc.df, contrasts = list(Brand = "contr.sum"))[1:10,]

# coding the chocolate data
choc.contrasts <- list(Brand = "contr.sum", Type = "contr.sum")
choc.coded <- model.matrix(~ Brand + Type, data = choc.df, contrasts = choc.contrasts)
choc.coded <- choc.coded[,2:ncol(choc.coded)] # remove intercept
# Fix the bad labels from contr.sum
choc.names <- c("BrandDove", "BrandGhirardelli", "BrandGodiva", "BrandHersheys", 
                "TypeDark", "TypeDarkNuts", "TypeMilk", "TypeMilkNuts")
colnames(choc.coded) <- choc.names
choc.df <- cbind(choc.df, choc.coded)
head(choc.df)

# munging into stan format
unique(choc.df$Ind)
unique(choc.df$Trial)

R <- length(unique(choc.df$Ind))
S <- length(unique(choc.df$Trial))
Y <- rep(NA, R*S)
X <- vector("list", R*S)
n <- 1
for (r in unique(choc.df$Ind)) { # respondents
  for (s in unique(choc.df$Trial)){ # choice scenarios
     scenario <- choc.df[choc.df$Ind==r & choc.df$Trial==s,]
     X[[n]] <- data.matrix(scenario[,c(7, 9:16)]) # price and coded brand and type
     Y[n] <- scenario$Alt[as.logical(scenario$Chosen)]
     n <- n + 1
  }
}

str(Y)
str(X)

choc.standata <- list(N=length(Y), C=3, K=9, Y=Y, X=X)
rm(Y, X, n, R, S, r, s, choc.contrasts, scenario, choc.coded)

# ========== Estimate mnl model for Chocolate data in Stan ==========
choc.stan <- stan(file="mnl.stan", data=choc.standata) # default: chains=4, iter=2000

summary(choc.stan)$summary

plot(choc.stan, plotfun="trace")

plot(choc.stan)

data.frame(params=c("Price", choc.names), check.names=FALSE, 
           summary(choc.stan, pars=c("beta"))$summary)

# ========= Explore posterior using ShinyStan ==========
launch_shinystan(choc.stan)

# ========= Share prediction based on model ==========
# computing shares from a point estimate of beta (not good practice)
shares.mnl.point <- function(beta,  # vector of parameters (part-worths)
                             X) {   # attribute matrix X for scenario (coded)
  if (length(beta) != ncol(X)) 
    stop("length of beta doesn't match columns in X")
  V <- exp(X %*% beta)
  data.frame(shares=V/sum(V), X)
}


beta.mean <- summary(choc.stan)$summary[1:9,1]
beta.mean
shares.mnl.point(beta.mean, choc.standata$X[[1]]) 

# computing the posterior draws of shares 
shares.mnl.post <- function(draws,  # matrix of draws (use extract())
                            X) {    # attribute matrix X for scenario
  shares <- matrix(NA, nrow=nrow(draws), ncol=nrow(X))
  for (draw in 1:nrow(draws)) {
     shares[draw,] <- shares.mnl.point(draws[draw,], X)[,1]
  }
  shares
}

beta.draws <- extract(choc.stan, pars="beta")$beta
head(beta.draws)
shares.draw <- shares.mnl.post(beta.draws, choc.standata$X[[1]])

summary(shares.draw)

apply(shares.draw, 2, quantile, probs=c(0.5, 0.025, 0.975))

