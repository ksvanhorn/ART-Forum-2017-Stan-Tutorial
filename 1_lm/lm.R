# Fitting a linear model using Stan
# Elea McDonnell Feit, efeit@drexel.edu
# 22 June 2017

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

# Make sure you have the following packages installed (only need to to once)
if (FALSE) {
  install.packages("MASS")
  install.packages("rstan")
  install.packages("bayesplot")
  install.pacakges("MCMCpack")
}

# ========= Load and inspect computer conjoint study ==========
cc.df <- read.csv("ComputerConjointData.csv")
head(cc.df)

summary(cc.df)

plot(LIKE ~ Price, data = cc.df, ylab="Rating (0-10)") 

# ========= Specify Stan model =========
lm.stan <- "
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
"

# ======== Prepare data for Stan =========
y <- cc.df[!is.na(cc.df$LIKE),"LIKE"]
str(y)

# wrong
X <- cc.df[!is.na(cc.df$LIKE), 4:16] # attribute cols
head(X)
  
X <- model.matrix(~ HotLine + RAM + Screen + CPU + HardDisk + CD + Cache + Color + 
                    Channel + Warranty + Software + Guarantee + Price, 
                  data=cc.df[!is.na(cc.df$LIKE),])
X <- X[,2:ncol(X)] # remove intercept
head(X)

cc.standata <- list(N=length(y), K=ncol(X), y=y, X=X)
str(cc.standata)
rm(y, X)

# ========= Run Stan Model ==========
library(rstan)  # load interface to Stan
m.stan <- stan(model_code = lm.stan, data = cc.standata, iter = 1000, chains=1)

str(m.stan)

mydraws <- extract(m.stan, pars=c("beta"))
head(mydraws$beta)

summary(m.stan)$summary

cc.lab <- c("Intercept", dimnames(cc.standata$X)[[2]], "sigma", "lp__")
data.frame(label=cc.lab, summary(m.stan)$summary)

mean(mydraws$beta[,6]>0)
mean(mydraws$beta[,8]>0)

plot(m.stan, pars=c("beta0", "beta"))

plot(m.stan, plotfun="trace", pars=c("beta0", "beta", "sigma"))

plot(m.stan, plotfun="hist", pars=c("beta0", "beta", "sigma"))

get_stanmodel(m.stan)


