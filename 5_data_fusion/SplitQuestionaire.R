# Data Fusion / split quetionaire example in Stan
# Elea McDonnell Feit, eleafeit@gmail.com
# 11 March 2016

# From Feit, Elea McDonnel and Eric T. Bradlow (2017) Fusion Modeling
# Invited chapter in Handbook of Market Research by 
# Christian Homburg, Martin Klarmann and Arnd Vomberg

# Copyright 2017, Elea McDonnell Feit 
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

rm(list=ls())
setwd("~/GitHub/artforum2017/5_data_fusion")
list.files()

library(MASS)
library(coda)
library(beanplot)
library(rstan)

# Functions ===============================================================
data.mvn.split <- function(K1=2, K2=2, Kb=3, N1=100, N2=100,
                           mu=rep(0, K1+K2+Kb), Sigma=diag(1, K1+K2+Kb)) 
{
  y <- mvrnorm(n=N1+N2, mu=mu, Sigma=Sigma)
  list(data=list(K1=K1, K2=K2, Kb=Kb, N1=N1, N2=N2, 
                 y1=as.matrix(y[1:N1, 1:K1], col=K1), 
                 y2=as.matrix(y[N1+1:N2, K1+1:K2], col=K2),
                 yb=as.matrix(y[,K1+K2+1:Kb], col=Kb)), 
       true=list(mu=mu, Sigma=Sigma, 
                 y1mis=y[1:N1, K1+1:K2], 
                 y2mis=y[N1+1:N2, 1:K1]))
}   

data.mvp.split <- function(K1=2, K2=2, Kb=3, N1=100, N2=100,
                           mu=rep(0, K1+K2+Kb), Sigma=diag(1, K1+K2+Kb)) 
{
  z <- mvrnorm(n=N1+N2, mu=mu, Sigma=Sigma)
  y <- z
  y[y>0] <- 1
  y[y<0] <- 0
  y1mis <- y[1:N1, K1+1:K2]
  y2mis <- y[N1+1:N2, 1:K1]
  y[1:N1, K1+1:K2] <- NA
  y[N1+1:N2, 1:K1] <- NA
  true=list(mu=mu, Sigma=Sigma, z=z, y=y, y1mis=y1mis, y2mis=y2mis)
  y[is.na(y)] <- 0
  data=list(K1=K1, K2=K2, Kb=Kb, N1=N1, N2=N2, y=y) 
  list(data=data, true=true)
} 

plot.post.density <- function(m.stan, pars, true, prefix=NULL){
  for (i in 1:length(pars)) {
    draws <- As.mcmc.list(m.stan, pars=pars[i])
    if (!is.null(prefix)) {
      filename <- paste(prefix, "Post",  pars[i], ".png", sep="")
      png(filename=filename, width=600, height=400)
    }
    beanplot(data.frame(draws[[1]]), 
             horizontal=TRUE, las=1, what=c(0, 1, 1, 0), side="second",
             main=paste("Posterior Density of", pars[[i]]))
    if (!is.null(prefix)) dev.off()
  }
}

plot.true.v.est <- function(m.stan, pars, true, prefix=NULL){
  for (i in 1:length(pars)) {
    draws <- As.mcmc.list(m.stan, pars=pars[i]) 
    est <- summary(draws)
    if (!is.null(prefix)) {
      filename <- paste(prefix, "TrueVEst", pars[i], ".png", sep="")
      png(filename=filename, width=600, height=400)
    }
    plot(true[[i]], est$quantiles[,3], col="blue", 
         xlab=paste("True", pars[i]), 
         ylab=paste("Estiamted", pars[i], "(posterior median)"))
    abline(a=0, b=1)
    arrows(true[[i]], est$quantiles[,3], true[[i]], est$quantiles[,1], 
           col="gray90", length=0)
    arrows(true[[i]], est$quantiles[,3], true[[i]], est$quantiles[,5], 
           col="gray90", length=0)
    points(true[[i]], est$quantiles[,3], col="blue")
    if (!is.null(prefix)) dev.off()
  }
}

# Example 1a: MVN ================================================================
# Generate synthetic data
set.seed(20030601)
Sigma <- matrix(c(1, 0.3, -0.2, 0.7, 0.3, 1, -0.6, 0.4, 
                  -0.2, -0.6, 1, 0.1, 0.7, 0.4, 0.1, 1), nrow=4)
d1 <- data.mvn.split(K1=1, K2=1, Kb=2, N1=100, N2=100, mu=rep(0,4), Sigma=Sigma)
str(d1$data)
# Call to Stan to generate posterior draws
m1 <- stan(file="SplitQuestionaire.stan", data=d1$data, 
           iter=10000, warmup=2000, chains=1, seed=12)
# Summaries of posterior draws for population-level parameters
summary(m1, par=c("mu"))  
summary(m1, par=c("tau"))
summary(m1, par=c("Omega"))
plot.post.density(m1, pars=c("mu", "tau"),
                  true=list(d1$true$mu, sqrt(diag(d1$true$Sigma)), 
                            cov2cor(d1$true$Sigma)))
draws <- As.mcmc.list(m1, pars=c("Omega"))
beanplot(data.frame(draws[[1]][,c(2:4, 7:8, 12)]), 
         horizontal=TRUE, las=1, what=c(0, 1, 1, 0), side="second",
         main=paste("Posterior Density of Omega (correlations)", log=""), 
         cex.axis=0.5)
# Summaries of posterior draws for missing data
summary(extract(m1, par=c("y1mis"))$y1mis[,3,]) 
plot(density(extract(m1, par=c("y1mis"))$y1mis[,3,]), 
     main="Posterior of Unobserved y_1", xlab="y_1")

summary(m1, par=c("y"))  # posteriors of observed data place a point mass at the observed value
plot.true.v.est(m1, pars=c("y1mis", "y2mis"), 
                true=list(d1$true$y1mis, d1$true$y2mis))

# Example 1b: MVN with zero correlations =======================================
# Generate synthetic data
set.seed(20030601)
Sigma <- matrix(0, nrow=4, ncol=4)
diag(Sigma) <- 1
# Call to Stan to generate posterior draws
d2 <- data.mvn.split(K1=1, K2=1, Kb=2, N1=100, N2=100, mu=rep(0,4), Sigma=Sigma)
m2 <- stan(file="SplitQuestionaire.stan", data=d2$data, 
           iter=10000, warmup=2000, chains=1, seed=12)
# Summarize posteriors of population-level parameters
summary(m2, par=c("mu"))  
summary(m2, par=c("tau"))
summary(m2, par=c("Omega"))
plot.post.density(m2, pars=c("mu", "tau"), 
                  true=list(d1$true$mu, sqrt(diag(d1$true$Sigma)), 
                            cov2cor(d1$true$Sigma)))
draws <- As.mcmc.list(m2, pars=c("Omega"))
beanplot(data.frame(draws[[1]][,c(2:4, 7:8, 12)]), 
         horizontal=TRUE, las=1, what=c(0, 1, 1, 0), side="second",
         main=paste("Posterior Density of Omega", log=""), cex.axis=0.5)

# Summaries of posterior draws for missing data
plot.true.v.est(m2, pars=c("y1mis", "y2mis"), 
                true=list(d2$true$y1mis, d2$true$y2mis))

# Example 1c: MVN with strong positive  correlations ===========================
# Generate synthetic data
set.seed(20030601)
Sigma <- matrix(0.9, nrow=4, ncol=4)
diag(Sigma) <- 1
# Call to Stan to generate posterior draws
d3 <- data.mvn.split(K1=1, K2=1, Kb=2, N1=100, N2=100, mu=rep(0,4), Sigma=Sigma)
m3 <- stan(file="SplitQuestionaire.stan", data=d3$data, 
           iter=10000, warmup=2000, chains=1, seed=12)
# Summaries of population-level parameters
summary(m3, par=c("mu"))  
summary(m3, par=c("tau"))
summary(m3, par=c("Omega"))
plot.post.density(m3, pars=c("mu", "tau"), 
                  true=list(d1$true$mu, sqrt(diag(d1$true$Sigma))))
draws <- As.mcmc.list(m3, pars=c("Omega"))
beanplot(data.frame(draws[[1]][,c(2:4, 7:8, 12)]), 
         horizontal=TRUE, las=1, what=c(0, 1, 1, 0), side="second",
         main=paste("Posterior Density of Omega", log=""))

# Summaries of posterior draws for missing data
plot.true.v.est(m3, pars=c("y1mis", "y2mis"), 
                true=list(d3$true$y1mis, d3$true$y2mis))

# Example 2: MVP =================================================================
# Generate synthetic data
set.seed(20030601)
Sigma <- matrix(c(1, 0.3, -0.2, 0.7, 0.3, 1, -0.6, 0.4, 
                  -0.2, -0.6, 1, 0.1, 0.7, 0.4, 0.1, 1), nrow=4)
d1 <- data.mvp.split(K1=1, K2=1, Kb=2, N1=100, N2=100, mu=rep(0,4), Sigma=Sigma)
# Call to Stan to generate posterior draws
m1 <- stan(file="SplitQuestionaireMVP.stan", data=d1$data, 
           iter=10000, warmup=2000, chains=1, seed=35)
# Check that the z's are consistent with their observed values
summary(m1, par=c("z"))
d1$true$y
# Summaries of posteriors of population-level parameters
summary(m1, par=c("mu", "Omega"))
plot.post.density(m1, pars=c("mu"), true=list(d1$true$mu))
draws <- As.mcmc.list(m1, pars=c("Omega"))
beanplot(data.frame(draws[[1]][,c(2:4, 7:8, 12)]), 
         horizontal=TRUE, las=1, what=c(0, 1, 1, 0), side="second",
         main=paste("Posterior Density of Omega", log=""))

# Summarize posterior of missing data
y1mis.est <- summary(m1, par=c("y1mis"))$summary
y1mis.est[1,]
xtabs(~y1mis.est[,"50%"]+d1$true$y1mis)
y2mis.est <- summary(m1, par=c("y1mis"))$summary
xtabs(~y2mis.est[,"50%"]+d1$true$y2mis)
z.est <- data.frame(z.true=as.vector(t(d1$true$z)), y=as.vector(t(d1$true$y)), 
                    z.postmed=summary(m1, pars=c("z"))$summary[,"50%"])
plot(z.est[,c(1,3)], xlab="True Latent Variable", ylab="Posterior Mean of Latent Variable")
points(z.est[is.na(z.est$y), c(1,3)], col="red")
abline(h=0, v=0)


Sigma <- matrix(0, nrow=4, ncol=4)
diag(Sigma) <- 1
d2 <- data.mvp.split(K1=1, K2=1, Kb=2, N1=100, N2=100, mu=rep(0,4), Sigma=Sigma)
m2 <- stan(file="SplitQuestionaireMVP.stan", data=d2$data, 
           iter=10000, warmup=2000, chains=1, seed=35)
print(m2, par=c("mu", "Omega"))
plot.post.density(m2, pars=c("mu"), true=list(d2$true$mu))
draws <- As.mcmc.list(m2, pars=c("Omega"))
beanplot(data.frame(draws[[1]][,c(2:4, 7:8, 12)]), 
         horizontal=TRUE, las=1, what=c(0, 1, 1, 0), side="second",
         main=paste("Posterior Density of Omega", log=""))

y1mis.est <- summary(m2, par=c("y1mis"))$summary
xtabs(~y1mis.est[,"50%"]+d2$true$y1mis)
y2mis.est <- summary(m2, par=c("y1mis"))$summary
xtabs(~y2mis.est[,"50%"]+d2$true$y2mis)
z.est <- data.frame(z.true=as.vector(t(d2$true$z)), y=as.vector(t(d2$true$y)), 
                    z.postmed=summary(m2, pars=c("z"))$summary[,"50%"])
plot(z.est[,c(1,3)], xlab="True Latent Variable", ylab="Posterior Mean of Latent Variable")
points(z.est[is.na(z.est$y), c(1,3)], col="red")
abline(h=0, v=0)



Sigma <- matrix(0.9, nrow=4, ncol=4)
diag(Sigma) <- 1
d3 <- data.mvp.split(K1=1, K2=1, Kb=2, N1=100, N2=100, mu=rep(0,4), Sigma=Sigma)
m3 <- stan(file="SplitQuestionaireMVP.stan", data=d3$data, 
           iter=10000, warmup=2000, chains=1, seed=12)
print(m3, par=c("mu", "Omega"))
plot.post.density(m3, pars=c("mu"), true=list(d2$true$mu))
draws <- As.mcmc.list(m3, pars=c("Omega"))
beanplot(data.frame(draws[[1]][,c(2:4, 7:8, 12)]), 
         horizontal=TRUE, las=1, what=c(0, 1, 1, 0), side="second",
         main=paste("Posterior Density of Omega", log=""))

y1mis.est <- summary(m3, par=c("y1mis"))$summary
xtabs(~y1mis.est[,"50%"]+d3$true$y1mis)
y2mis.est <- summary(m3, par=c("y1mis"))$summary
xtabs(~y2mis.est[,"50%"]+d3$true$y2mis)
z.est <- data.frame(z.true=as.vector(t(d3$true$z)), y=as.vector(t(d3$true$y)), 
                    z.postmed=summary(m3, pars=c("z"))$summary[,"50%"])
plot(z.est[,c(1,3)], xlab="True Latent Variable", ylab="Posterior Mean of Latent Variable")
points(z.est[is.na(z.est$y), c(1,3)], col="red")
abline(h=0, v=0)


Sigma <- matrix(0.5, nrow=7, ncol=7)
diag(Sigma) <- 1
d4 <- data.mvp.split(K1=2, K2=2, Kb=3, N1=100, N2=100, mu=rep(0,7), Sigma=Sigma)
m4 <- stan(file="SplitQuestionaireMVP20160312.stan", data=d4$data, 
           iter=10000, warmup=2000, chains=1, seed=12)
print(m4, pars=c("mu", "Omega"))
