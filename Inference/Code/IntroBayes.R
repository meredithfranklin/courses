library(R.utils)
library(lattice) # plotting
library(rstan)
library(R2jags)  # Must have JAGS installed  http://www-ice.iarc.fr/~martyn/software/jags/
library(manipulate) # manipulating plots in RStudio (requires RStudio)
library(MCMCpack)
library(R2WinBUGS)

# Introduction to Bayesian inference. Three equivalent statements:
# What we think about the world after seeing data = What we thought about the world before seeing data x Likelihood of our data under various assumptions about the world
# Pr(world|data) = Pr(world) x Pr(data|world)
# Posterior distribution = Prior distribution x Likelihood
# We use the posterior distribution for inference

# A simple example: Estimating a proportion from dichotomous 0/1 data

# For the prior distribution, choose "hyperparameters" that describe our belief about a quantity of interest (here, p) before seeing data.
# The Beta distribution is "conjugate" to the binomial likelihood, with hyperparameters alpha and beta.

p <- seq(from=0.005, to=0.995, by=0.005)

manipulate(
{plot(p, dbeta(p, alpha.hyper, beta.hyper), 
      col="blue", lwd=2, type="l", las=1, bty="n", 
      ylim=c(0, 8), ylab="density", 
      main="Beta prior distribution")
 polygon(c(p, rev(p)), c(dbeta(p, alpha.hyper, beta.hyper), 
                         rep(0, length(p))), col=rgb(0, 0, 1, 0.2), border=NA)}, 
alpha.hyper=slider(0.1, 10, step=0.1, initial=1), 
beta.hyper=slider(0.1, 10, step=0.1, initial=1))

# Now we observe some data
p.true <- 0.7
N <- 30
y <- rbinom(N, size=1, prob=p.true)
table(y)/N

# Likelihood of the data at each possible value of p

likelihood <- sapply(p, function(p) { prod(p^y * (1-p)^(1-y)) } )
plot(p, likelihood, lwd=2, las=1, bty="n", type="l")

# (To help with visibility)
like.rescale <- 4 * likelihood/max(likelihood)

# To get the posterior, multiply Prior x Likelihood at each value of p
# Or easier: Prior is conjugate, so posterior is Beta distributed with
#  alpha = alpha + k
#  beta = beta + N - k
#  Where N = sample size, k = number of "successes".
#  The prior is most influential when data are sparse.

manipulate(
{plot(p, like.rescale, lwd=2, las=1, bty="n", 
      ylim=c(0,8), type="l", ylab="density", 
      main="Beta prior (blue) x Likelihood (black) = Beta posterior (red)")
 alpha.hyper.post <- alpha.hyper + sum(y)
 beta.hyper.post <- beta.hyper + N - sum(y)
 lines(p, dbeta(p, alpha.hyper, beta.hyper), col="blue", lwd=2)
 polygon(c(p, rev(p)), c(dbeta(p, alpha.hyper, beta.hyper), 
                         rep(0, length(p))), col=rgb(0, 0, 1, 0.2), border=NA)
 lines(p, dbeta(p, alpha.hyper.post, beta.hyper.post), col="red", lwd=2)
 polygon(c(p, rev(p)), c(dbeta(p, alpha.hyper.post, beta.hyper.post), 
                         rep(0, length(p))), col=rgb(1, 0, 0, 0.2), border=NA)
 lines(p, like.rescale, lwd=2)}, 
alpha.hyper=slider(0.1, 10, step=0.1, initial=1), 
beta.hyper=slider(0.1, 10, step=0.1, initial=1))

#  Inferences about the quantity of interest come from the posterior:
#  A best guess (posterior mean or mode), as well as a statement of
#  uncertainty, typically expressed as a "credible interval" or range
#  representing percent area of highest posterior density.


# What if an analytical solution isn't practical?
# Use fast computers and Markov chain Monte Carlo (MCMC) sampling
# to find the shape of the posterior distribution.

# JAGS = "Just Another Gibbs Sampler" to take random draws from the posterior
# Very similar to BUGS = "Bayesian inference Using Gibbs Sampling"

# Specify the Beta-Bernoulli model
jags.bb <- function() {
  for (i in 1:N) {
    y[i] ~ dbin(p, 1)
  }
  # prior on p
  p ~ dbeta(alpha.hyper, beta.hyper)
}

#write.model(jags.bb, "jagsbb.txt")
#file.show("jagsbb.txt")

# Specify hyperpriors and initial values
alpha.hyper <- 1
beta.hyper <- 1
inits <- function() { list(p=runif(1)) }

# Fit the model
jagsfit <- jags.parallel(data=c("y", "N", "alpha.hyper", "beta.hyper"), 
                         inits=inits, 
                         parameters.to.save=c("p"), 
                         model.file=jags.bb, 
                         n.chains=3,
                         n.iter=10000)
jagsfit
plot(jagsfit)
traceplot(jagsfit, mfrow=c(1,1), "p")

## There are additional MCMC diagnostics in the coda package
jagsfit.mcmc <- as.mcmc(jagsfit)
xyplot(jagsfit.mcmc)
densityplot(jagsfit.mcmc)
autocorr.plot(jagsfit.mcmc)
codamenu()

## Posterior distribution of parameters of interest
jags.mod <- jagsfit$BUGSoutput
names(jags.mod)
jags.mod$mean
hist(jags.mod$sims.list$p, breaks=30, col="gray", xlim=c(0,1), main="", freq=F)

# Should match analytical solution:
alpha.hyper.post <- alpha.hyper + sum(y)
beta.hyper.post <- beta.hyper + N - sum(y)
lines(p, dbeta(p, alpha.hyper.post, beta.hyper.post), col="red", lwd=2)

