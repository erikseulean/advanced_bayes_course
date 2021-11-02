#writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
library(nimble)
library(igraph)
library(coda)
library(R6)
#
# HEIGHTS DATA (lecture 15)
#
# Generate the data #####
set.seed(15)
HeightData=c(rnorm(550,185,8),rnorm(700,170,8))# heights in cm.
# #######################
#
# Specify the statistical model
MixNormCode <- nimbleCode({
  
  # Specify the likelihood:
  for (i in 1:N){
    x[i] ~ dnorm(mu[t[i]],s)
    t[i] <- z[i]+1
    z[i] ~ dbin(pi,1)# auxiliary variable
  }
  # Prior specification:
  mu[1] ~ dnorm(phi,lambda)
  mu[2] ~ dnorm(phi,lambda)
  pi ~ dbeta(1,1)
  
})
#
# Values for some constants in the model
MixNormConsts <- list(N = 700+550, s=1/8^2, phi=177, lambda=1/15^2 ) 
#
# The data values
MixNormData <- list(x=HeightData)

# one set of initial values before building the model                 
MixNormInits <- list(pi=0.5,mu=c(177,177),z=rbinom(1250,1,0.5)) 

# to build the model
MixNorm <- nimbleModel(code = MixNormCode, name = "MixNorm", constants = MixNormConsts,
                        data = MixNormData, inits= MixNormInits)
#
# To compile the model
CMixNorm <- compileNimble(MixNorm)

# set up the monitored quantities. Default is all of the random quantities
MixNormConf <- configureMCMC(MixNorm, monitors = c('pi','mu','z'), print = TRUE)
# binary sampler= gibbs sampler for binary-valued obs
#
# build the MCMC algorithm
MixNormMCMC <- buildMCMC(MixNormConf)
# compile the MCMC chain 
CMixNormMCMC <- compileNimble(MixNormMCMC, project = MixNorm)
#
set.seed(15)
#set.seed(100)
MixNormInits <- list(list(pi = 0.5, mu=c(20,500),z=rbinom(1250,1,0.5)),
                      list(pi=0.5, mu=c(500,20),z=rbinom(1250,1,0.5)))
#
posterior <- runMCMC(CMixNormMCMC, niter = 10000, thin=1, nburnin=1, 
                     summary = TRUE, WAIC = FALSE, samples = TRUE, nchains=2, 
                     samplesAsCodaMCMC=TRUE, inits = MixNormInits) 
#
#
combinedchains <- mcmc.list(posterior$samples$chain1[,c("pi","mu[1]","mu[2]")],
                            posterior$samples$chain2[,c("pi","mu[1]","mu[2]")])
windows()# or quartz() on a MacOS
plot(combinedchains)
windows();gelman.plot(combinedchains)# 
autocorr.plot(posterior$samples$chain1[,c("pi","mu[1]","mu[2]")])
autocorr.plot(posterior$samples$chain2[,c("pi","mu[1]","mu[2]")])
effectiveSize(combinedchains)
posterior$summary$all.chains
summary(combinedchains)# returns Time-series SE
#####################################
# However...
#####################################

set.seed(100)# by just changing the seed..
#
MixNormInits <- list(list(pi = 0.5, mu=c(20,500),z=rbinom(1250,1,0.5)),
                     list(pi=0.5, mu=c(500,20),z=rbinom(1250,1,0.5)))
#
posterior <- runMCMC(CMixNormMCMC, niter = 10000, thin=1, nburnin=1, 
                     summary = TRUE, WAIC = FALSE, samples = TRUE, nchains=2, 
                     samplesAsCodaMCMC=TRUE, inits = MixNormInits) 
#
#
combinedchains <- mcmc.list(posterior$samples$chain1[,c("pi","mu[1]","mu[2]")],
                            posterior$samples$chain2[,c("pi","mu[1]","mu[2]")])
windows()# or quartz() on a MacOS
plot(combinedchains)
#
# The posterior distributions for mu and pi have actually two modes (we saw just
# one mode in the first run -line 60- only by chance). 
# The marginal distributions for mu[1] and mu[2] are identical unless we relabel
# the samples (manually or using postprocessing algorithms) or put constraints 
# to the priors to make the separate mixture components distinguishible (not
# covered in this course) 