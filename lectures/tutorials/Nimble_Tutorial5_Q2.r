#writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
library(nimble)
library(igraph)
library(coda)
library(R6)

# Tutorial 3 Question 2

################################################################
##### BINOMIAL LIKELIHOOD - BETA PRIOR
################################################################

# Specify the statistical model
BinBetaCode <- nimbleCode({
  
  #likelihood
  x~dbin(p,n)
  
  #priors
  p~dbeta(0.5,0.5)
  
})

# The data values
BinBetaData <- list(n=70,x=34)

# one set of initial values before building the model                 
BinBetaInits <- list(p=0.1)   

# to build the model
BinBeta <- nimbleModel(code = BinBetaCode, name = "BinBeta",
                    data = BinBetaData, inits<-BinBetaInits)

# To compile the model
CBinBeta <- compileNimble(BinBeta)

# set up the monitored quantities. Default is all of the random quantities
BinBetaConf <- configureMCMC(BinBeta, monitors = c('p'), print = TRUE) 
# build the MCMC algorithm
BinBetaMCMC <- buildMCMC(BinBetaConf)
# compile the MCMC chain 
CBinBetaMCMC <- compileNimble(BinBetaMCMC, project = BinBeta)

####################################################################################
####### POSTERIOR SAMPLES IN CODA FORMAT TO GET MORE EASILY PLOTS AND DIAGNOSTICS  #
####################################################################################
set.seed(10)
BinBetaInits <- list(list(p=1), list(p=0.9))
posterior <- runMCMC(CBinBetaMCMC, niter = 20000, thin=1, nburnin=2000, 
                     summary = TRUE, WAIC = FALSE, samples = TRUE, nchains=2, 
                     samplesAsCodaMCMC=TRUE, inits = BinBetaInits) 

combinedchains <- mcmc.list(posterior$samples$chain1, posterior$samples$chain2)
plot(combinedchains)
autocorr.plot(posterior$samples$chain1)
autocorr.plot(posterior$samples$chain2)
gelman.diag(combinedchains)
gelman.plot(combinedchains)
posterior$summary$all.chains


###############################################################################
###############################################################################
################################################################
##### NORMAL LIKELIHOOD - NORMAL AND IG PRIORS
################################################################

# Specify the statistical model
NormalCode <- nimbleCode({
  
  #likelihood
  for(i in 1:n){
    x[i]~dnorm(mu,precA)
  }
  sigma2<-1/precA
  #priors
  mu~dnorm(phi,precB)
  tau2<-1/precB
  precA~dgamma(alpha,beta)
  
})


# Constant values in the model
NormalConsts <- list(n = 8) # Need to define this as a constant as it is the number of loops

# The data values
NormalData <- list(phi=0,precB=0.01,alpha=0.001,beta=0.001,x=c(36, 67, 44, 39, 56, 65, 43, 49))

# one set of initial values before building the model                 
NormalInits <- list(mu=1, precA=2)  # Nimble can do this by sampling from the priors. 

# to build the model
Normal <- nimbleModel(code = NormalCode, name = "Normal", constants = NormalConsts,
                       data = NormalData, inits<-NormalInits)

# To compile the model
CNormal <- compileNimble(Normal)

# set up the monitored quantities. Default is all of the random quantities
NormalConf <- configureMCMC(Normal, monitors = c('mu','sigma2'), print = TRUE) 
# build the MCMC algorithm
NormalMCMC <- buildMCMC(NormalConf)
# compile the MCMC chain 
CNormalMCMC <- compileNimble(NormalMCMC, project = Normal)

####################################################################################
####### POSTERIOR SAMPLES IN CODA FORMAT TO GET MORE EASILY PLOTS AND DIAGNOSTICS  #
####################################################################################
set.seed(10)
BinBetaInits <- list(list(mu=1, precA=2), list(mu=-5, precA=4))
posterior <- runMCMC(CNormalMCMC, niter = 20000, thin=1, nburnin=2000, 
                     summary = TRUE, WAIC = FALSE, samples = TRUE, nchains=2, 
                     samplesAsCodaMCMC=TRUE, inits = NormalInits) 

combinedchains <- mcmc.list(posterior$samples$chain1, posterior$samples$chain2)
windows();plot(combinedchains)
autocorr.plot(posterior$samples$chain1)
autocorr.plot(posterior$samples$chain2)
gelman.diag(combinedchains)
gelman.plot(combinedchains)
posterior$summary$all.chains

