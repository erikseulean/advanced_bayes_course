#writeLines('PATH="${RTOOLS40_HOME}\\usr\\bin;${PATH}"', con = "~/.Renviron")
library(nimble)
library(igraph)
library(coda)
library(R6)

# Specify the statistical model
ratsCode <- nimbleCode({
  
  # Specify the likelihood:
  for(i in 1 : N) {
    for(j in 1 : T) {
      Y[i,j] ~ dnorm(mu[i,j],tau)
      mu[i,j] <- alpha + beta*(x[j] - mean(x[]))/sd(x[])
    }		
  }
  
  # Prior specification:
  alpha ~ dnorm(0,0.00001)
  beta ~ dnorm(0,0.00001)
  tau ~ dgamma(0.001,0.001)
  sigma2 <- 1 / tau
  
})

# Constant values in the model
ratsConsts <- list(N = 30, T = 5)

# Data values
ratsData <- list(x = c(8.0, 15.0, 22.0, 29.0, 36.0),	
                      Y = t(matrix( c(151, 199, 246, 283, 320, # reading data in matrix form
                                    145, 199, 249, 293, 354,   
                                    147, 214, 263, 312, 328,
                                    155, 200, 237, 272, 297,
                                    135, 188, 230, 280, 323,
                                    159, 210, 252, 298, 331,
                                    141, 189, 231, 275, 305,
                                    159, 201, 248, 297, 338,
                                    177, 236, 285, 350, 376,
                                    134, 182, 220, 260, 296,
                                    160, 208, 261, 313, 352,
                                    143, 188, 220, 273, 314,
                                    154, 200, 244, 289, 325,
                                    171, 221, 270, 326, 358,
                                    163, 216, 242, 281, 312,
                                    160, 207, 248, 288, 324,
                                    142, 187, 234, 280, 316,
                                    156, 203, 243, 283, 317,
                                    157, 212, 259, 307, 336,
                                    152, 203, 246, 286, 321,
                                    154, 205, 253, 298, 334,
                                    139, 190, 225, 267, 302,
                                    146, 191, 229, 272, 302,
                                    157, 211, 250, 285, 323,
                                    132, 185, 237, 286, 331,
                                    160, 207, 257, 303, 345,
                                    169, 216, 261, 295, 333,
                                    157, 205, 248, 289, 316,
                                    137, 180, 219, 258, 291,
                                    153, 200, 244, 286, 324),
                        5,30)))
# one set of initial values before building the model                 
ratsInits <- list(alpha = 1, beta = 1, tau = 1)
# to build the model
rats <- nimbleModel(code = ratsCode, name = "rats", constants = ratsConsts,
                    data = ratsData, inits = ratsInits)
rats$getNodeNames()
rats$Y
# To compile the model
Crats <- compileNimble(rats)

# set up the monitored quantities. Default is all of the random quantities
ratsConf <- configureMCMC(rats, monitors = c('alpha','beta', 'tau', 'sigma2'), print = TRUE) 
# build the MCMC algorithm
ratsMCMC <- buildMCMC(ratsConf)
# compile the MCMC chain 
CratsMCMC <- compileNimble(ratsMCMC, project = rats)


####### POSTERIOR SAMPLES IN CODA FORMAT TO GET MORE EASILY PLOTS AND DIAGNOSTICS
set.seed(10)
ratsInits <- list(list(alpha = 1, beta = 1, tau = 1), list(alpha = 10, beta = 10, tau = 10))
posterior <- runMCMC(CratsMCMC, niter = 5000, thin=1, nburnin=1000, 
                     summary = TRUE, WAIC = FALSE, samples = TRUE, nchains=2, 
                     samplesAsCodaMCMC=TRUE, inits = ratsInits) 

combinedchains <- mcmc.list(posterior$samples$chain1, posterior$samples$chain2)
plot(combinedchains)
autocorr.plot(posterior$samples$chain1)
autocorr.plot(posterior$samples$chain2)
gelman.diag(combinedchains)
gelman.plot(combinedchains)
posterior$summary$all.chains
summary(combinedchains)# returns Time-series SE
#
####### POSTERIOR SAMPLES AS R MATRICES   ##############################

par(mfrow = c(1, 3))
plot(posterior$samples$chain1[ , "alpha"], type = "l", xlab = "iteration",
     ylab = expression(alpha))
lines(posterior$samples$chain2[ , "alpha"], type = "l", xlab = "iteration", col=2, 
      ylab = expression(alpha))
plot(posterior$samples$chain1[ , "beta"], type = "l", xlab = "iteration",
     ylab = expression(beta))
lines(posterior$samples$chain2[ , "beta"], type = "l", xlab = "iteration", col=2, 
      ylab = expression(beta))
plot(posterior$samples$chain1[ , "tau"], type = "l", xlab = "iteration",
     ylab = expression(tau))
lines(posterior$samples$chain2[ , "tau"], type = "l", xlab = "iteration", col=2, 
      ylab = expression(tau))
par(mfrow = c(1, 3))
acf(posterior$samples$chain1[, "alpha"]) # plot autocorrelation of alpha sample
acf(posterior$samples$chain1[, "beta"]) # plot autocorrelation of beta sample
acf(posterior$samples$chain1[, "tau"]) # plot autocorrelation of tau sample
par(mfrow = c(1, 3))
acf(posterior$samples$chain2[, "alpha"]) # plot autocorrelation of alpha sample
acf(posterior$samples$chain2[, "beta"]) # plot autocorrelation of beta sample
acf(posterior$samples$chain2[, "tau"]) # plot autocorrelation of tau sample
par(mfrow = c(1, 3))
plot(density(posterior$samples$chain1[ , "alpha"]), type = "l", xlab = expression(alpha),
     ylab = "Posterior density")
plot(density(posterior$samples$chain1[ , "beta"]), type = "l", xlab = expression(beta),
     ylab = "Posterior density")
plot(density(posterior$samples$chain1[ , "tau"]), type = "l", xlab = expression(tau),
     ylab = "Posterior density")
posterior$summary$all.chains
