#Example Gibbs sampler, for X~N(mu,sigma^2) when mu and sigma are unknown
# normal(phi,tau^2) prior on mu and IG on sigma^2(alpha,beta)

#Uncomment this line if you want to get the same answer each time you run it
#set.seed(1293820)

#Priors
phi<-0
tau2<-1
alpha<-0.1
beta<-0.01

#Simulate some data
n<-10
x<-rnorm(n,0,1)
xbar<-mean(x)

#Number of samples to take in the chain
T<-100
#Create space to store the Markov chain in 
mu<-sigma2<-numeric(T)

#Set starting values
mu[1]<-4
sigma2[1]<-0.1

#Run the Gibbs sampler 
for(t in 1:(T-1)){
  mu[t+1]<-rnorm(1,(tau2*n*xbar + sigma2[t]*phi)/(tau2*n + sigma2[t]),
            sqrt(sigma2[t]*tau2/(tau2*n+sigma2[t])))
  sigma2[t+1]<-1/rgamma(1,shape=n/2 + alpha,
                          rate=1/2*sum((x-mu[t+1])^2) + beta)
}

#Produce trace plots of the outputs
par(mfrow=c(1,2))
plot(1:T,mu,xlab="Iteration",ylab="mu",type="l")
plot(1:T,sigma2,xlab="Iteration",ylab="sigma2",type="l")
par(mfrow=c(1,1))

#Produce density plots of the outputs
#(Note - you'd want to remove some initial samples as burn-in
# before using the samples for inference)
par(mfrow=c(1,2))
plot(density(mu))
plot(density(sigma2))
par(mfrow=c(1,1))

#Produce bivariate plot of samples
plot(mu,sigma2,type="p")

posteriormean_mu=mean(mu)
sortsample=sort(mu)
left95interval=(sortsample[2]+sortsample[3])/2
right95interval=(sortsample[97]+sortsample[98])/2
cat("posterior mean for mu=",posteriormean_mu, sep=c(""))
cat(" ", sep=c("\n"))
cat("95% credible interval=(",left95interval,",",right95interval,")", sep=c(""))

posteriormean_sigma2=mean(sigma2)
sortsample=sort(sigma2)
left95interval=(sortsample[2]+sortsample[3])/2
right95interval=(sortsample[97]+sortsample[98])/2
cat("posterior mean for sigma^2=",posteriormean_sigma2, sep=c(""))
cat(" ", sep=c("\n"))
cat("95% credible interval=(",left95interval,",",right95interval,")", sep=c(""))

# Of course, if you were to write your own code not as a very basic example, 
# but to properly perform some analysis,
# you would have to allow for a burn-in, and to also obtain additional inferences and 
# diagnostic plots 
