# Example of Metropolis sampling from section 2.3 of the lecture notes.
# Sampling from the standard Normal distribution. 

T<-500
#Variability of the normal distribution proposal
std<-sqrt(1) #0.1,1,100 

#store samples
theta<-numeric(T)
#store number of acceptances
n.accept<-0

#starting value

theta[1]<-0 #Set start value at 0 when exploring proposals
for (i in 1:(T-1)){
  phi<-rnorm(1,theta[i],std)  
  alpha<-min(c(1,exp(-0.5*(phi^2-theta[i]^2))))
  if(runif(1)<alpha) {
    theta[i+1]<-phi
    n.accept<-n.accept+1
  } else {
    theta[i+1]<-theta[i]
  }
}
p.accept<-n.accept/T
p.accept
#
#
#Calculate effective sample size
library(coda)
chain<-mcmc(theta)
ess<-effectiveSize(chain)
ess

#Plot using home-grown code
#
windows()
par(mfrow=c(2,2))
plot(1:T,theta,type="l",main=paste("p.accept=",format(p.accept,digits=2),", ESS=",format(ess,digits=1),sep=""))
qqnorm(theta)
hist(theta)
acf(theta)
summary(theta)
sd(theta)

# For an analysis (rather than an illustration of the MH sampler properties)
#  remember to discard the burn-in sample. 
