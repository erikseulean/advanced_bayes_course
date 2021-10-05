# FOR THE UPDATING OF BELIEFS FOR THE PROPORTION OF SUCCESS UNDER BINOMIAL SAMPLING

#**********************************************************************
#**********************************************************************
# Simple R code for Beta prior-Binomial likelihood Bayesian analysis
#**********************************************************************
prior_to_posterior<- function(priormean,priora,n,x){

# priormean is the mean of the Beta prior
# priora is the first parameter of the Beta prior (it will define the variance of the prior Beta)
# n is the number of trials
# x is the number of successes

x11(w=10,h=10)# open a new window: it also works with windows() (or quartz() if you are a Mac user)        
par(mfrow=c(2,1))

priorb=(priora-priora*priormean)/priormean

priorsample=rbeta(50000,priora,priorb)
sortsample=sort(priorsample)
left95interval=sortsample[250]
right95interval=sortsample[9750]

plot(density(priorsample), main = "Prior Beta distribution of p", 
     xlab=paste("Prior is Beta(", priora,", ",round(priorb,3), ")"), 
     sub=paste("Prior mean=",priormean, 
     " , 95% cred. int.=(",round(left95interval,3),",",round(right95interval,3),")"))

posteriora=priora+x
posteriorb=priorb+n-x
posteriorsample=rbeta(50000,posteriora,posteriorb)
posteriormean=round(mean(posteriorsample),3)
sortsample=sort(posteriorsample)
left95interval=sortsample[250]
right95interval=sortsample[9750]

plot(density(posteriorsample), main = "Posterior Beta distribution of p", 
     xlab=paste("Posterior is Beta(", round(posteriora,3),", ",round(posteriorb,3), ")"), 
     sub=paste("Posterior mean=",round(posteriormean,3), 
     " , 95% cred. int.=(",round(left95interval,3),",",round(right95interval,3),")"))
}

# try different values for the function arguments to see how the Bayesian updating works.
# For example, priormean=0.3, priora=1, n=70, x=34 
# i.e. 		prior_to_posterior(0.3,1,70,34)

# Note: You can use the code above to investigate which Beta distribution could match your 
#   prior beliefs about a proportion p. Just input n=0 and x=0 and try different Betas!
