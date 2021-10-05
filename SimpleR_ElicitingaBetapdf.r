
#***********************************************************************************************************
# Simple R code (naive?) for eliciting a Beta(a,b) pdf for the probability of success in a Binomial experiment
#***********************************************************************************************************

elicit.beta <- function(alpha=0.05,n.sim.samples=10000){

cat("Provide a value in (0,1) that you expect to be the true probability of success", sep="\n")
priormean<-scan() # NOTE: press ENTER again when 2: appears

cat(paste("Provide the left bound of a single interval within (0,1) you believe is ",(1-alpha)*100,"% likely to contain the true probability",sep=""), sep="\n")
leftpriorint<-scan()

cat(paste("Provide the right bound of a single interval within (0,1) you believe is ",(1-alpha)*100,"% likely to contain the true probability", sep=""), sep="\n")
rightpriorint<-scan()

# loop to find the parameters of the Beta distribution 
priora<-0
e<-10
while(abs(e)>0.0001) {
priora<-priora+0.001
priorb=(priora-priora*priormean)/priormean
e<-(pbeta(rightpriorint,priora,priorb)-pbeta(leftpriorint,priora,priorb))-(1-alpha)
if (priora>1000) {
cat("priora>1000. Likely no Beta density exists that satisfies these beliefs", sep="\n")
stop("No Beta density found")} 
}

# examine/check the derived the distribution 
cat("prior distribution is Beta(",priora,",",priorb,")", sep=c(""))
cat("", sep=c("\n"))
cat("prior mean=",priora/(priora+priorb), sep=c("","\n"))
leftint<-qbeta((alpha/2),priora,priorb)
rightint<-qbeta((1-alpha/2),priora,priorb)
cat(paste((1-alpha)*100,"% credible interval =(",leftint,",",rightint,")",sep=""),sep=c("\n"))
priorsample=rbeta(n.sim.samples,priora,priorb)
hist(priorsample)

# check probabilities for elicitation feedback
rprob1<-sum(priorsample>0.2 & priorsample<0.8)/n.sim.samples
intprob1<-pbeta(0.8,priora,priorb)-pbeta(0.2,priora,priorb)
cat("probability within (0.2,0.8)=",rprob1," (using simulation)", sep=c("","\n"))
cat("probability within (0.2,0.8)=",intprob1," (using integration)", sep=c("","\n"))

rprob2<-sum(priorsample>0.9)/n.sim.samples
intprob2<-1-pbeta(0.9,priora,priorb)
cat("probability within (0.9,1)=",rprob2," (using simulation)", sep=c("","\n"))
cat("probability within (0.9,1)=",intprob2," (using integration)", sep=c("","\n"))
}

# The code above can find a Beta distribution that fits lower and upper quartiles by setting
# elicit.beta(0.5,10000)

# There are beliefs for which this code works well, e.g.
# for prior mean=0.5 and prior 50% CI (0.4,0.6) or
# for priormean=0.12, prior 95% CI (0.001,0.5))

# There are other scenarios (e.g. priormean=0.12, 95% CI (0.1,0.5))
# where the code above is not satisfactory. 