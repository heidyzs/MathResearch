library(sn) # sampling from the skew-normal distribution is accomplished via the 'sn' package in R
library(e1071) # needed to calculate sample skewness
library(resample)
library(tidyverse)
library(msm)

#FUNCTIONS#

common.sd = 12.5

#Likelihood
likelihood<-function(mu,component,dat){
  sn_mean<-mu
  cp<-c(mean=sn_mean, s.d.=common.sd,gamma1= .5*(-1)^component)
  dp<-cp2dp(cp,family="SN")
  LB<-psn(0,dp=dp)
  UB<-psn(100,dp=dp)
  dsn(dat, dp = dp)/(UB-LB)
}


#Log Likelihood
log.lik<-function(mu, dat){
  if ((length(dat) == 0)) {
    return(0)
  } else {
    sn_mean<-mu
    cp<-c(mean=mu, s.d.=sig, gamma1= gam)
    dp<-cp2dp(cp,family="SN")
    LB<-psn(0,dp=dp)
    UB<-psn(100,dp=dp)
    sum(dsn(dat,dp=dp,log=T)-log(UB-LB))
  }
}

#Log likelihood for truncated normal distribution
log.lik.tnorm<-function(mu, sig, dat){
  if ((length(dat) == 0)) {
    return(0)
  } else {
    sum(dtnorm(dat, mean = mu, sd = sig, lower = 0, upper = 100, log=T))
  }
}

estBetaParams <- function(mu.scaled, var.scaled) {
  mu<-mu.scaled/100
  var<-var.scaled/100^2
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

#Prior For Mu
prior<-function(mu,muprior,sdprior=20){
  varprior<-sdprior^2
  convert<-estBetaParams(muprior,sdprior)
  alpha<-convert$alpha
  beta<-convert$beta
  dbeta(mu/100,alpha,beta,log=TRUE)
}

#Likelihood function used for MLE

lik2<-function(parvec, dat){
  if ((length(dat) == 0)) {
    return(0)
  } else {
    cp<-c(mean=parvec[1], s.d.=parvec[2], gamma1= parvec[3])
    dp<-cp2dp(cp,family="SN")
    LB<-psn(0,dp=dp)
    UB<-psn(100,dp=dp)
    sum(dsn(dat,dp[1], dp[2], dp[3], log = T)-log(UB-LB)) #need to simulate data before using function
  }
}



###################################################################################

#Start Simulation
set.seed(8)

N.bootstrap <- 1000
n <- 25


mean1.range <- 40 #seq(from = 45, to = 20, by = -5)
mean2.range <- 60 #seq(from = 55, to = 80, by = 5 )
N.trials <-100

results <- data.frame(matrix(0, ncol = 24, nrow = N.trials*length(mean1.range))) #*length(mean1.range)
x <- c("trial", "mean1", "mean2", "median.sampsmu1", "lb.mu1", "ub.mu1", "width.mu1", "cover.mu1",
       "median.samps.mu2", "lb.mu2", "ub.mu2", "width.mu2", "cover.mu2", "perc.correct", "lb.tdist",
       "ub.tdist", "width.tdist","cover.tdist.mean1", "cover.tdist.mean2",
       "lb.bca.exp", "ub.bca.exp", "width.bca.exp", "cover.bca.exp.mean1", "cover.bca.exp.mean2")
colnames(results) <- x

common.sd <- 12.5

start_time <- Sys.time()
count <- 1


for (j in 1:length(mean1.range)) {
  
  #Distribution 1
  sn_mean1 <- mean1.range[j]
  sn_stdev1 <- common.sd 
  sn_skew1 <- 0.5  
  cp1 <- c(mean=sn_mean1, s.d.=sn_stdev1, gamma1=sn_skew1)
  dp1 <- cp2dp(cp1, family="SN") # the cp2dp function converts the centered parameterization to the direct parameterization 
  LB.for.prob1<-psn(0, dp=dp1) # finds the probability associated with the minimum SUS
  UB.for.prob1<-psn(100, dp=dp1) # finds the probability associated with the maximum SUS
  
  #Distribution 2
  sn_mean2 <- mean2.range[j]
  sn_stdev2 <- common.sd
  sn_skew2 <- -0.5
  cp2 <- c(mean=sn_mean2, s.d.=sn_stdev2, gamma1=sn_skew2)
  dp2 <- cp2dp(cp2, family="SN") # the cp2dp function converts the centered parameterization to the direct parameterization 
  LB.for.prob2<-psn(0, dp=dp2) # finds the probability associated with the minimum SUS
  UB.for.prob2<-psn(100, dp=dp2) # finds the probability associated with the maximum SUS
  
  func1<-function(x){
    x*dsn(x,dp1[1],dp1[2],dp1[3])*(1/(UB.for.prob1-LB.for.prob1))
  }
  true.mean1<-integrate(func1,0,100)[1]$value
  
  func2<-function(x){
    x*dsn(x,dp2[1],dp2[2],dp2[3])*(1/(UB.for.prob2-LB.for.prob2))
  }
  true.mean2<-integrate(func2,0,100)[1]$value
  
  for (l in 1:N.trials) {
    
    w <- rbinom(n, 1, 0.5) #samples to create a probability for each cluster
    data = vector(length = n)
    
    a = 1
    #Loop to generate data from a mixture distribution
    for(i in w){
      if (i == 1) {
        data[a] = qsn(runif(n=1,min=LB.for.prob1,max=UB.for.prob1), dp=dp1, solver="RFB")
        a <- a +1
      }
      else{
        data[a] = qsn(runif(n=1,min=LB.for.prob2,max=UB.for.prob2), dp=dp2, solver="RFB")
        a <- a + 1
      }
    }
    
    
    #Begin Sampling
    b <- bootstrap(data, mean(data), R = N.bootstrap)  
    # generate 95% CIs from the bootstrap sample
    lb_tdist <- mean(data) - qt(0.975, df = length(data)-1)*sd(data)/sqrt(length(data))
    ub_tdist <- mean(data) + qt(0.975, df = length(data)-1)*sd(data)/sqrt(length(data))
    ci_bca_expanded <- CI.bca(b, probs = c(0.025, 0.975), expand = TRUE)
    
    M<-5000
    
    #Assume we know mus and p and simulate zs
    samps.mu1<-c()
    samps.mu2<-c()
    samps.p<-c()
    samps.mu1[1]=30
    samps.mu2[1]=70
    samps.p[1]=.5
    samps.zs<-matrix(0,nrow=length(data),ncol=M)
    
    #Outer Loop for MCMC
    for(i in 1:M){
      
      ##Get p's and z's
      dist1 <- likelihood(samps.mu1[i], 1, data)
      dist2 <- likelihood(samps.mu2[i], 2, data)
      prob1 <- dist1/(dist1+dist2)
      samps.zs[,i] <- rbinom(length(data), 1, prob1)
      
      ##Separate 0s and 1s in zs
      data.with.class<-data.frame(value=data, group=samps.zs[,i])
      first.group<-data.with.class%>%filter(group==1)%>%select(value)
      second.group<-data.with.class%>%filter(group==0)%>%select(value)
      
      #Get MLE for fit of a truncated normal
      MLE1 <- optim(par = c(0.1, 0.1, 0.1), # initial values for mu, sigma, and gamma
                    fn = lik2, #function to maximize
                    method = "L-BFGS-B", # this method lets set lower bounds
                    lower = c(0.1, 0.1, -0.99527), # lower limit for parameters
                    upper = c(100, 100, 0.99527),
                    control = list(fnscale = -1), #maximize the function
                    dat = first.group$value
      )
      sig1 <- MLE1$par[2]
      
      MLE2 <- optim(par = c(0.1, 0.1, 0.1), # initial values for mu, sigma, and gamma
                    fn = lik2, #function to maximize
                    method = "L-BFGS-B", # this method lets set lower bounds
                    lower = c(0.1, 0.1, -0.99527), # lower limit for parameters
                    upper = c(100, 100, 0.99527),
                    control = list(fnscale = -1), #maximize the function
                    dat = second.group$value
      )
      
      sig2 <- MLE2$par[2]
      
      ##Simulate a p
      alpha.param<-nrow(first.group)
      beta.param<-nrow(data.with.class)
      samps.p[i]<-rbeta(1,alpha.param,beta.param)
      
      ##Estimate mu1
      mu1.old<-samps.mu1[i]
      mu1.prop<-rnorm(1,mu1.old,1)
      if(mu1.prop>100|mu1.prop<0){
        samps.mu1[length(samps.mu1)+1] <-mu1.old
      }else{
        ratio<-log.lik.tnorm(mu1.prop, sig = sig1, first.group$value)+prior(mu1.prop,muprior=50)-log.lik.tnorm(mu1.old,sig = sig1, first.group$value)-prior(mu1.old,muprior=50)
        alph<-log(runif(1))
        if(alph<ratio){
          samps.mu1[length(samps.mu1)+1]<-mu1.prop
        }else{
          samps.mu1[length(samps.mu1)+1]<-mu1.old
        }
      }
      
      ##Estimate mu2
      mu2.old<-samps.mu2[i]
      mu2.prop<-rnorm(1,mu2.old,1)
      if(mu2.prop>100|mu2.prop<0){
        samps.mu2[length(samps.mu2)+1] <-mu1.old
      }else{
        ratio<-log.lik.tnorm(mu2.prop,sig = sig2, second.group$value)+prior(mu2.prop,muprior=80)-log.lik.tnorm(mu2.old,sig = sig2, second.group$value)-prior(mu2.old,muprior=80)
        alph<-log(runif(1))
        if(alph<ratio){
          samps.mu2[length(samps.mu2)+1]<-mu2.prop
        }else{
          samps.mu2[length(samps.mu2)+1]<-mu2.old
        }
      }
      if(samps.mu1[i]>samps.mu2[i]){
        ph<-samps.mu2[i]
        ph2<-samps.mu1[i]
        samps.mu1[i]<-ph
        samps.mu2[i]<-ph2
      }
    }
    post_mu1 <- samps.mu1[-c(1:500)]
    post_mu2 <- samps.mu2[-c(1:500)]
    
    z_results <- data.frame(matrix(0, ncol = 5, nrow = length(w)))
    finalz <- apply(samps.zs, 1, mean)
    for (k in 1:length(w)) {
      z_results[k, 1] <- data[k]
      z_results[k, 2] <-w[k]
      z_results[k, 3] <- finalz[k]
      if (z_results[k,3] >= 0.5) {z_results[k, 4] <- 1}
      if (z_results[k,4] == w[k]) {z_results[k, 5] <- 1}
    }
    
    results[count, 1] <- l
    results[count, 2] <- mean1.range[j]
    results[count, 3] <- mean2.range[j]
    results[count, 4] <- median(post_mu1)
    results[count, 5] <- sort(post_mu1)[.025*length(post_mu1)] #LB
    results[count, 6] <- sort(post_mu1)[.975*length(post_mu1)] #UB
    results[count, 7] <- results[count,6] - results[count,5] #width
    if (true.mean1 >= results[count, 5] && true.mean1 <= results[count, 6]) {results[count,8] <- 1} #truth
    results[count, 9] <- median(post_mu2)
    results[count, 10] <- sort(post_mu2)[.025*length(post_mu2)] #LB
    results[count, 11] <- sort(post_mu2)[.975*length(post_mu2)] #UB
    results[count, 12] <- results[count,11] - results[count,10] #width
    if (true.mean2 >= results[count, 10] && true.mean2 <= results[count, 11]) {results[count,13] <- 1}
    results[count, 14] <- mean(z_results[,5])
    results[count, 15] <- lb_tdist
    results[count, 16] <- ub_tdist
    results[count, 17] <- ub_tdist - lb_tdist #width
    if (true.mean1 >= lb_tdist && true.mean1 <= ub_tdist) {results[count,18] <- 1}
    if (true.mean2 >= lb_tdist && true.mean2 <= ub_tdist) {results[count,19] <- 1}
    results[count, 20] <- ci_bca_expanded[1]
    results[count, 21] <- ci_bca_expanded[2]
    results[count, 22] <- ci_bca_expanded[2] - ci_bca_expanded[1] #width
    if (true.mean1 >= ci_bca_expanded[1] && true.mean1 <= ci_bca_expanded[2]) {results[count,23] <- 1}
    if (true.mean2 >= ci_bca_expanded[1] && true.mean2 <= ci_bca_expanded[2]) {results[count,24] <- 1}
    
    print(paste0("Done with trial ", l))
    count <- count + 1
  } #end k loop for trials
  
  print(paste0("done with mean group ", j))
  
}# end j loop for means

end_time <- Sys.time()
total_time <- end_time - start_time


write.csv(results, "mixture_mle_mean40and60.csv")

