library(sn) # sampling from the skew-normal distribution is accomplished via the 'sn' package in R
library(e1071) # needed to calculate sample skewness


##################################################################################


common.sd<-12.5
sn_stdev<-common.sd


#FUNCTIONS#
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
log.lik<-function(mu,component,dat){
  if ((length(dat) == 0)) {
    return(0)
  } else {
    sn_mean<-mu
    cp<-c(mean=sn_mean, s.d.=common.sd,gamma1= .5*(-1)^(component-1))
    dp<-cp2dp(cp,family="SN")
    LB<-psn(0,dp=dp)
    UB<-psn(100,dp=dp)
    sum(dsn(dat,dp=dp,log=T)-log(UB-LB))
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

look_at_prior<-function(muprior,sdprior){
  varprior<-sdprior^2
  pars<-estBetaParams(muprior,varprior)
  alph<-pars$alpha
  bet<-pars$beta
  prior.vals.2<-rbeta(10000,alph,bet)*100
  return(plot(density(prior.vals.2)))
}

#Start Simulation to Gather Data
set.seed(8)

n <- 25


#mean1.range <- seq(from = 45, to = 10, by = -5)
#mean2.range <- seq(from = 55, to = 90, by = 5 )
N.trials <-100

results <- data.frame(matrix(0, ncol = 14, nrow = N.trials))
x <- c("trial", "mean1", "mean2", "median.sampsmu1", "lb.mu1", "ub.mu1", "width.mu1", "cover.mu1",
       "median.samps.mu2", "lb.mu2", "ub.mu2", "width.mu2", "cover.mu2", "perc.correct")
colnames(results) <- x


start_time <- Sys.time()
count <- 1


#for (j in 1:length(mean1.range)) {
  
  #Distribution 1
  sn_mean1 <- 60
  sn_stdev1 <- common.sd 
  sn_skew1 <- 0.5  
  cp1 <- c(mean=sn_mean1, s.d.=sn_stdev1, gamma1=sn_skew1)
  dp1 <- cp2dp(cp1, family="SN") # the cp2dp function converts the centered parameterization to the direct parameterization 
  LB.for.prob1<-psn(0, dp=dp1) # finds the probability associated with the minimum SUS
  UB.for.prob1<-psn(100, dp=dp1) # finds the probability associated with the maximum SUS
  
  #Distribution 2
  sn_mean2 <- 80
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
    
  for (l in 1:N.trials) {
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
      dist1 <- likelihood(samps.mu1[i],1,data)
      dist2 <- likelihood(samps.mu2[i],2,data)
      prob1 <- dist1/(dist1+dist2)
      samps.zs[,i] <- rbinom(length(data), 1, prob1)
      
      ##Separate 0s and 1s in zs
      data.with.class<-data.frame(value=data, group=samps.zs[,i])
      first.group<-data.with.class%>%filter(group==1)%>%select(value)
      second.group<-data.with.class%>%filter(group==0)%>%select(value)
      
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
        ratio<-log.lik(mu1.prop,1,first.group$value)+prior(mu1.prop,muprior=50)-log.lik(mu1.old,1,first.group$value)-prior(mu1.old,muprior=50)
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
        ratio<-log.lik(mu2.prop,2,second.group$value)+prior(mu2.prop,muprior=80)-log.lik(mu2.old,2,second.group$value)-prior(mu2.old,muprior=80)
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
    results[count, 2] <- sn_mean1
    results[count, 3] <- sn_mean2
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
    
    print(paste0("Done with trial ", l))
    count <- count + 1
  } #end k loop for trials
  
  #print(paste0("done with mean group ", j))
  
#}# end j loop for means

end_time = Sys.time()
end_time - start_time 


write.csv(results, "./Mixture/mixture_mean60_mean_80_n.pop25_trial100.csv")



#Analyzing the results 
names(results)
results

results %>% 
  group_by(mean1, mean2) %>% 
  summarise(
    avgcover1 = mean(cover.mu1),
    avgcover2 = mean(cover.mu2),
    #avgwidth1 = mean(width.mu1), 
    #avgwidth2 = mean(width.mu2),
    #avglb1 = mean(lb.mu1), 
    #avglb2 = mean(lb.mu2),
    avgmu1 = mean(median.sampsmu1),
    avgmu2 = mean(median.samps.mu2),
    avg.correct = mean(perc.correct))


#Look at the distributions
d1 <- density(post_mu1) 
d2 <- density(post_mu2)

plot(range(d1$x, d2$x), range(d1$y, d2$y), type = "n", 
     xlab = "x", ylab = "Density")
abline(v = true.mean1)
lines(d1, col = "red") 
lines(d2, col = "blue")


######Playing with mixture results

mixture <- read.csv("./Mixture/mixture_mean70_trial100_muprior80and50.csv")
names(mixture)

mixture %>% 
  group_by(mean1, mean2) %>% 
  summarise(
    avgcover1 = mean(cover.mu1),
    avgcover2 = mean(cover.mu2),
    #avgwidth1 = mean(width.mu1), 
    #avgwidth2 = mean(width.mu2),
    #avglb1 = mean(lb.mu1), 
    #avglb2 = mean(lb.mu2),
    avgmu1 = mean(median.sampsmu1),
    avgmu2 = mean(median.samps.mu2),
    avg.correct = mean(perc.correct))





