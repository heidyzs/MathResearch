library(sn) # sampling from the skew-normal distribution is accomplished via the 'sn' package in R
library(e1071) # needed to calculate sample skewness
library(resample)
library(tidyverse)
library(msm) #used for a truncated normal
library(rstan)
#FUNCTIONS#
#None at this time
####GLOBAL DECLARATIONS###

set.seed(8)
N.trials <- 100 # number of Monte Carlo trials to run (e.g., to estimate the CI coverage probabilities); 1000 default
N.bootstrap <- 1000 # size of bootstrap sample for trial i; 1000 default
N.pop <- 100
#N.pop.range<-seq(from=3, to=30, by = 5) # N.pop to try (was3 to 30)
sn.skew.range <-c(0) # sn_skew to try (was -0.9 to 0.9)

sn_mean <- 50 # the mean of the parent SN; based on Table 8.6 in Sauro & Lewis (2016, p. 205)
sn_stdev <- 5 # the standard deviation of the parent SN; based on Table 8.6 in Sauro & Lewis (2016, p. 205)


sd.prior.range <- seq(from = 5, to = 30, by = 5) #standard deviation of the prior distribution, with muprior set to 70

### START EXPERIMENTATION #####

#Build data frame to hold final results for each trial
results <- data.frame(matrix(0, ncol = 15, nrow = (N.trials*length(sn.skew.range))))

colnames(results) <- c("trial", "N.pop", "skew", "lb.tdist", "ub.tdist", "cover.tdist", "width.tdist",
                       "lb.bca.exp", "ub.bca.exp", "cover.bca.exp", "width.bca.exp",
                       "lb.bayes.sdprior", "ub.bayes.sdprior", "cover.bayes.sdprior", "width.bayes.sdprior")

count <- 1 #used to move to the next row in the final results table

for(k in sn.skew.range) {
  
  sn_skew <- k
  cp <- c(mean=sn_mean, s.d.=sn_stdev, gamma1=sn_skew)
  dp <- cp2dp(cp, family="SN") # the cp2dp function converts the centered parameterization to the direct parameterization 
  LB.for.prob<-psn(0, dp=dp) # finds the probability associated with the minimum SUS
  UB.for.prob<-psn(100, dp=dp) # finds the probability associated with the maximum SUS
  
  func<-function(x){
    x*dsn(x,dp[1],dp[2],dp[3])*(1/(UB.for.prob-LB.for.prob))
  }
  
  true.mean<-integrate(func,0,100)[1]$value
  
  for (i in 1:N.trials) {
    
    # generate N.pop random variates from the SN parent within SUS' bounds [0,100]
    y<-qsn(runif(n=N.pop,min=LB.for.prob,max=UB.for.prob), dp=dp, solver="RFB")
    # generate bootstrap sample from the random variates
    b <- bootstrap(y, mean(y), R = N.bootstrap)  
    # generate 95% CIs from the bootstrap sample
    lb.tdist <- mean(y) - qt(0.975, df = length(y)-1)*sd(y)/sqrt(length(y))
    ub.tdist <- mean(y) + qt(0.975, df = length(y)-1)*sd(y)/sqrt(length(y))
    ci_bca_expanded <- CI.bca(b, probs = c(0.025, 0.975), expand = TRUE)
    
    #Format for Stan
    dat<-list(N=length(y),
              y=y,
              sigma=sd(y))
    fit <- stan(file = 'stanmod.stan', data = dat,refresh=0)
    bayes.lb<-summary(fit, pars = c("mu"), probs = c(0.05, 0.95))$summary[[4]]
    bayes.ub<-summary(fit, pars = c("mu"), probs = c(0.05, 0.95))$summary[[5]]
    
    
    #summarize results per skew
    results[count,1] <- i #trial
    results[count,2] <- N.pop
    results[count,3] <- k #skew
    results[count,4] <- lb.tdist
    results[count,5] <- ub.tdist
    if ((true.mean >= lb.tdist) && (true.mean <=ub.tdist)) {results[count,6] <- 1}
    results[count,7] <- ub.tdist-lb.tdist #width of ci
    results[count,8] <- ci_bca_expanded[1]
    results[count,9] <- ci_bca_expanded[2]
    if ((true.mean >= ci_bca_expanded[1]) && (true.mean <= ci_bca_expanded[2])) {results[count,10] <- 1}
    results[count,11] <- ci_bca_expanded[2] - ci_bca_expanded[1] #width of ci
    results[count,12] <- bayes.lb #lowerbound for bayes
    results[count,13] <- bayes.ub #uppperbound for bayes
    if ((true.mean >= bayes.lb) && (true.mean <= bayes.ub)) {results[count,14] <- 1}
    results[count,15] <- bayes.ub-bayes.lb #Width of bayes
    
    count <- count + 1
    print(paste("Done with trial = ", i))
    
  }#end i for trials
  
  print(paste("Done with skew = ", k))
  
  
}#end k for skew


write.csv(results, "C:/Users/Nicholas.Clark/Desktop/Research/MathResearch/Results/Bayesbadmustan.csv")

