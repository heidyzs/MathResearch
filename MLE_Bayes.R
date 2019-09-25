library(sn) # sampling from the skew-normal distribution is accomplished via the 'sn' package in R
library(e1071) # needed to calculate sample skewness
library(resample)
library(tidyverse)
library(msm) #used for a truncated normal

#FUNCTIONS#

#Log Likelihood
log.lik<-function(mu, sig, gam, dat){
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

#Log Likelihood for truncated normal distribution
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


#Used to integrate and find true.mean
func<-function(x){
  x*dsn(x,dp[1],dp[2],dp[3])*(1/(UB.for.prob-LB.for.prob))
}


look_at_prior2<-function(muprior,sdprior){
  varprior<-sdprior^2
  pars<-estBetaParams(muprior,varprior)
  alph<-pars$alpha
  bet<-pars$beta
  prior.vals.2<-rbeta(10000,alph,bet)*100
  return(prior.vals.2)
}

#Likelihood function used for MLE

lik2<-function(parvec, dat){
  cp<-c(mean=parvec[1], s.d.=parvec[2], gamma1= parvec[3])
  dp<-cp2dp(cp,family="SN")
  LB<-psn(0,dp=dp)
  UB<-psn(100,dp=dp)
  sum(dsn(dat,dp[1], dp[2], dp[3], log = T)-log(UB-LB))
}

####GLOBAL DECLARATIONS###

set.seed(8)
N.trials <- 100 # number of Monte Carlo trials to run (e.g., to estimate the CI coverage probabilities); 1000 default
N.bootstrap <- 1000 # size of bootstrap sample for trial i; 1000 default
N.pop <- 10
#N.pop.range<-seq(from=3, to=30, by = 5) # N.pop to try (was3 to 30)
sn.skew.range <-c(-0.9, -0.5, 0, 0.5, 0.9) # sn_skew to try (was -0.9 to 0.9)

sn_mean <- 68 # the mean of the parent SN; based on Table 8.6 in Sauro & Lewis (2016, p. 205)
sn_stdev <- 12.5 # the standard deviation of the parent SN; based on Table 8.6 in Sauro & Lewis (2016, p. 205)


sd.prior.range <- seq(from = 5, to = 30, by = 5) #standard deviation of the prior distribution, with muprior set to 70

### START EXPERIMENTATION #####

#Build data frame to hold final results for each trial
results <- data.frame(matrix(0, ncol = 36, nrow = (N.trials*length(sn.skew.range))))

colnames(results) <- c("trial", "N.pop", "skew", "lb.tdist", "ub.tdist", "cover.tdist", "width.tdist",
                       "lb.bca.exp", "ub.bca.exp", "cover.bca.exp", "width.bca.exp",
                       "lb.bayes.sdprior5", "ub.bayes.sdprior5", "cover.bayes.sdprior5", "width.bayes.sdprior5",
                       "lb.bayes.sdprior10", "ub.bayes.sdprior10", "cover.bayes.sdprior10", "width.bayes.sdprior10",
                       "lb.bayes.sdprior15", "ub.bayes.sdprior15", "cover.bayes.sdprior15", "width.bayes.sdprior15", 
                       "lb.bayes.sdprior20", "ub.bayes.sdprior20", "cover.bayes.sdprior20","width.bayes.sdprior20",
                       "lb.bayes.sdprior25", "ub.bayes.sdprior25", "cover.bayes.sdprior25", "width.bayes.sdprior25",
                       "lb.bayes.sdprior30", "ub.bayes.sdprior30", "cover.bayes.sdprior30", "width.bayes.sdprior30")

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
    
    #Get MLE for bayes
    MLE <- optim(par = c(0.1, 0.1, 0.1), # initial values for mu, sigma, and gamma
                 fn = lik2, #function to maximize
                 method = "L-BFGS-B", # this method lets set lower bounds
                 lower = c(0.1, 0.1, -0.99527), # lower limit for parameters
                 upper = c(100, 100, 0.99527),
                 control = list(fnscale = -1), #maximize the function
                 dat = y
    )
    
    sig <- MLE$par[2]
    #if (abs(MLE$par[3]) == 0.99527) {
    #  gam <- 0
    #} else {
    #gam <- MLE$par[3]
    #}
    
    #buildiing dataframe to hold only the bayes results by sdprior
    bayes_results <- data.frame(matrix(0, ncol = 5, nrow = length(sd.prior.range)))
    bayes_count <- 1 #used to move to the next row in the bayes results table
    
    #Bayes Loop
    for (l in sd.prior.range) {
      samps<-c()
      samps[1]<-50
      for(j in 2:5000){
        mu.old<-samps[j-1]
        mu.prop<-rnorm(1,mu.old,1)
        if(mu.prop>100|mu.prop<0){
          samps[j]<-mu.old
        }else{
          ratio<-log.lik.tnorm(mu.prop, sig, y)+prior(mu.prop, muprior = 70, sdprior = l)-log.lik.tnorm(mu.old, sig, y)-prior(mu.old, muprior = 70, sdprior = l)
          alph<-log(runif(1))
          if(alph<ratio){
            samps[j]<-mu.prop
          }else{
            samps[j]<-mu.old
          }
        }
      }
      posterior.samp <- samps[-(1:500)] #eliminate first 500 "bad" samples
      bayes_results[bayes_count,1] <- l #sdprior
      bayes_results[bayes_count,2] <- sort(posterior.samp)[.025*length(posterior.samp)]
      bayes_results[bayes_count,3] <-sort(posterior.samp)[.975*length(posterior.samp)]
      if (true.mean >= bayes_results[bayes_count,2] && true.mean <= bayes_results[bayes_count,3]) {bayes_results[bayes_count,4] <-1}
      bayes_results[bayes_count,5] <- bayes_results[bayes_count,3]-bayes_results[bayes_count,2] #width of bayes
      
      bayes_count <- bayes_count + 1
      
      
    } #ends "l" loop for sdprior.range
    
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
    results[count,12] <- bayes_results[1,2] #lowerbound for bayes
    results[count,13] <- bayes_results[1,3] #uppperbound for bayes
    results[count,14] <- bayes_results[1,4] #if covered
    results[count,15] <- bayes_results[1,5] #Width of bayes for sdprior = 5
    results[count,16] <- bayes_results[2,2] 
    results[count,17] <- bayes_results[2,3] 
    results[count,18] <- bayes_results[2,4] 
    results[count,19] <- bayes_results[2,5] #width of bayes for sdprior = 10
    results[count,20] <- bayes_results[3,2] 
    results[count,21] <- bayes_results[3,3] 
    results[count,22] <- bayes_results[3,4] 
    results[count,23] <- bayes_results[3,5] #width of bayes for sdprior = 15
    results[count,24] <- bayes_results[4,2] 
    results[count,25] <- bayes_results[4,3] 
    results[count,26] <- bayes_results[4,4] 
    results[count,27] <- bayes_results[4,5] #Width of bayes for sdprior = 20
    results[count,28] <- bayes_results[5,2] 
    results[count,29] <- bayes_results[5,3] 
    results[count,30] <- bayes_results[5,4] 
    results[count,31] <- bayes_results[5,5] #Width of bayes for sdprior = 25
    results[count,32] <- bayes_results[6,2] 
    results[count,33] <- bayes_results[6,3] 
    results[count,34] <- bayes_results[6,4] 
    results[count,35] <- bayes_results[6,5] #Width of bayes for sdprior = 30
    results[count,36] <- 70 #muprior
    
    count <- count + 1
    print(paste("Done with trial = ", i))
    
  }#end i for trials
  
  print(paste("Done with skew = ", k))
  
  
}#end k for skew


#write.csv(results, "./MLE/bayes_fixedgam_muprior70_trial100.csv")



#Analyze Results
mle_results <- read.csv("./MLE/bayes_fixedgam_muprior70_trial100.csv")
width_results <- mle_results %>% select("trial", "N.pop", "skew", starts_with("width"))

names(width_results) <- substring(names(width_results), 7)
new_col_names <- c("Trial", "N.pop", "Skew", "T Dist", "BCaExp", "Bayes SDP5", 
                   "Bayes SDP10", "Bayes SDP15", "Bayes SDP20", "Bayes SDP25",
                   "Bayes SDP30")
colnames(width_results) <- new_col_names

width_results <- width_results %>% gather("T Dist", "BCaExp", "Bayes SDP5", 
                                          "Bayes SDP10", "Bayes SDP15", "Bayes SDP20", 
                                          "Bayes SDP25", "Bayes SDP30",
                                          key = "CI.type", value = "width")

#Generate boxplot for different width CIs

means <- aggregate(width ~ CI.type + Skew, width_results, mean) #calculate means by CI.type to generate labels

width_results %>% group_by(CI.type) %>% arrange(desc(CI.type)) %>% 
  ggplot(aes(x = as.factor(CI.type), y = width)) +
  geom_boxplot() +
  geom_hline(yintercept = 6.2) +
  geom_hline(yintercept = 20) +
  stat_summary(fun.y = mean, colour = "darkred", 
               geom="point", shape = 18, size = 3, show.legend = FALSE)+
  geom_text_repel(data = means, aes(label = round(width, digits = 2))) +
  facet_wrap(~Skew, labeller = label_both) +
  coord_flip()

#Cover Percentages

cover_percentages <- mle_results %>% group_by(skew) %>% 
  select("skew", starts_with("cover")) %>% 
  summarise_each(mean)

colnames(cover_percentages) <- c("skew", "Tdist", "BCaExp", "BayesSDP5", 
                                 "BayesSDP10", "BayesSDP15",
                                 "BayesSD20", "BayesSDP25", "BayesSDP30")
