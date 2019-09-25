library(sn) # sampling from the skew-normal distribution is accomplished via the 'sn' package in R
library(e1071) # needed to calculate sample skewness

###Global Declarations###

sn_mean <- 68 
sn_stdev <- 12.5 
sn_skew <- -0.5  
cp <- c(mean=sn_mean, s.d.=sn_stdev, gamma1=sn_skew)
dp <- cp2dp(cp, family="SN") # the cp2dp function converts the centered parameterization to the direct parameterization 
LB.for.prob<-psn(0, dp=dp) # finds the probability associated with the minimum SUS
UB.for.prob<-psn(100, dp=dp) # finds the probability associated with the maximum SUS

#play with n here
dat<-qsn(runif(n=10,min=LB.for.prob,max=UB.for.prob), dp=dp, solver="RFB")

###########Functions created to use in simulation#####################

log.lik<-function(x,mu, sdev, sk){
  bayes_mean <- mu 
  cp <- c(mean=bayes_mean, s.d.=sdev, gamma1 = sk)
  dp <- cp2dp(cp, family="SN")
  sum(dsn(x,dp=dp,log=T)-log(UB.for.prob-LB.for.prob))
}

#This function takes your proposed mu and sd, and makes alpha and beta parameters to build a beta distribution for your prior
estBetaParams <- function(mu.scaled, var.scaled) {
  mu<-mu.scaled/100
  var<-var.scaled/100^2
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}


#Play with sd prior
prior<-function(mu,muprior=70,sdprior=20){
  varprior<-sdprior^2
  convert<-estBetaParams(muprior,varprior)
  alpha<-convert$alpha
  beta<-convert$beta
  dbeta(mu/100,alpha,beta,log=TRUE)
}

##View what your prior distribution looks like
look_at_prior<-function(muprior,sdprior){
  varprior<-sdprior^2
  pars<-estBetaParams(muprior,varprior)
  alph<-pars$alpha
  bet<-pars$beta
  prior.vals.2<-rbeta(10000,alph,bet)*100
  return(plot(density(prior.vals.2)))
}

#Simulation#

num.runs<-50000
sample_results <- data.frame(matrix(0, ncol = 4, nrow = num.runs-1))
colnames(sample_results) <- c("sample", "mu", "sig", "gamma")

samps.mu<-c()
samps.sig <-c()
samps.gam <- c()
samps.mu[1]<-50 #mu ranges from 0 to 100
samps.sig[1]<-0.5 #sigma ranges from 0 to 1
samps.gam[1] <- 0 #gamma ranges from -1 to 1
count <- 1
for(j in 2:num.runs){
  sample_results[count, 1] <- j
  mu.old<-samps.mu[j-1]
  mu.prop<-rnorm(1,mu.old,1)
  if(mu.prop>100|mu.prop<0){
    samps.mu[j]<-mu.old
  }else{
    ratio<-log.lik(dat,mu.prop,samps.sig[j-1],samps.gam[j-1])+prior(mu.prop)-log.lik(dat,mu.old,samps.sig[j-1],samps.gam[j-1])-prior(mu.old)
    alph<-log(runif(1))
    if(alph<ratio){
      samps.mu[j]<-mu.prop
    }else{
      samps.mu[j]<-mu.old
    }
  }
  sample_results[count, 2] <- samps.mu[j]
  sig.old<-samps.sig[j-1]
  sig.prop<-rnorm(1,sig.old, .1)
  if(sig.prop>100|sig.prop<0){
    samps.sig[j]<-sig.old
  }else{
    ratio<-log.lik(dat,samps.mu[j],sig.prop,samps.gam[j-1])-log.lik(dat,samps.mu[j], sig.old, samps.gam[j-1])
    alph<-log(runif(1))
    if(alph<ratio){
      samps.sig[j]<-sig.prop
    }else{
      samps.sig[j]<-sig.old
    }
  }
  sample_results[count, 3] <- samps.sig[j]
  gam.old<-samps.gam[j-1]
  gam.prop<-rnorm(1, gam.old, .05)
  if(gam.prop>0.99527|gam.prop< -0.99527){
    samps.gam[j]<-gam.old
  }else{
    ratio<-log.lik(dat,samps.mu[j],samps.sig[j],gam.prop)-log.lik(dat,samps.mu[j], samps.sig[j], gam.old)
    alph<-log(runif(1))
    if(alph<ratio){
      samps.gam[j]<-gam.prop
    }else{
      samps.gam[j]<-gam.old
    }
  }
  sample_results[count, 4] <- samps.gam[j]
  count <- count + 1
}


hist(samps.gam)
hist(samps.mu)
hist(samps.sig)

posterior_samples <- sample_results[-c(1:5000),]

plot(samps.mu,samps.sig)
median(posterior_samples$mu)
median(posterior_samples$sig)
median(posterior_samples$gamma)


sort(posterior_samples$mu)[.975*length(posterior_samples$mu)]-sort(posterior_samples$mu)[.025*length(posterior_samples$mu)]

#TRUTH


gah<-cp2dp(c(68,12.5,.5),"SN")
ub<-psn(100,gah[1],gah[2],gah[3])
lb<-psn(0,gah[1],gah[2],gah[3])

first.moment<-function(x){
  x*dsn(x,gah[1],gah[2],gah[3])*(1/(ub-lb))
}


second.moment<-function(x){
  x^2*dsn(x,gah[1],gah[2],gah[3])*(1/(ub-lb))
}
true.mean<-integrate(first.moment,0,100)[1]$value

true.var<-integrate(second.moment,0,100)[1]$value-(true.mean)^2


######Experimentation#####

####GLOBAL DECLARATIONS###

set.seed(8)
N.trials <- 100 # number of Monte Carlo trials to run (e.g., to estimate the CI coverage probabilities); 1000 default
N.bootstrap <- 1000 # size of bootstrap sample for trial i; 1000 default
N.pop <- 10
#N.pop.range<-seq(from=3, to=30, by = 5) # N.pop to try (was3 to 30)
sn.skew.range <-c(-0.9, -0.5, 0.5, 0.9) # sn_skew to try (was -0.9 to 0.9)

sn_mean <- 68 # the mean of the parent SN; based on Table 8.6 in Sauro & Lewis (2016, p. 205)
sn_stdev <- 12.5 # the standard deviation of the parent SN; based on Table 8.6 in Sauro & Lewis (2016, p. 205)


sd.prior = 20 #standard deviation of the prior distribution, with muprior set to 70

### START EXPERIMENTATION #####

#Build data frame to hold final results for each trial
results <- data.frame(matrix(0, ncol = 18, nrow = (N.trials*length(sn.skew.range))))

colnames(results) <- c("trial", "N.pop", "skew", "lb.tdist", "ub.tdist", "cover.tdist", "width.tdist",
                       "lb.bca.exp", "ub.bca.exp", "cover.bca.exp", "width.bca.exp", "median.posterior.mu",
                       "lb.bayes.sdprior20", "ub.bayes.sdprior20", "width.bayes.sdprior20", "cover.bayes.sdprior20",
                       "median.posterior.sig","median.posterior.gam")
                       

count <- 1 #used to move to the next row in the final results table
start_time = Sys.time()

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

    #Bayes Loop
    num.runs<-10000
    samps.mu<-c()
    samps.sig <-c()
    samps.gam <- c()
    samps.mu[1]<-50 #mu ranges from 0 to 100
    samps.sig[1]<-0.5 #sigma ranges from 0 to 1
    samps.gam[1] <- 0 #gamma ranges from -0.99527 to 0.99527
    for(j in 2:num.runs){
      mu.old<-samps.mu[j-1]
      mu.prop<-rnorm(1,mu.old,1)
      if(mu.prop>100|mu.prop<0){
        samps.mu[j]<-mu.old
      }else{
        ratio<-log.lik(y,mu.prop,samps.sig[j-1],samps.gam[j-1])+prior(mu.prop)-log.lik(y,mu.old,samps.sig[j-1],samps.gam[j-1])-prior(mu.old)
        alph<-log(runif(1))
        if(alph<ratio){
          samps.mu[j]<-mu.prop
        }else{
          samps.mu[j]<-mu.old
        }
      }
      sig.old<-samps.sig[j-1]
      sig.prop<-rnorm(1,sig.old, .1)
      if(sig.prop>100|sig.prop<0){
        samps.sig[j]<-sig.old
      }else{
        ratio<-log.lik(y,samps.mu[j],sig.prop,samps.gam[j-1])-log.lik(y,samps.mu[j], sig.old, samps.gam[j-1])
        alph<-log(runif(1))
        if(alph<ratio){
          samps.sig[j]<-sig.prop
        }else{
          samps.sig[j]<-sig.old
        }
      }
      gam.old<-samps.gam[j-1]
      gam.prop<-rnorm(1, gam.old, .05)
      if(gam.prop>0.99527|gam.prop< -0.99527){
        samps.gam[j]<-gam.old
      }else{
        ratio<-log.lik(y,samps.mu[j],samps.sig[j],gam.prop)-log.lik(y,samps.mu[j], samps.sig[j], gam.old)
        alph<-log(runif(1))
        if(alph<ratio){
          samps.gam[j]<-gam.prop
        }else{
          samps.gam[j]<-gam.old
        }
      }
    }
    posterior.mu <- samps.mu[-c(1:100)]
    posterior.sig <- samps.sig[-c(1:100)]
    posterior.gam <- samps.gam[-c(1:100)]
    
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
    results[count,12] <- median(posterior.mu)
    results[count,13] <- sort(posterior.mu)[.025*length(posterior.mu)]
    results[count,14] <- sort(posterior.mu)[.975*length(posterior.mu)]
    results[count,15] <- results[count,14] - results[count,13] #width of ci
    if (true.mean >= results[count,13] && true.mean <= results[count,14]) {results[count,16] <-1} #cover
    results[count,17] <- median(posterior.sig) 
    results[count,18] <- median(posterior.gam) 

    count <- count + 1
    print(paste("Done with trial = ", i))
    
  }#end i for trials
  
  print(paste("Done with skew = ", k))
  
  
}#end k for skew

end_time = Sys.time()
end_time - start_time

#write.csv(results, "full_bayes_with_sdprior_20.csv")

#Analyze Results
library(tidyverse)
library(ggplot2)

full_bayes <- read.csv("full_bayes_with_sdprior_20_v2.csv")

cover_percentages <- full_bayes %>% select("skew", "cover.bayes.sdprior20", "cover.tdist", "cover.bca.exp") %>% 
  group_by(skew) %>% summarize(
    cover.perc.bayes = sum(cover.bayes.sdprior20),
    cover.perc.tdist = sum(cover.tdist),
    cover.perc.bca.exp = sum(cover.bca.exp))


#TRUTH
gah<-cp2dp(c(68,12.5, 0.9),"SN")
ub<-psn(100,gah[1],gah[2],gah[3])
lb<-psn(0,gah[1],gah[2],gah[3])

first.moment<-function(x){
  x*dsn(x,gah[1],gah[2],gah[3])*(1/(ub-lb))
}

second.moment<-function(x){
  x^2*dsn(x,gah[1],gah[2],gah[3])*(1/(ub-lb))
}
true.mean<-integrate(first.moment,0,100)[1]$value
true.var<-integrate(second.moment,0,100)[1]$value-(true.mean)^2


#full_bayes %>% filter(skew == 0.9) %>% summarize(mean(median.posterior.mu))
#full_bayes %>% filter(skew == 0.9) %>% summarize(mean(median.posterior.sig))
#full_bayes %>% filter(skew == 0.9) %>% summarize(mean(median.posterior.gam))


#Making Data Tidy To Plot######
width_results <- full_bayes %>% select("trial", "N.pop", "skew", starts_with("width"))

names(width_results) <- substring(names(width_results), 7)
new_col_names <- c("Trial", "N.pop", "Skew", "T Dist", "BCaExp", "Bayes SDP20")
colnames(width_results) <- new_col_names

width_results <- width_results %>% gather("T Dist", "BCaExp", "Bayes SDP20", 
                                          key = "CI.type", value = "width")


#Generate boxplot for different width CIs

means <- aggregate(width ~ CI.type + Skew, width_results, mean) #calculate means by CI.type to generate labels

width_results %>% group_by(CI.type) %>% arrange(desc(CI.type)) %>% 
  ggplot(aes(x = as.factor(CI.type), y = width)) +
  geom_boxplot() +
  geom_hline(yintercept = 6.2) +
  geom_hline(yintercept = 20) +
  stat_summary(fun.y = mean, colour = "darkred", geom="point", shape = 18, size = 3, show.legend = FALSE)+
  geom_text_repel(data = means, aes(label = round(width, digits = 2))) +
  facet_wrap(~Skew, labeller = label_both) +
  coord_flip()

