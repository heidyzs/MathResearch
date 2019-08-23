library(sn) # sampling from the skew-normal distribution is accomplished via the 'sn' package in R
library(e1071) # needed to calculate sample skewness
library(resample) # a package that allows for the expanded version of the percentile bootstrap

#Functions##

#Used to integrate and find true.mean
func<-function(x){
  x*dsn(x,dp[1],dp[2],dp[3])*(1/(UB.for.prob-LB.for.prob))
}

###########Functions created to use in simulation#####################

log.lik<-function(x,mu){
  bayes_mean <- mu 
  cp <- c(mean=bayes_mean, s.d.=sn_stdev, gamma1=sn_skew)
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
prior<-function(mu,muprior=70,sdprior=25){
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


look_at_prior2<-function(muprior,sdprior){
  varprior<-sdprior^2
  pars<-estBetaParams(muprior,varprior)
  alph<-pars$alpha
  bet<-pars$beta
  prior.vals.2<-rbeta(10000,alph,bet)*100
  return(prior.vals.2)
}



####GLOBAL DECLARATIONS###

set.seed(8)
N.trials <- 500 # number of Monte Carlo trials to run (e.g., to estimate the CI coverage probabilities); 1000 default
N.bootstrap <- 1000 # size of bootstrap sample for trial i; 1000 default
N.pop <- 10
#N.pop.range<-seq(from=3, to=30, by = 5) # N.pop to try (was3 to 30)
sn.skew.range <-c(-0.9, -0.5, 0.5, 0.9) # sn_skew to try (was -0.9 to 0.9)

sn_mean <- 68 # the mean of the parent SN; based on Table 8.6 in Sauro & Lewis (2016, p. 205)
sn_stdev <- 12.5 # the standard deviation of the parent SN; based on Table 8.6 in Sauro & Lewis (2016, p. 205)


sd.prior.range <- seq(from = 5, to = 30, by = 5) #standard deviation of the prior distribution, with muprior set to 70

### START EXPERIMENTATION #####

#Build data frame to hold final results for each trial
results <- data.frame(matrix(0, ncol = 35, nrow = (N.trials*length(sn.skew.range))))

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

    #buildiing dataframe to hold only the bayes results by sdprior
    bayes_results <- data.frame(matrix(0, ncol = 5, nrow = length(sdprior.range)))
    bayes_count <- 1 #used to move to the next row in the bayes results table
    
    #Bayes Loop
    for (l in sdprior.range) {
      samps<-c()
      samps[1]<-50
        for(j in 2:5000){
          mu.old<-samps[j-1]
          mu.prop<-rnorm(1,mu.old,1)
          if(mu.prop>100|mu.prop<0){
            samps[j]<-mu.old
          }else{
            ratio<-log.lik(y,mu.prop)+prior(mu.prop, muprior = 70, sdprior = l)-log.lik(y,mu.old)-prior(mu.old, muprior = 70, sdprior = l)
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
    
    count <- count + 1
    print(paste("Done with trial = ", i))
    
    }#end i for trials
  
  print(paste("Done with skew = ", k))


}#end k for skew


#write.csv(results, "bayes_vs_bootstrap_full_data_set.csv")



#####Analyze Results#########
library(tidyverse)
library(gridExtra)


#width results needs to be tidy by CI.Type

combined_results <- read.csv("bayes_vs_bootstrap_full_data_set.csv")
width_results <- combined_results %>% select("trial", "N.pop", "skew", starts_with("width"))

names(width_results) <- substring(names(width_results), 7)
new_col_names <- c("Trial", "N.pop", "Skew", "T Dist", "BCaExp", "Bayes SDP5", 
                   "Bayes SDP10", "Bayes SDP15", "Bayes SDP20", "Bayes SDP25",
                   "Bayes SDP30")
colnames(width_results) <- new_col_names

width_results <- width_results %>% gather("T Dist", "BCaExp", "Bayes SDP5", 
                                          "Bayes SDP10", "Bayes SDP15", "Bayes SDP20", "Bayes SDP25",
                                          "Bayes SDP30",
                                          key = "CI.type", value = "width")


#Generate boxplot for different width CIs

means <- aggregate(width ~ CI.type + Skew, width_results, mean) #calculate means by CI.type to generate labels

width_plot <- width_results %>% group_by(CI.type) %>% arrange(desc(CI.type)) %>% 
  ggplot(aes(x = as.factor(CI.type), y = width)) +
  geom_boxplot() +
  geom_hline(yintercept = 6.2) +
  geom_hline(yintercept = 20) +
  stat_summary(fun.y = mean, colour = "darkred", geom="point", shape = 18, size = 3, show.legend = FALSE)+
  geom_text_repel(data = means, aes(label = round(width, digits = 2))) +
  facet_wrap(~Skew, labeller = label_both) +
  coord_flip()

##Filtered out T distribution and BCaExp to focus just on Bayes
width_results_bayes <- width_results %>% 
  group_by(CI.type) %>% filter(CI.type != "T Dist") %>% 
  filter(CI.type != "BCaExp")

means_bayes <- aggregate(width ~ CI.type + Skew, width_results_bayes, mean) #calculate means by CI.type to generate labels

width_results_bayes %>% group_by(CI.type) %>% arrange(desc(CI.type)) %>% 
  ggplot(aes(x = as.factor(CI.type), y = width)) +
  geom_boxplot() +
  geom_hline(yintercept = 6.2) +
  geom_hline(yintercept = 20) +
  stat_summary(fun.y = mean, colour = "darkred", geom="point", shape = 18, size = 3, show.legend = FALSE)+
  geom_text_repel(data = means_bayes, aes(label = round(width, digits = 2))) +
  facet_wrap(~Skew, labeller = label_both) +
  coord_flip()

#Plot of Different prior distributions
layout.matrix <- matrix(c(0,1,2,7,3,4,0,5,6), nrow = 3, ncol = 3)
layout(mat = layout.matrix, heights = c(2, 1, 1), widths = c(1, 1, 1))
layout.show(7)


par(mfrow = c(2, 3))
par(mar = c(5.1, 4.1, 4.1, 2.1))
plot(density(look_at_prior2(70, 5)), main = "Prior with Muprior = 70 and SDprior = 5")
plot(density(look_at_prior2(70, 10)), main = "Prior with Muprior = 70 and SDprior = 10")
plot(density(look_at_prior2(70, 15)), main = "Prior with Muprior = 70 and SDprior = 15")
plot(density(look_at_prior2(70, 20)), main = "Prior with Muprior = 70 and SDprior = 20")
plot(density(look_at_prior2(70, 25)), main = "Prior with Muprior = 70 and SDprior = 25")
plot(density(look_at_prior2(70, 30)), main = "Prior with Muprior = 70 and SDprior = 30")



grid.arrange(width_plot, combined_priors, nrow = 2)
width_plot
