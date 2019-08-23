library(sn) # sampling from the skew-normal distribution is accomplished via the 'sn' package in R
library(e1071) # needed to calculate sample skewness
library(Cairo) # Used for exporting high quality graphics


#####################################################################

### GLOBAL DECLARATIONS

set.seed(8)
N.trials <- 1000 # number of Monte Carlo trials to run (e.g., to estimate the CI coverage probabilities); 1000 default
N.bootstrap <- 1000 # size of bootstrap sample for trial i; 1000 default
N.pop.range<-seq(from=3, to=30, by = 1) # N.pop to try
sn.skew.range <-1 # can't adjust this parameter here

# used for experimentation on obtaining the desired coverage probability by reducing alpha
# e.g., if equal to 0, then you get a 95% CI; if equal to 0.005, then you get a 96% CI
ci.hw.increase <- 0

sn_mean <- 68 # the mean of the parent SN; based on Table 8.6 in Sauro & Lewis (2016, p. 205)
sn_stdev <- 12.5 # the standard deviation of the parent SN; based on Table 8.6 in Sauro & Lewis (2016, p. 205)
x.bar<-sn_mean/100
s.sq<-12.5^2*(1/100)^2
#Moment Based Estimates for Beta Distribution##
alph.mom<-x.bar*(x.bar*(1-x.bar)/s.sq-1)
bet.mom<-alph.mom*(1-x.bar)/x.bar


results.summary <- data.frame(matrix(0, ncol = 15, nrow = (length(N.pop.range)*length(sn.skew.range)))) # build a preallocated data frame to hold the results
colnames(results.summary) <- c("N.pop","cover.perc","width.perc", "cover.perc.exp","width.perc.exp","cover.perc.clt","width.clt", "cover.perc.tdist", "width.tdist",
                               "cover.bca", "width.bca", "cover.bca.exp", "width.bca.exp", "cover.bootstrapT", "width.bootstrapT")
count <- 1

for (j in N.pop.range) {
  N.pop <- j # sample from parent population (e.g., skew-normal to replicate the respondent's SUS scores)
  
  results <- data.frame(matrix(0, ncol = 29, nrow = N.trials)) # build a preallocated data frame to hold the results
  x <- c("trial", "LB.perc", "RB.perc", "cover.perc","width.perc",
         "LB.perc.exp", "RB.perc.exp", "cover.perc.exp","width.perc.exp","LB.clt","RB.clt","cover.perc.clt","width.clt", "LB.tdist","RB.tdist","cover.perc.tdist","width.tdist",
         "LB.bca", "RB.bca", "cover.bca", "width.bca", "LB.bca.exp", "RB.bca,exp", "width.bca.exp", "cover.bca.exp",
         "LB.bootstrapT", "RB.bootstrapT", "cover.bootstrapT", "width.bootstrapT")
  
  colnames(results) <- x
  
  # start the simulation
  for (i in 1:N.trials) {
    # generate N.pop random variates from the SN parent within SUS' bounds [0,100]
    #y<-qsn(runif(n=N.pop,min=LB.for.prob,max=UB.for.prob), dp=dp, solver="RFB")
    y<-rbeta(n=N.pop,alph.mom,bet.mom)*100
    # generate bootstrap sample from the random variates
    #b <- resample::bootstrap(y, mean(y), R = N.bootstrap)
    #bootsrap sample using mean and standard deviation for bootstrapT
    b2 <- bootstrap(y, c(mean = mean(y), sd = sd(y)), R = N.bootstrap)
    # generate 95% CIs from the bootstrap sample
    ci_percentile <- CI.percentile(b2, probs = c(0.025-ci.hw.increase, 0.975+ci.hw.increase), expand = FALSE)
    ci_percentile_expanded <- CI.percentile(b2, probs = c(0.025-ci.hw.increase, 0.975+ci.hw.increase), expand = TRUE)
    lb.clt<-mean(y)-1.96*sd(y)/sqrt(length(y))
    ub.clt<-mean(y)+1.96*sd(y)/sqrt(length(y))
    lb.tdist <- mean(y) - qt(0.975, df = length(y)-1)*sd(y)/sqrt(length(y))
    ub.tdist <- mean(y) + qt(0.975, df = length(y)-1)*sd(y)/sqrt(length(y))
    ci_bca <- CI.bca(b2, probs = c(0.025-ci.hw.increase, 0.975+ci.hw.increase), expand = FALSE)
    ci_bca_expanded <- CI.bca(b2, probs = c(0.025-ci.hw.increase, 0.975+ci.hw.increase), expand = TRUE)
    ci_bootstrapT <- CI.bootstrapT(b2, probs = c(0.025-ci.hw.increase, 0.975+ci.hw.increase))
    results[i,1] <- i
    results[i,2] <- ci_percentile[1,1]
    results[i,3] <- ci_percentile[1,2]
    if ((sn_mean >= ci_percentile[1,1]) && (sn_mean <= ci_percentile[1,2])) {results[i,4] <- 1}
    results[i,5] <- ci_percentile[1,2] - ci_percentile[1,1]
    results[i,6] <- ci_percentile_expanded[1,1]
    results[i,7] <- ci_percentile_expanded[1,2]
    if ((sn_mean >= ci_percentile_expanded[1,1]) && (sn_mean <= ci_percentile_expanded[1,2])) {results[i,8] <- 1}
    results[i,9] <- ci_percentile_expanded[1,2] - ci_percentile_expanded[1,1]
    results[i,10]<- lb.clt
    results[i,11]<-ub.clt
    if ((sn_mean >= lb.clt) && (sn_mean <=ub.clt)) {results[i,12] <- 1}
    results[i,13] <- ub.clt-lb.clt
    results[i,14]<- lb.tdist
    results[i,15]<-ub.tdist
    if ((sn_mean >= lb.tdist) && (sn_mean <=ub.tdist)) {results[i,16] <- 1}
    results[i,17] <- ub.tdist-lb.tdist
    results[i,18] <- ci_bca[1,1]
    results[i,19] <- ci_bca[1,2]
    if ((sn_mean >= ci_bca[1,1]) && (sn_mean <= ci_bca[1,2])) {results[i,20] <- 1}
    results[i,21] <- ci_bca[1,2] - ci_bca[1,1]
    results[i,22] <- ci_bca_expanded[1,1]
    results[i,23] <- ci_bca_expanded[1,2]
    if ((sn_mean >= ci_bca_expanded[1,1]) && (sn_mean <= ci_bca_expanded[1,2])) {results[i,24] <- 1}
    results[i,25] <- ci_bca_expanded[1,2] - ci_bca_expanded[1,1]
    results[i,26] <- ci_bootstrapT[1]
    results[i,27] <- ci_bootstrapT[2]
    if ((sn_mean >= ci_bootstrapT[1]) && (sn_mean <= ci_bootstrapT[2])) {results[i,28] <- 1}
    results[i,29] <- ci_bootstrapT[2] - ci_bootstrapT[1]
  } # end for (i in 1:N.trials)
  
  # summarize results
  results.summary[count,1] <- j
  results.summary[count,2] <- mean(results[,4])
  results.summary[count,3] <- mean(results[,5])
  results.summary[count,4] <- mean(results[,8])
  results.summary[count,5] <- mean(results[,9])
  results.summary[count,6] <- mean(results[,12])
  results.summary[count,7] <- mean(results[,13])
  results.summary[count,8] <- mean(results[,16])
  results.summary[count,9] <- mean(results[,17])
  results.summary[count,10] <- mean(results[,20])
  results.summary[count,11] <- mean(results[,21])
  results.summary[count,12] <- mean(results[,24])
  results.summary[count,13] <- mean(results[,25])
  results.summary[count,14] <- mean(results[,28])
  results.summary[count,15] <- mean(results[,29])
  
  count <- count + 1
  print(paste("Done with N.pop =",j))
  
}


#test.name <- paste0("~/SN_Bootstrap_Sim_Results_hw.increase=",ci.hw.increase,".csv")
write.csv(results.summary, file = "betatest_bootstrap_b2.csv")

