library(sn) # sampling from the skew-normal distribution is accomplished via the 'sn' package in R
library(e1071) # needed to calculate sample skewness
library(resample) # a package that allows for the expanded version of the percentile bootstrap
library(foreach) # allows the code to run in a vectorized fashion [Note: There is still a pesky inner for loop.]
library(doParallel) # allows parallel execution

### set the number of cores
registerDoParallel(cores=4)

###Functions###
func<-function(x){
  x*dsn(x,dp[1],dp[2],dp[3])*(1/(ub-lb))
}

### GLOBAL DECLARATIONS

set.seed(8)
N.trials <- 1000 # number of Monte Carlo trials to run (e.g., to estimate the CI coverage probabilities); 1000 default
N.bootstrap <- 1000 # size of bootstrap sample for trial i; 1000 default
N.pop.range<-seq(from=3, to=30, by = 1) # N.pop to try
sn.skew.range <-seq(from=-0.9, to=0.9, by = 0.05) # sn_skew to try

# used for experimentation on obtaining the desired coverage probability by reducing alpha
# e.g., if equal to 0, then you get a 95% CI; if equal to 0.005, then you get a 96% CI
ci.hw.increase <- 0

sn_mean <- 68 # the mean of the parent SN; based on Table 8.6 in Sauro & Lewis (2016, p. 205)
sn_stdev <- 12.5 # the standard deviation of the parent SN; based on Table 8.6 in Sauro & Lewis (2016, p. 205)

### START EXPERIMENTATION

start_time <- Sys.time()



results.summary2 <- foreach(k=sn.skew.range) %:% 
  foreach(j=N.pop.range, .combine='rbind', .packages=c('sn','resample')) %dopar% {
    
    sn_skew <- k  
    
    # See Azzalini's "How to sample from the SN and related distributions when we want to fix skewness and other cumulants"
    # cp is the centered parameterization of the SN, where s.d. => standard deviation, and gamma1 => skewness
    # Note that |gamma1| < 0.99527
    cp <- c(mean=sn_mean, s.d.=sn_stdev, gamma1=sn_skew)
    dp <- cp2dp(cp, family="SN") # the cp2dp function converts the centered parameterization to the direct parameterization 
    LB.for.prob<-psn(0, dp=dp) # finds the probability associated with the minimum SUS
    UB.for.prob<-psn(100, dp=dp) # finds the probability associated with the maximum SUS
    
    N.pop <- j # sample from parent population (e.g., skew-normal to replicate the respondent's SUS scores)
    
    func<-function(x){
      x*dsn(x,dp[1],dp[2],dp[3])*(1/(UB.for.prob-LB.for.prob))
    }
    
    true.mean<-integrate(func,0,100)[1]$value
    
    results <- data.frame(matrix(0, ncol = 26, nrow = N.trials)) # build a preallocated data frame to hold the results
    
    for (i in 1:N.trials) {
      
      # generate N.pop random variates from the SN parent within SUS' bounds [0,100]
      y<-qsn(runif(n=N.pop,min=LB.for.prob,max=UB.for.prob), dp=dp, solver="RFB")
      
      # generate bootstrap sample from the random variates
      b <- bootstrap(y, mean(y), R = N.bootstrap)  
      
      # generate 95% CIs from the bootstrap sample
      ci_percentile <- CI.percentile(b, probs = c(0.025-ci.hw.increase, 0.975+ci.hw.increase), expand = FALSE)
      ci_percentile_expanded <- CI.percentile(b, probs = c(0.025-ci.hw.increase, 0.975+ci.hw.increase), expand = TRUE)
      lb.clt<-mean(y)-1.96*sd(y)/sqrt(length(y))
      ub.clt<-mean(y)+1.96*sd(y)/sqrt(length(y))
      lb.tdist <- mean(y) - qt(0.975, df = length(y)-1)*sd(y)/sqrt(length(y))
      ub.tdist <- mean(y) + qt(0.975, df = length(y)-1)*sd(y)/sqrt(length(y))
      ci_bca <- CI.bca(b, probs = c(0.025-ci.hw.increase, 0.975+ci.hw.increase), expand = FALSE)
      ci_bca_expanded <- CI.bca(b, probs = c(0.025-ci.hw.increase, 0.975+ci.hw.increase), expand = TRUE)
      results[i,1] <- i
      if ((true.mean <= ci_percentile[1])) {results[i,2] <- 1} #Misses and is too low
      if ((true.mean >= ci_percentile[2])) {results[i,3] <- 1} #Misses and is too high
      if ((true.mean >= ci_percentile[1]) && (true.mean <= ci_percentile[2])) {results[i,4] <- 1} # Acutally inside the interval
      results[i,5] <- ci_percentile[2] - ci_percentile[1] #Width of the interval
      if ((true.mean <= ci_percentile_expanded[1])) {results[i,6] <- 1}
      if ((true.mean >= ci_percentile_expanded[2])) {results[i,7] <- 1}
      if ((true.mean >= ci_percentile_expanded[1]) && (true.mean <= ci_percentile_expanded[2])) {results[i,8] <- 1}
      results[i,9] <- ci_percentile_expanded[2] - ci_percentile_expanded[1]
      if ((true.mean <= lb.clt)){results[i,10]<- 1}
      if ((true.mean >= ub.clt)) {results[i,11]<- 1}
      if ((true.mean >= lb.clt) && (true.mean <=ub.clt)) {results[i,12] <- 1}
      results[i,13] <- ub.clt-lb.clt
      if ((true.mean <= lb.tdist)) {results[i,14]<- 1}
      if ((true.mean >= ub.tdist)) {results[i,15]<- 1}
      if ((true.mean >= lb.tdist) && (true.mean <=ub.tdist)) {results[i,16] <- 1}
      results[i,17] <- ub.tdist-lb.tdist
      if ((true.mean <= ci_bca[1])){results[i,18] <- 1}
      if ((true.mean >= ci_bca[2])) {results[i,19] <- 1}
      if ((true.mean >= ci_bca[1]) && (true.mean <= ci_bca[2])) {results[i,20] <- 1}
      results[i,21] <- ci_bca[2] - ci_bca[1]
      if ((true.mean <= ci_bca_expanded[1])) {results[i,22] <- 1}
      if ((true.mean >= ci_bca_expanded[2])){results[i,23] <- 1}
      if ((true.mean >= ci_bca_expanded[1]) && (true.mean <= ci_bca_expanded[2])) {results[i,24] <- 1}
      results[i,25] <- ci_bca_expanded[2] - ci_bca_expanded[1]
    } # end (i in N.trials)
    
    # summarize results
    
    c(j,k,mean(results[,2]),mean(results[,3]),mean(results[,4]),mean(results[,5]),mean(results[,6]),
      mean(results[,7]), mean(results[,8]), mean(results[,9]), mean(results[,10]), mean(results[,11]),
      mean(results[,12]), mean(results[,13]), mean(results[,14]), mean(results[,15]), mean(results[,16]), 
      mean(results[,17]), mean(results[,18]), mean(results[,19]), mean(results[,20]), mean(results[,21]),
      mean(results[,22]), mean(results[,23]), mean(results[,24]), mean(results[,25]))
    #print(paste("Done with N.pop =",j,"Skew =",k))
    
  } # end foreach (k in sn.skew.range) %:% foreach(j=N.pop.range, .combine='rbind')

end_time <- Sys.time()
end_time - start_time

results.summary2 <- as.data.frame(do.call(rbind, results.summary2))
colnames(results.summary2) <- c("N.pop","skew", "low.perc", "high.perc","cover.perc", "width.perc", 
                                "low.perc.exp", "high.perc.exp", "cover.perc.exp", "width.perc.exp",
                                "low.clt", "high.clt","cover.clt", "width.clt", 
                                "low.tdist", "high.tdist", "cover.tdist",  "width.tdist",
                                "low.bca", "high.bca","cover.bca", "width.bca",
                                "low.bca.exp", "high.bca.exp", "cover.bca.exp", "width.bca.exp")

# test.name <- paste0("~/Dabkowski/ORCEN Projects/NAG (Me)/Analyzing the SUS in R/Assessing Sample Size/SN_Bootstrap_Sim_Results_hw.increase=",ci.hw.increase,".csv")
# write.csv(results.summary, file = test.name)

write.csv(results.summary2, "results_0hw.csv")



##############################################
#Building individual summmary tables

#Builds a summary table for Percentages Below LB CI

for (x in sn.skew.range) {
  low_results <- data.frame(matrix(0, ncol = length(sn.skew.range +1), nrow = 6*length(N.pop.range))) # build a preallocated data frame to hold the results
  colnames(low_results)  <- sn.skew.range
  
  for (y in 6*length(N.pop.range)) {
    if (y <= length(N.pop.range)) {low_results[] <- results.summary2[,3]}
    if(y > length(N.pop.range) && y <= 2*length(N.pop.range)) {low_results[] <- results.summary2[,7]}
    if(y > 2*length(N.pop.range) && y <= 3*length(N.pop.range)) {low_results[] <- results.summary2[,11]}
    if(y > 3*length(N.pop.range) && y <= 4*length(N.pop.range)) {low_results[] <- results.summary2[,15]}
    if(y > 4*length(N.pop.range)&& y <= 5*length(N.pop.range)) {low_results[] <- results.summary2[,19]}
    if(y > 5*length(N.pop.range)) {low_results[] <- results.summary2[,23]}
  }
}

#Builds a summary table for Percentages Above UB CI
for (x in sn.skew.range) {
  high_results <- data.frame(matrix(0, ncol = length(sn.skew.range +1), nrow = 6*length(N.pop.range))) # build a preallocated data frame to hold the results
  colnames(high_results)  <- sn.skew.range
  
  for (y in 6*length(N.pop.range)) {
    if (y <= length(N.pop.range)) {high_results[] <- results.summary2[,4]}
    if(y > length(N.pop.range) && y <= 2*length(N.pop.range)) {high_results[] <- results.summary2[,8]}
    if(y > 2*length(N.pop.range) && y <= 3*length(N.pop.range)) {high_results[] <- results.summary2[,12]}
    if(y > 3*length(N.pop.range) && y <= 4*length(N.pop.range)) {high_results[] <- results.summary2[,16]}
    if(y > 4*length(N.pop.range)&& y <= 5*length(N.pop.range)) {high_results[] <- results.summary2[,20]}
    if(y > 5*length(N.pop.range)) {high_results[] <- results.summary2[,24]}
  }
}

#Builds a summary table for all cover percentages
for (x in sn.skew.range) {
  cover_results <- data.frame(matrix(0, ncol = length(sn.skew.range +1), nrow = 6*length(N.pop.range))) # build a preallocated data frame to hold the results
  colnames(cover_results)  <- sn.skew.range
  
  for (y in 6*length(N.pop.range)) {
    if (y <= length(N.pop.range)) {cover_results[] <- results.summary2[,5]}
    if(y > length(N.pop.range) && y <= 2*length(N.pop.range)) {cover_results[] <- results.summary2[,9]}
    if(y > 2*length(N.pop.range) && y <= 3*length(N.pop.range)) {cover_results[] <- results.summary2[,13]}
    if(y > 3*length(N.pop.range) && y <= 4*length(N.pop.range)) {cover_results[] <- results.summary2[,17]}
    if(y > 4*length(N.pop.range)&& y <= 5*length(N.pop.range)) {cover_results[] <- results.summary2[,21]}
    if(y > 5*length(N.pop.range)) {cover_results[] <- results.summary2[,25]}
  }
}


#Builds a summary table for all Widths for each type of CI
for (x in sn.skew.range) {
  width_results <- data.frame(matrix(0, ncol = length(sn.skew.range +1), nrow = 6*length(N.pop.range))) # build a preallocated data frame to hold the results
  colnames(width_results)  <- sn.skew.range
  
  for (y in 6*length(N.pop.range)) {
    if (y <= length(N.pop.range)) {width_results[] <- results.summary2[,6]}
    if(y > length(N.pop.range) && y <= 2*length(N.pop.range)) {width_results[] <- results.summary2[,10]}
    if(y > 2*length(N.pop.range) && y <= 3*length(N.pop.range)) {width_results[] <- results.summary2[,14]}
    if(y > 3*length(N.pop.range) && y <= 4*length(N.pop.range)) {width_results[] <- results.summary2[,18]}
    if(y > 4*length(N.pop.range)&& y <= 5*length(N.pop.range)) {width_results[] <- results.summary2[,22]}
    if(y > 5*length(N.pop.range)) {width_results[] <- results.summary2[,26]}
  }
}


write.csv(low_results, "skewnormal_low.csv")
write.csv(high_results, "skewnormal_high.csv")
write.csv(cover_results, "skewnormal_cover.csv")
write.csv(width_results, "skewnormal_width.csv")

##Extra Stuff trying to build summary tables##

#BUILD A SUMMARY TABLE FOR Bootstrap
for (x in sn.skew.range) {
  boot_perc_table <- data.frame(matrix(0, ncol = length(sn.skew.range+1), nrow = length(N.pop.range))) # build a preallocated data frame to hold the results
  colnames(boot_perc_table)  <- sn.skew.range
  rownames(boot_perc_table) <- N.pop.range
  
  for (y in N.pop.range) {
    boot_perc_table[] <- results.summary2[,6]
  }
}

write.csv(boot_perc_table, "dataframe.csv")


############################
#########################BOOOT EXP##################
#######################################################

for (x in sn.skew.range) {
  boot_exp_low <- data.frame(matrix(0, ncol = length(sn.skew.range+1), nrow = length(N.pop.range))) # build a preallocated data frame to hold the results
  colnames(boot_exp_low)  <- sn.skew.range
  rownames(boot_exp_low) <- N.pop.range
  
  for (y in N.pop.range) {
    boot_exp_low[] <- results.summary2[,23]
  }
}

for (x in sn.skew.range) {
  boot_exp_high <- data.frame(matrix(0, ncol = length(sn.skew.range+1), nrow = length(N.pop.range))) # build a preallocated data frame to hold the results
  colnames(boot_exp_high)  <- sn.skew.range
  rownames(boot_exp_high) <- N.pop.range
  
  for (y in N.pop.range) {
    boot_exp_high[] <- results.summary2[,24]
  }
}

for (x in sn.skew.range) {
  boot_exp_cover <- data.frame(matrix(0, ncol = length(sn.skew.range+1), nrow = length(N.pop.range))) # build a preallocated data frame to hold the results
  colnames(boot_exp_cover)  <- sn.skew.range
  rownames(boot_exp_cover) <- N.pop.range
  
  for (y in N.pop.range) {
    boot_exp_cover[] <- results.summary2[,25]
  }
}

for (x in sn.skew.range) {
  boot_exp_width <- data.frame(matrix(0, ncol = length(sn.skew.range+1), nrow = length(N.pop.range))) # build a preallocated data frame to hold the results
  colnames(boot_exp_width)  <- sn.skew.range
  rownames(boot_exp_width) <- N.pop.range
  
  for (y in N.pop.range) {
    boot_exp_width[] <- results.summary2[,26]
  }
}

write.csv(boot_exp_low, "boot_exp_low.csv")
write.csv(boot_exp_high, "boot_exp_high.csv")
write.csv(boot_exp_cover, "boot_exp_cover.csv")
write.csv(boot_exp_width, "boot_exp_width.csv")

#BUILD A SUMMARY TABLE FOR BOOTSTRAP EXPANDED

for (x in sn.skew.range) {
  boot_perc_exp_table <- data.frame(matrix(0, ncol = length(sn.skew.range), nrow = length(N.pop.range))) # build a preallocated data frame to hold the results
  colnames(boot_perc_exp_table)  <- sn.skew.range
  rownames(boot_perc_exp_table) <- N.pop.range
  
  for (y in N.pop.range) {
    boot_perc_exp_table[] <- results.summary2[,5]
  }
}

#BUILD A SUMMARY TABLE FOR CLT
for (x in sn.skew.range) {
  clt_table <- data.frame(matrix(0, ncol = length(sn.skew.range), nrow = length(N.pop.range))) # build a preallocated data frame to hold the results
  colnames(clt_table)  <- sn.skew.range
  rownames(clt_table) <- N.pop.range
  
  for (y in N.pop.range) {
    clt_table[] <- results.summary2[,7]
  }
}

#BUILD A SUMMARY TABLE FOR TDISTRIBUTION
for (x in sn.skew.range) {
  tdist_table <- data.frame(matrix(0, ncol = length(sn.skew.range), nrow = length(N.pop.range))) # build a preallocated data frame to hold the results
  colnames(tdist_table)  <- sn.skew.range
  rownames(tdist_table) <- N.pop.range
  
  for (y in N.pop.range) {
    tdist_table[] <- results.summary2[,9]
  }
}

#BUILD A SUMMMARY TABLE FOR BCa
for (x in sn.skew.range) {
  bca_table <- data.frame(matrix(0, ncol = length(sn.skew.range), nrow = length(N.pop.range))) # build a preallocated data frame to hold the results
  colnames(bca_table)  <- sn.skew.range
  rownames(bca_table) <- N.pop.range
  
  for (y in N.pop.range) {
    bca_table[] <- results.summary2[,11]
  }
}

#BUILD A SUMMARY TABLE FOR BCa EXPANDED
for (x in sn.skew.range) {
  bca_exp_table <- data.frame(matrix(0, ncol = length(sn.skew.range), nrow = length(N.pop.range))) # build a preallocated data frame to hold the results
  colnames(bca_exp_table)  <- sn.skew.range
  rownames(bca_exp_table) <- N.pop.range
  
  for (y in N.pop.range) {
    bca_exp_table[] <- results.summary2[,13]
  }
}



##########
########## extra stuff below
# plot the results to visualize the results
hist(b,prob=TRUE)
abline(v=sn_mean,lty=2,col="red") # the population mean
abline(v=c(ci_percentile[1],ci_percentile[2]),lty=1,col="blue")
abline(v=c(ci_percentile_expanded[1],ci_percentile_expanded[2]),lty=1,col="green")

# checks to ensure the sample cp is close to the parent
mean(y)
sd(y)
skewness(y)
