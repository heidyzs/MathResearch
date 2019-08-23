library(sn) # sampling from the skew-normal distribution is accomplished via the 'sn' package in R
library(e1071) # needed to calculate sample skewness
library(resample) # a package that allows for the expanded version of the percentile bootstrap


#Functions##
func<-function(x){
  x*dsn(x,dp[1],dp[2],dp[3])*(1/(UB.for.prob-LB.for.prob))
}

### GLOBAL DECLARATIONS

set.seed(8)
N.trials <- 100 # number of Monte Carlo trials to run (e.g., to estimate the CI coverage probabilities); 1000 default
N.bootstrap <- 1000 # size of bootstrap sample for trial i; 1000 default
N.pop.range<-seq(from=3, to=30, by = 5) # N.pop to try (was3 to 30)
sn.skew.range <-seq(from=-0.9, to=0.9, by = 0.1) # sn_skew to try (was -0.9 to 0.9)

# used for experimentation on obtaining the desired coverage probability by reducing alpha
# e.g., if equal to 0, then you get a 95% CI; if equal to 0.005, then you get a 96% CI
ci.hw.increase <- 0 

sn_mean <- 68 # the mean of the parent SN; based on Table 8.6 in Sauro & Lewis (2016, p. 205)
sn_stdev <- 12.5 # the standard deviation of the parent SN; based on Table 8.6 in Sauro & Lewis (2016, p. 205)

### START EXPERIMENTATION

results.summary <- data.frame(matrix(0, ncol = 26, nrow = (length(N.pop.range)*length(sn.skew.range)))) # build a preallocated data frame to hold the results
colnames(results.summary) <- c("N.pop","skew", "low.perc", "high.perc","cover.perc", "width.perc", 
                               "low.perc.exp", "high.perc.exp", "cover.perc.exp", "width.perc.exp",
                               "low.clt", "high.clt","cover.clt", "width.clt", 
                               "low.tdist", "high.tdist", "cover.tdist",  "width.tdist",
                               "low.bca", "high.bca","cover.bca", "width.bca",
                               "low.bca.exp", "high.bca.exp", "cover.bca.exp", "width.bca.exp")
count <- 1

for(k in sn.skew.range) {
  
  sn_skew <- k 
  
  # See Azzalini's "How to sample from the SN and related distributions when we want to fix skewness and other cumulants"
  # cp is the centered parameterization of the SN, where s.d. => standard deviation, and gamma1 => skewness
  # Note that |gamma1| < 0.99527
  cp <- c(mean=sn_mean, s.d.=sn_stdev, gamma1=sn_skew)
  dp <- cp2dp(cp, family="SN") # the cp2dp function converts the centered parameterization to the direct parameterization 
  LB.for.prob<-psn(0, dp=dp) # finds the probability associated with the minimum SUS
  UB.for.prob<-psn(100, dp=dp) # finds the probability associated with the maximum SUS
  
  func<-function(x){
    x*dsn(x,dp[1],dp[2],dp[3])*(1/(UB.for.prob-LB.for.prob))
  }
  
  true.mean<-integrate(func,0,100)[1]$value
  
  for (j in N.pop.range) {
    N.pop <- j # sample from parent population (e.g., skew-normal to replicate the respondent's SUS scores)
    
    results <- data.frame(matrix(0, ncol = 25, nrow = N.trials)) # build a preallocated data frame to hold the results
    x <- c("trial", "LB.perc", "RB.perc", "cover.perc","width.perc",
           "LB.perc.exp", "RB.perc.exp", "cover.perc.exp","width.perc.exp", "LB.clt", "RB.clt", "cover.perc.clt", "width.clt", 
           "LB.tdist","RB.tdist","cover.perc.tdist","width.tdist",
           "LB.bca", "RB.bca", "cover.bca", "width.bca", "LB.bca.exp", "RB.bca,exp", "width.bca.exp", "cover.bca.exp")
    colnames(results) <- x
    
    # start the simulation
    for (i in 1:N.trials) {
      # generate N.pop random variates from the SN parent within SUS' bounds [0,100]
      y<-qsn(runif(n=N.pop,min=LB.for.prob,max=UB.for.prob), dp=dp, solver="RFB")
      
      # generate bootstrap sample from the random variates
      b <- bootstrap(y, mean(y), R = N.bootstrap)  
      #b2 <- bootstrap(y, c(mean = mean(y), sd = sd(y)), R = N.bootstrap)
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
    } # end for (i in 1:N.trials)
    
    # summarize results
    results.summary[count,1] <- j
    results.summary[count,2] <- k
    results.summary[count,3] <- mean(results[,2])
    results.summary[count,4] <- mean(results[,3])
    results.summary[count,5] <- mean(results[,4])
    results.summary[count,6] <- mean(results[,5])
    results.summary[count,7] <- mean(results[,6])
    results.summary[count,8] <- mean(results[,7])
    results.summary[count,9] <- mean(results[,8])
    results.summary[count,10] <- mean(results[,9])
    results.summary[count,11] <- mean(results[,10])
    results.summary[count,12] <- mean(results[,11])
    results.summary[count,13] <- mean(results[,12])
    results.summary[count,14] <- mean(results[,13])
    results.summary[count,15] <- mean(results[,14])
    results.summary[count,16] <- mean(results[,15])
    results.summary[count,17] <- mean(results[,16])
    results.summary[count,18] <- mean(results[,17])
    results.summary[count,19] <- mean(results[,18])
    results.summary[count,20] <- mean(results[,19])
    results.summary[count,21] <- mean(results[,20])
    results.summary[count,22] <- mean(results[,21])
    results.summary[count,23] <- mean(results[,22])
    results.summary[count,24] <- mean(results[,23])
    results.summary[count,25] <- mean(results[,24])
    results.summary[count,26] <- mean(results[,25])
    
    count <- count + 1
    print(paste("Done with N.pop =",j,"Skew =",k))
    
  } # end for (j in N.pop.range)
} # end for (k in sn.skew.range)

# test.name <- paste0("~/Dabkowski/ORCEN Projects/NAG (Me)/Analyzing the SUS in R/Assessing Sample Size/SN_Bootstrap_Sim_Results_hw.increase=",ci.hw.increase,".csv")
# write.csv(results.summary, file = test.name)


write.csv(results.summary, "results_summary_test2.csv")

#BUILD A SUMMARY TABLE FOR Bootstrap
for (x in sn.skew.range) {
  boot_perc_table_low <- data.frame(matrix(0, ncol = length(sn.skew.range +1), nrow = 6*length(N.pop.range))) # build a preallocated data frame to hold the results
  colnames(boot_perc_table_low)  <- sn.skew.range
  #rownames(boot_perc_table_low) <- rep(N.pop.range, times = 6)
  
  for (y in 6*length(N.pop.range)) {
    if (y <= length(N.pop.range)) {
      boot_perc_table_low[] <- results.summary[,3]
      
    }
    if(y > length(N.pop.range) && y <= 2*length(N.pop.range)) {
      boot_perc_table_low[] <- results.summary[,7]
      
    }
    if(y > 2*length(N.pop.range) && y <= 3*length(N.pop.range)) {
      boot_perc_table_low[] <- results.summary[,11]
      
    }
    if(y > 3*length(N.pop.range) && y <= 4*length(N.pop.range)) {
      boot_perc_table_low[] <- results.summary[,15]

    }
    if(y > 4*length(N.pop.range)&& y <= 5*length(N.pop.range)) {
      boot_perc_table_low[] <- results.summary[,19]
     
    }
    if(y > 5*length(N.pop.range)) {
      boot_perc_table_low[] <- results.summary[,23]
      
    }
  }
}


for (x in sn.skew.range) {
  boot_perc_table_high <- data.frame(matrix(0, ncol = length(sn.skew.range+1), nrow = length(N.pop.range))) # build a preallocated data frame to hold the results
  colnames(boot_perc_table_high)  <- sn.skew.range
  rownames(boot_perc_table_high) <- N.pop.range
  
  for (y in N.pop.range) {
    boot_perc_table_high[] <- results.summary[,4]
  }
}

for (x in sn.skew.range) {
  boot_perc_table <- data.frame(matrix(0, ncol = length(sn.skew.range+1), nrow = length(N.pop.range))) # build a preallocated data frame to hold the results
  colnames(boot_perc_table)  <- sn.skew.range
  rownames(boot_perc_table) <- N.pop.range
  
  for (y in N.pop.range) {
    boot_perc_table[] <- results.summary[,5]
  }
}

#BUILD A SUMMARY TABLE FOR BOOTSTRAP EXPANDED

for (x in sn.skew.range) {
  boot_perc_exp_table <- data.frame(matrix(0, ncol = length(sn.skew.range), nrow = length(N.pop.range))) # build a preallocated data frame to hold the results
  colnames(boot_perc_exp_table)  <- sn.skew.range
  rownames(boot_perc_exp_table) <- N.pop.range
  
  for (y in N.pop.range) {
    boot_perc_exp_table[] <- results.summary[,9]
  }
}

#BUILD A SUMMARY TABLE FOR CLT
for (x in sn.skew.range) {
  clt_table <- data.frame(matrix(0, ncol = length(sn.skew.range), nrow = length(N.pop.range))) # build a preallocated data frame to hold the results
  colnames(clt_table)  <- sn.skew.range
  rownames(clt_table) <- N.pop.range
  
  for (y in N.pop.range) {
    clt_table[] <- results.summary[,7]
  }
}

#BUILD A SUMMARY TABLE FOR TDISTRIBUTION
for (x in sn.skew.range) {
  tdist_table <- data.frame(matrix(0, ncol = length(sn.skew.range), nrow = length(N.pop.range))) # build a preallocated data frame to hold the results
  colnames(tdist_table)  <- sn.skew.range
  rownames(tdist_table) <- N.pop.range
  
  for (y in N.pop.range) {
    tdist_table[] <- results.summary[,9]
  }
}

for (x in sn.skew.range) {
  bca_table <- data.frame(matrix(0, ncol = length(sn.skew.range), nrow = length(N.pop.range))) # build a preallocated data frame to hold the results
  colnames(bca_table)  <- sn.skew.range
  rownames(bca_table) <- N.pop.range
  
  for (y in N.pop.range) {
    bca_table[] <- results.summary[,11]
  }
}

for (x in sn.skew.range) {
  bca_exp_table <- data.frame(matrix(0, ncol = length(sn.skew.range), nrow = length(N.pop.range))) # build a preallocated data frame to hold the results
  colnames(bca_exp_table)  <- sn.skew.range
  rownames(bca_exp_table) <- N.pop.range
  
  for (y in N.pop.range) {
    bca_exp_table[] <- results.summary[,13]
  }
}



write.csv(boot_perc_table, "skewnormal_boot.csv")
write.csv(boot_perc_exp_table, "skewnormal_boot_exp.csv")
write.csv(clt_table, "skewnormal_clt.csv")
write.csv(tdist_table, "skewnormal_tdist.csv")
write.csv(bca_table, "skewnormal_bca.csv")
write.csv(bca_exp_table, "skewnormal_bca_exp.csv")

#write.csv(boot_perc_table_low, "testtable.csv")


