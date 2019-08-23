library(tidyverse)
library(ggplot2)
library(sn)

sdprior <- read.csv("./Bayes/bayes_changing_n.pop_and_sdprior.csv") #N.pop.range 3:20, muprior = 70, sdprior 1:10
muprior <- read.csv("./Bayes/bayes_changing_n.pop_and_muprior.csv") #N.pop.range 3:20, muprior.range 10:90 by 5, sdprior = 10
musdprior <- read.csv("./Bayes/bayes_changing_n.pop_muprior_and_sdprior.csv") #N.pop.range 3:20, muprior.range 10:90 by 5, sdprior.range 1:10
test <- read.csv("./Bayes/bayes_test1.csv") #N.pop.range 3:20, muprior.range 10:80 by 10, sdprior.range 10:20 by 5
test2 <- read.csv("./Bayes/bayes_test2.csv") #N.pop.range 3:20, muprior.range 10:70 by 10, sdprior.range 20:30 by 5
test3 <- read.csv("./Bayes/bayes_test3.csv") #N.pop.range 3:20, muprior.range 10:90 by 5, sdprior.range 0:30 by 5
skewsdprior <- read.csv("./Bayes/bayes_changingN.pop_skew_and_sdprior.csv") #N.pop range 3:20, muprior set to 70, sdprior.range 5:30 by 5, skew = (-.9, -.5, .5, and .9)

#Establishing the true.mean for comparison


#Adding "cover" factor to data frames

musdprior <- musdprior %>%
  mutate(cover = ifelse(true.mean >= lowerbound.post & true.mean <= upperbound.post, 1,0))

sdprior %>% 
  mutate(cover = ifelse(true.mean >= lowerbound.post & true.mean <= upperbound.post, 1,0))



#Cleans up test3 dataframe to elimnate muprior and sdprior combinations outside the bouunds of 0:100
test3 <- test3 %>% 
  mutate(cover = ifelse(true.mean >= lowerbound.post & true.mean <= upperbound.post, 1, 0)) %>% 
  filter(muprior-sdprior >= 0) %>% 
  filter(muprior + sdprior <= 100) %>% 
  filter(sdprior != 0)

#Playing around with plots to try and figure out how to look at the data

test3 %>% filter(cover == 1) %>% #only want to look at results that covered true.mean
  ggplot(aes(x = width.ci.post, y = mean.post)) +
  geom_point(aes(color = muprior)) +
  geom_hline(yintercept = true.mean) +
  facet_wrap(~ sdprior, labeller = label_both)

#Generated plot that was kept and named (Bayes Width vs Muprior by SDprior)
test3 %>% filter(cover == 1) %>% #only want to look at results that covered true.mean
  ggplot(aes(x = width.ci.post, y = muprior)) +
  geom_point(aes(color = N.pop)) +
  geom_hline(yintercept = true.mean) +
  facet_wrap(~ sdprior, labeller = label_both)


###########################################################################################
#######################Changing Skew, sdprior, and N.pop ##################################
##############################################################################################

skewsdprior %>% filter(cover == 1) %>% 
  ggplot (aes(x = width.ci.post, y = N.pop)) +
  geom_point(aes(color = as.factor(skew))) +
  geom_vline(xintercept = 6.2) +
  facet_wrap(~ sdprior, labeller = label_both)





test %>% mutate(cover = ifelse(true.mean >= lowerbound.post & true.mean <= upperbound.post, 1, 0)) %>% 
  filter(cover ==1) %>% 
  #filter(muprior <= 50) %>% 
  ggplot(aes(x = width.ci.post, y = mean.post)) +
  geom_point(aes(color = muprior, size = N.pop)) + 
  geom_hline(yintercept = true.mean) +
  facet_wrap(~ sdprior, labeller = label_both)

test2 %>% mutate(cover = ifelse(true.mean >= lowerbound.post & true.mean <= upperbound.post, 1, 0)) %>% 
  filter(cover ==1) %>%
  filter(muprior - sdprior >= 0) %>% 
  ggplot(aes(x = width.ci.post, y = mean.post)) +
  geom_point(aes(color = muprior, size = N.pop)) + 
  geom_hline(yintercept = true.mean) +
  facet_wrap(~ sdprior, labeller = label_both)





#########Looking at How Changing SD of prior distribution affects aspects of the posterior

sdprior %>% 
  mutate(upper = mean.post+sd.post) %>% 
  mutate(lower = mean.post - sd.post) %>% 
  ggplot(aes(x = factor(sdprior)))+
  geom_boxplot(aes(lower = lower, middle = mean.post, upper = upper, ymax = upperbound.post, ymin = lowerbound.post), stat = "identity") +
  facet_wrap(~ N.pop)

sdprior %>% 
  ggplot(aes(x = sdprior, y = mean.post))+
  geom_point() +
  geom_errorbar(aes(ymin = lowerbound.post, ymax = upperbound.post)) +
  facet_wrap(~ N.pop)

#########Looking at How Changing Mu of prior distribution affects aspects of the posterior
muprior %>% 
  ggplot(aes(x = muprior, y = mean.post)) +
  geom_point() +
  geom_errorbar(aes(ymin = lowerbound.post, ymax = upperbound.post)) +
  geom_hline(yintercept = true.mean) +
  facet_wrap(~ N.pop)

####Looking at all Changes#################
musdprior %>% filter(cover == 1) %>% filter (width.ci.post <= 6) %>% 
  ggplot(aes(x = sd.post, y = mean.post)) +
  geom_point(aes(color = muprior)) +
  ##geom_errorbar(aes(ymin = lowerbound.post, ymax = upperbound.post)) +
  geom_hline(yintercept = true.mean) +
  facet_wrap(~ sdprior, labeller = label_both)
