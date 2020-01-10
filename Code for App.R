library(rstan)
#Generates fake data
y<-runif(10,60,80)
dat<-list(N=length(y),
          y=y,
          sigma=sd(y))
fit <- stan(file = 'stanmod.stan', data = dat,refresh=0)

#Compiles the stan model in C++


#Once it's compiled all we need to run is:

##Get User Data from practitioner stored as user.dat


#Fit user dat to (precompiled) stan model
#new.fit <- stan(file = 'stanmod.stan', data = user.dat,refresh=0)

#Extract Credible Interval
#bayes.lb<-summary(fit, pars = c("mu"), probs = c(0.05, 0.95))$summary[[4]]
#bayes.ub<-summary(fit, pars = c("mu"), probs = c(0.05, 0.95))$summary[[5]]
#print(c(bayes.lb,bayes.ub))
