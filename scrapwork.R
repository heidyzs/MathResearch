dat<-list(N=length(y),
          y=y,
          sigma=sd(y))
fit <- stan(file = 'stanmod.stan', data = dat)
summary(fit, pars = c("mu"), probs = c(0.05, 0.95))$summary[[4]]
summary(fit, pars = c("mu"), probs = c(0.05, 0.95))$summary[[5]]


cp <- c(mean=70, s.d.=15, gamma1=-.5)
dp <- cp2dp(cp, family="SN")

library(sn) 

x<-seq(0,100,1)

trunc.dsn<-function(x){
  trunc.x<-1/(psn(100,dp=dp)-psn(0,dp=dp))*dsn(x,dp=dp)
  return(trunc.x)
}


my.df<- data.frame(x=x,y=trunc.dsn(x))
my.df %>% ggplot(aes(x=x,y=y))+geom_line(lwd=2)+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())



samp.mu<-sum(x*trunc.dsn(x))
samp.mu
