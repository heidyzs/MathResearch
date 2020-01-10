library(tidyverse)
x<-seq(0,100,.1)

y<- dsn(x, xi=83.78,omega=31.03,alpha=-1.81)
y1<-func(x)

my.df<-data.frame(x=x,y=y1)

my.df %>% ggplot(aes(x=x,y=y1))+geom_line(lwd=2)+
  theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.x=element_blank(),axis.title.y=element_blank())
