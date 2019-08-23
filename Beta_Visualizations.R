library(ggplot2)
library(tidyverse)
library(gridExtra)

beta <- read.csv("betatest_bootstrap_b2.csv")

#Making Data TIDY

cover <- beta %>% select("N.pop", starts_with("cover"))
names(cover) <- substring(names(cover), 7)
colnames(cover)[1] <- "N.pop"
cover <- cover %>% gather("perc", "perc.exp", "clt", "tdist", "bca", "bca.exp", "bootstrapT",
                          key = "CI.type", value = "cover.percentage")

width <- beta %>% select("N.pop", starts_with("width"))
names(width) <- substring(names(width), 7)
colnames(width)[1] <- "N.pop"
width <- width %>% gather("perc", "perc.exp", "clt","tdist", "bca", "bca.exp", "bootstrapT",
                          key = "CI.type", value = "width")

tidybeta <- cover %>% right_join(width, by = c("N.pop", "CI.type"))


####GENERATE PLOTS####
p1 <- tidybeta %>%  filter(N.pop ==5) %>% 
  ggplot(aes(x = width, y = cover.percentage)) +
  geom_point(aes(color = CI.type)) + 
  geom_hline(yintercept = 0.95, size = 1) +
  ggtitle("N.pop = 5")

p2 <- tidybeta %>%  filter(N.pop ==10) %>% 
  ggplot(aes(x = width, y = cover.percentage)) +
  geom_point(aes(color = CI.type)) +
  geom_hline(yintercept = 0.95, size = 1) +
  ggtitle("N.pop = 10")

p3<- tidybeta %>%  filter(N.pop ==15) %>% 
  ggplot(aes(x = width, y = cover.percentage)) +
  geom_point(aes(color = CI.type)) +
  geom_hline(yintercept = 0.95, size = 1) +
  ggtitle("N.pop = 15")

grid.arrange(p1, p2, p3)


#What happens as you increase the population? Sorted for each CI type?
#N.pop is filtered greater than 4, because the bootstrapT for small samples generates a huge width for CI
tidybeta %>% filter(N.pop > 5) %>% 
  ggplot(aes(x = width, y = cover.percentage)) +
  geom_point(aes(color = N.pop)) +
  scale_color_gradient(low = "red", high = "green") +
  geom_hline(yintercept = 0.95, size = 1) +
  facet_wrap(~ CI.type)


#What happens if you change the CI type, based on population?
tidybeta %>% filter(N.pop > 5) %>% 
  ggplot(aes(x = width, y = cover.percentage)) +
  geom_point(aes(color = CI.type)) +
  geom_hline(yintercept = 0.95, size = 1) +
  facet_wrap(~N.pop)


#Overall visualization of different aspects, combined
#N.pop is filtered to be greater than 5 (3, 4, and 5 for BootstrapT are on the extremes)
tidybeta %>% filter(N.pop >= 5) %>% 
  ggplot(aes(x = width, y = cover.percentage)) +
  geom_point(aes(color = CI.type, size = N.pop)) +
  geom_hline(yintercept = 0.95, size = 1)






  