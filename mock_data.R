###################
# mock arapima data
###################

# libraries
library(dplyr)
library(ggplot2)

#######
# data
#######

#40 lakes
lake.id <- rep(1:40)

#assign the status for each lake using rbinom() function 
# 0 for fished and 1 for protected or vice versa
status <- c(rep(0, 20), rep(1,20))

#Assign 20 lakes to whitewater and the other 20 to blackwater
lake.type <- rep(c("Blackwater","Whitewater"), times=20)

#abundance in year 1 and year 2
n1 <- rnbinom(40, 1, 0.02)

#Factor by which to change the abundance, per status and water type
sfactor  <- c(
  rep(c(0.3,0.9), times = 10),
  rep(c(1.1,1.5), times = 10)
)


#Put all variables data together
arapima <- data.frame(lake.id, status, lake.type, n1, sfactor)

#Apply treatment effect per status and water type
arapima <- arapima %>% 
  mutate(n2 = n1*sfactor + sample(-3:3, 40, replace = T)) %>%
  mutate(n2 = ifelse(n2 < 0, 0, n2)) %>%
  mutate(lambda = 1+(n2-n1)/n1) %>%
  mutate(status = ifelse(status == 0, "Fished", "Protected"))


########
# stats
########

#compute descriptive statistics per status and water
stats <- arapima %>%
  group_by(status, lake.type) %>%
  summarise(
    count = n(),
    mean_n1 = mean(n1,na.rm=TRUE),
    se_n1 = sd(n1, na.rm = TRUE)/sqrt(count),
    mean_n2 = mean(n2,na.rm=TRUE),
    se_n2 = sd(n2, na.rm = TRUE)/sqrt(count),
    mean_lambda = mean(lambda,na.rm=TRUE),
    se_lambda = sd(lambda, na.rm = TRUE)/sqrt(count),
    ci95lower = mean_lambda - se_lambda*1.96,
    ci95upper = mean_lambda + se_lambda*1.96
  )

stats



#forest plot
p1 <- ggplot(stats, aes(x=mean_lambda, y=interaction(status, lake.type), color=status, shape=lake.type)) +
  geom_errorbar(aes(xmin = ci95lower, xmax = ci95upper), width = 0.5) +
  geom_point(size = 4) + 
  labs(x="Population Growth Rate", y = "") +
  xlim(0, 2) + 
  scale_color_manual(values = c("#ff6d6d", "#4e8ca5")) +
  scale_shape_manual(values=c(15,16)) +
  theme_classic() +
  theme(legend.position = "none",
        #remove y-axis line and ticks
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()) +
  geom_vline(xintercept = 1, linetype = "longdash") 

p1


