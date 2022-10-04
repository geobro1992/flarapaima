###################
# mock arapima data
###################

# libraries
library(dplyr)
library(ggplot2)
library(MASS)

#######
# data
#######

#60 lakes
lake.id <- rep(1:60)

# lake area
area <-  round(rnorm(60, 200, 50),1)

#assign the status for each lake
status <- as.factor(c(rep("Commercial", 20), rep("Subsistance",20), rep("Protected",20)))

#abundance in year 1
n1 <- rnbinom(60, 1, 0.02)

#Factor by which to change the abundance, per status
sfactor  <- c(
  rep(c(0.2), times = 20),
  rep(c(1), times = 20),
  rep(c(4), times = 20)
)


#Put all variables data together
arapima <- data.frame(lake.id, area, status, n1, sfactor)

#Apply treatment effect per status and water type
arapima <- arapima %>% 
  mutate(n1 = n1*sfactor) %>%
  mutate(n1 = as.integer(n1*area/100))

########
# stats
########

#compute descriptive statistics per status and water
stats <- arapima %>%
  group_by(status) %>%
  summarise(
    count = n(),
    mean_n1 = mean(n1,na.rm=TRUE),
    se_n1 = sd(n1, na.rm = TRUE)/sqrt(count),
    ci95lower = mean_n1 - se_n1*1.96,
    ci95upper = mean_n1 + se_n1*1.96,
    mean_area = mean(area,na.rm=TRUE),
    se_area = sd(area, na.rm = TRUE)/sqrt(count),
    
  )

stats



#forest plot (pop size vs lake status)
p1 <- ggplot(stats, aes(x=mean_n1, y=status, color=status, shape=status)) +
  geom_errorbar(aes(xmin = ci95lower, xmax = ci95upper), width = 0.5) +
  geom_point(size = 4) + 
  labs(x="Population Size", y = "") +
  scale_color_manual(values = c("#ff6d6d", "#4e8ca5", "#66aa66", "#aa9966")) +
  scale_shape_manual(values=c(15,16, 17)) +
  theme_classic() +
  theme(legend.position = "none",
        #remove y-axis line and ticks
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank()) 
p1

#forest plot (pop size vs lake size)
p2 <- ggplot(arapima, aes(x=area, y=n1, color=status, shape=status)) +
  geom_point(size = 4) + 
  labs(x="Lake Area", y = "Pop Size") +
  scale_color_manual(values = c("#ff6d6d", "#4e8ca5", "#66aa66", "#aa9966")) +
  scale_shape_manual(values=c(15,16, 17)) +
  theme_classic() 

p2


##############
# nb model
m1 = glm.nb(n1 ~ status + area, data = arapima)
summary(m1)


# predictions
newdat <- data.frame(
  area = rep(seq(from = min(arapima$area), to = max(arapima$area), length.out = 100), 3),
  status = factor(rep(1:3, each = 100), levels = 1:3, labels =
                  levels(arapima$status)))

newdat <- cbind(newdat, predict(m1, newdat, type = "link", se.fit=TRUE))
newdat <- within(newdat, {
  N <- exp(fit)
  LL <- exp(fit -  se.fit)
  UL <- exp(fit +  se.fit)
})

ggplot(newdat, aes(area, N)) +
  geom_ribbon(aes(ymin = LL, ymax = UL, fill = status), alpha = .25) +
  geom_line(aes(colour = status), size = 2) +
  labs(x = "Lake Area", y = "Predicted Population Size")+
  scale_color_manual(values = c("#ff6d6d", "#4e8ca5", "#66aa66", "#aa9966")) +
  scale_shape_manual(values=c(15,16, 17)) +
  theme_classic() 
