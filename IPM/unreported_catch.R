############################################################################
# modified from Chapter 17 from Schaub and Kery 2021
#
# age structured pop using telemetry & count data 
# estimates proportion of unreported harvest data from catch and survey data
############################################################################

library(jagsUI)

######
# Data
######

# make fake catch data
catch = matrix(c(rep(0,20), 
                 sample(4:6, size = 10, replace = T),
                 sample(3:5, size = 10, replace = T),
                 sample(1:3, size = 10, replace = T),
                 sample(0:3, size = 10, replace = T),
                 sample(0:2, size = 10, replace = T),
                 sample(0:2, size = 10, replace = T)), 
               ncol = 10, byrow = T)

colnames(catch) = 2011:2020  # years
rownames(catch) = 1:8        # age of fish

# make fake survey data (to estimate rate of illegal fishing)
surveys = matrix(c(sample(40:80, size = 10, replace = T),
                   sample(20:40, size = 10, replace = T)),
                 ncol = 10, byrow = T)

colnames(surveys) = 2011:2020                                 # years
rownames(surveys) = c("Fishers surveyes", "Registered catch") # number of fishers vs number that report catch


# make fake telemetry data (to help separate natural mortlaity and fishing mortality)
telemetry = matrix(c(sample(0:2, size = 10, replace = T),
                     sample(0:2, size = 10, replace = T),
                     sample(20:25, size = 10, replace = T)),
                   ncol = 10, byrow = T)

colnames(telemetry) = 2011:2020
rownames(telemetry) = c("fished", "died naturally", "survived")

##########################################
# combine all 3 sources of data into a list
flounder = list(C = catch, H = surveys, R = telemetry)


###################
# demographic model
###################

#estimate of survival  
#M = c(1.341, 0.585, 0.441, 0.386, 0.359, 0.344, 0.336) # real data for flounder (need to modify model to be age-dependent survival)
#s= 1-exp(-M)
s = 0.4           # survival rate
a.mat = 2         # age at maturity
a.max = 8         # lifespan

# fecundity (need to modify model to be age-dependent fecundity)
f= 0.25 

# matrix population model (transition matrix a)
A= matrix(0, ncol=a.max, nrow=a.max) # flounder age cap at 8 w/ one year age classes

A[1,a.mat:a.max] <- f

for(a in 2:a.max){
  A[a,a-1] <- s
}

#compute stable age distribution (right eigenvector of A)
z= which.max(Re(eigen(A)$values))
revec= Re(eigen(A)$vectors[,z])


#Population size in first age class of first year 
N1= 200
#compute age specific population sizes in first year
n= N1 * revec / revec[1]

#Bundle data and produce overview 
jags.data <- with(flounder, list(C=C, a=H[1,], b=H[2,], R=R, total=colSums(R), n.years=ncol(C),
                            a.max=a.max, a.mat = a.mat, n=n))

#write JAGS model file
cat(file="model1.txt", " 
model {
#priors and linear models
for (t in 1: (n.years-1)) {
  f[t] ~ dunif(0,1)      # fecundity
}
for (t in 1:n.years) {
  s[t] ~ dunif(0,1)      # survival
  h[t] ~ dunif(0,1)      # hunting mortality
  r[t] ~ dunif(0,1)      # proportion of legal fishing
}

#age at harvest data (state space model)
#model for the initial population size: poisson priors
for (a in 1:a.max) {
  N[a,1] ~ dpois(n[a])
}

#Process model over time: our model of population dynamics
for (t in 1: (n.years-1)) {
N[1, t+1] ~ dpois((Ntot[t] - sum(N[1:(a.mat-1),t])) * f[t])

for (a in a.mat:a.max){
   N[a, t+1] ~ dbin( (1-h[t]) * s[t], N[a-1,t])
 } #a
} #t

#derived quantity: total year-specific population size 
for (t in 1:n.years){
  Ntot[t] <- sum(N[,t])
}

#observation model 
for (t in 1:n.years){
  for (a in a.mat:a.max){
    C[a,t] ~ dbin(h[t] * r[t], N[a,t])
  } #a
} #t
#hunter survey data (logistic regression)
for (t in 1:n.years) {
  b[t] ~ dbin(r[t], a[t])
}
#radio tracking data (multinomial)
for (t in 1:n.years) {
  R[,t] ~ dmulti(prt[,t], total[t])
  prt[1,t] <- h[t]
  prt[2,t] <- (1-h[t]) * (1-s[t])
  prt[3,t] <- (1-h[t]) * s[t]
}
}
")

#Initial values 
inits <- function() {list(s= runif(jags.data$n.years, 0.75, 0.85),  # guess of survival
                          h= runif(jags.data$n.years, 0.05, 0.15),  # guess of fishing mortality
                          r= runif(jags.data$n.years, 0.45, 0.55))} # guess of reporting rate (proportion of legal fishing)

#parameters monitored
parameters <- c("s", "h", "f", "r", "Ntot")

#MCMC settings
ni <- 250000; nb <- 50000; nc <- 3; nt <-200; na <- 5000

#call JAGS from R
out1= jags(data = jags.data, inits = inits, parameters.to.save = parameters, model.file = "model1.txt", 
           n.iter=ni, n.burnin = nb, n.chains= nc,
           n.thin= nt, n.adapt = na, parallel = TRUE)

traceplot(out1)
mean(out1$sims.list$r)
mean(out1$sims.list$s)
mean(out1$sims.list$h)
