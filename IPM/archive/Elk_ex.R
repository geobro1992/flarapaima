#Ch. 17 IPM book- elk example of an age structured pop'n using telemetry & count data in conjunction with harvest data
library(IPMbook); library(jagsUI)
data(elk)
str(elk)
#Define the priors with our "guesses" about population parameters 
#start with estimate of survival and fecundity 
s= 0.8
f= 0.25 

# matrix population model (transition matrix a)
A= matrix(0, ncol=17, nrow=17) # elk age cap at 17 w/ one year age classes
A[1,2:17] <- f
for(a in 2:17){
  A[a,a-1] <- s
}

#compute stable age distribution (right eigenvector of A)
z= which.max(Re(eigen(A)$values))
revec= Re(eigen(A)$vectors[,z])

#using report rate (hunter sightings) and hunting mortality, estimate popn size 
#Population size in first age class of first year 
r=0.5 #guess of reporting rate
h=0.1 #guess of hunting mortality
N1= elk$C[1,1] / (h*r)

#compute age specific population sizes in first year
n= N1 * revec / revec[1]

#Bundle data and produce overview 
jags.data <- with(elk, list(C=C, a=H[1,], b=H[2,], R=R, total=colSums(R), n.years=ncol(C),
                            n.age=nrow(C), n=n))

#write JAGS model file
cat(file="model1.txt", " 
model {
#priors and linear models
for (t in 1: (n.years-1)) {
  f[t] ~ dunif(0,1)
}
for (t in 1:n.years) {
  s[t] ~ dunif(0,1)
  h[t] ~ dunif(0,1)
  r[t] ~ dunif(0,1)
}
#age at harvest data (state space model)
#model for the initial population size: poisson priors
for (a in 1:n.age) {
  N[a,1] ~ dpois(n[a])
}

#Process model over time: our model of population dynamics
for (t in 1: (n.years-1)) {
N[1, t+1] ~ dpois((Ntot[t] - N[1,t]) * f[t])
for (a in 2:n.age){
N[a, t+1] ~ dbin( (1-h[t]) * s[t], N[a-1,t])
} #a
} #t

#derived quantity: total year-specific population size 
for (t in 1:n.years){
  Ntot[t] <- sum(N[,t])
}
#observation model 
for (t in 1:n.years){
  for (a in 1:n.age){
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
inits <- function() {list(s= runif(jags.data$n.years, 0.8, 1), h= runif(jags.data$n.years, 0.05, 0.15))}

#parameters monitored
parameters <- c("h", "s", "f", "r", "Ntot", "N")

#MCMC settings
ni <- 250000; nb <- 50000; nc <- 3; nt <-20; na <- 5000

#call JAGS from R
out1= jags(jags.data, inits, parameters, "model1.txt", n.iter=ni, n.burnin = nb, n.chains= nc,
           n.thin= nt, n.adapt = na, parallel = TRUE)
traceplot(out1)
mean(out1$sims.list$r)
