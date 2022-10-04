# code from:
# Zhao, Q., Royle, J.A. and Boomer, G.S., 2017. Spatially explicit dynamic N-mixture models. Population Ecology, 59(4), pp.293-300.

# data to provide
nsite = 30 # scalar, number of habitat patches
nyear = 10 # scalar, number of years
nreps = 3 # scalar, number of surveys within a primary sampling period
#adj =  matrix(c(rep(rep(c(0,0), 15), 15), rep(rep(c(1,0), 15), 15)), nrow = nsite) # a matrix[nsite, nsite] describing the adjacency among habitat patches
#nadj = c(rep(1, 15), rep(14, 15)) # a vector[nsite] containing the number of adjacent patches for each patch
y = array(dim = c(nsite, nyear, nreps)) # an array[snite, nyear, nreps] containing the observed counts

surv = 0.8
fec = 0.4
p = 0.8

for(i in 1:nsite){
  for(j in 1:nreps){
    
    y[i,1,j] = sample(150:160, 1, replace = TRUE)
    
    for(k in 2:nyear){
      
      y[i,k,j] = as.integer((y[i,k-1,j] * surv + y[i,k-1,j] * fec))
      
    }
  }
}

y = round(y*p)
head(y)
# parameters to estimate.
# lambda0: scalar, mean abundance in the first year
# omega: scalar, survival rate
# kappa: scalar, the rate of emigration
# gamma: scalar, reproductive rate
# pobs: scalar, detection probability

# generated quantities
# Ilambda[nsite, nyear-1]: a matrix containing the expectation of immigration
# S[nsite, nyear-1]: a matrix containing the number of survived individuals
# E[nsite, nyear-1]: a matrix containing the number of emigrated individuals
# I[nsite, nyear-1]: a matrix containing the number of immigrated individuals
# R[nsite, nyear-1]: a matrix containing the number of reproduced individuals.
# N[nsite, nyear]: a matrix containing the true local population size


cat(file = "Nmix.jags", "
model {
  
  # priors
  lambda0~dunif(150, 160)
  omega~dunif(0.7, 0.9)
  gamma~dunif(0.2, 0.5)
  pobs~dunif(0.6, 0.9)
  
  # process model
  
  for(i in 1:nsite) {
    
    N[i,1]~dpois(lambda0)
    
    for(t in 2:nyear) {
      

      S[i,t-1]~dbin(omega, N[i,t-1])
      
      R[i,t-1]~dpois(gamma*N[i,t-1])
      
      N[i,t] <- S[i,t-1]+R[i,t-1] # +I[i,t-1]-E[i,t-1]
      
    } # t
    
  } # i
  
  # observation model
  
  for (i in 1:nsite) {
    
    for (t in 1:nyear) {
      
      for (j in 1:nreps) {
        
        y[i,t,j] ~ dbin(pobs, N[i,t])
        
      } # j
      
    } # t
    
  } # i
  
} # model
    ")





###############################################################################################
# run model and store results
###############################################################################################

# bundle data
bugs.data <- list(nsite = nsite, nyear = nyear, nreps = nreps, y = y)
parameters = c("lambda0", "omega", "gamma", "pobs", "S", "R", "N")


# Initial values
inits <- function(){list(lambda0 = 155, omega = 0.8, gamma = 0.4, pobs = 0.7)}


# storage matrices
Ilambda = matrix(nrow = nsite, ncol = nyear-1) # a matrix containing the expectation of immigration
S = matrix(nrow = nsite, ncol = nyear-1) # a matrix containing the number of survived individuals
E = matrix(nrow = nsite, ncol = nyear-1) # a matrix containing the number of emigrated individuals
I = matrix(nrow = nsite, ncol = nyear-1) # a matrix containing the number of immigrated individuals
R = matrix(nrow = nsite, ncol = nyear-1) # a matrix containing the number of reproduced individuals.
N = matrix(nrow = nsite, ncol = nyear) # a matrix containing the true local population size

# sampling parameters
nc = 3
nt = 5
ni = 1000
nb = 500

# Call JAGS from R (jagsUI)
require(jagsUI)
m.Nmix <- jags(data = bugs.data, inits = inits, parameters, "Nmix.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

hist(m.Nmix$sims.list$omega, breaks = 50)
hist(m.Nmix$sims.list$gamma, breaks = 50)
hist(m.Nmix$sims.list$pobs, breaks = 50)

