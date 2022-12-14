#############
# code from:
# Zhao, Q., Royle, J.A. and Boomer, G.S., 2017. Spatially explicit dynamic N-mixture models. Population Ecology, 59(4), pp.293-300.
# and
# Zhao, Q. 2021. A simulation study of the age-structured spatially explicit dynamic N-mixture model.
###############

# data to provide
nsite = 1000 # scalar, number of habitat patches
nyear = 10 # scalar, number of years
nreps = 1 # scalar, number of surveys within a primary sampling period
#adj =  matrix(c(rep(rep(c(0,0), 15), 15), rep(rep(c(1,0), 15), 15)), nrow = nsite) # a matrix[nsite, nsite] describing the adjacency among habitat patches
#nadj = c(rep(1, 15), rep(14, 15)) # a vector[nsite] containing the number of adjacent patches for each patch

yA = array(dim = c(nsite, nyear, nreps)) # an array[snite, nyear, nreps] containing the observed counts
yJ = array(dim = c(nsite, nyear, nreps)) # an array[snite, nyear, nreps] containing the observed counts
SA = array(dim = c(nsite, nyear, nreps)) # an array[snite, nyear, nreps] containing the observed counts
SJ = array(dim = c(nsite, nyear, nreps)) # an array[snite, nyear, nreps] containing the observed counts
R = array(dim = c(nsite, nyear, nreps)) # an array[snite, nyear, nreps] containing the observed counts
M = array(dim = c(nsite, nyear, nreps)) # an array[snite, nyear, nreps] containing the observed counts

lambda0_a = 100
lambda0_j = 100
surv.a = 0.9
surv.j = 0.3
fec = 1.5
p.mat = 0.25
e = 0.1

for(i in 1:nsite){
  
  for(j in 1:nreps){
    
    yA[i,1,j] = rpois(1, lambda0_a)
    yJ[i,1,j] = rpois(1, lambda0_j)
    
    for(k in 2:nyear){
      
      SA[i,k-1, j] = yA[i,k-1, j] * surv.a
      SJ[i,k-1, j] = yJ[i,k-1, j] * surv.j
      
      R[i,k-1, j] = rpois(1, (fec*yA[i,k-1, j]))
      M[i,k-1, j] = SJ[i,k-1, j] * p.mat
      
      yA[i,k,j] = as.integer((SA[i,k-1,j] + M[i,k-1,j]))
      yJ[i,k,j] = as.integer((SJ[i,k-1,j] + R[i,k-1,j] - M[i,k-1,j]))
      
    }
  }
}


yA = round(rnorm(length(yA),0,yA*e) + yA)
yJ = round(rnorm(length(yJ),0,yJ*e) + yJ)
head(yA)
head(yJ)


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


cat(file = "Age_Structured_Nmix.jags", "
model {
  
  # priors
  lambda0.a ~ dgamma(0.01, 0.01)
  omega.a ~ dunif(0.5, 1)
  gamma ~ dgamma(0.01, 0.01)

  lambda0.j ~ dgamma(0.01, 0.01)
  omega.j ~ dunif(0, 0.5)

  # process model
  
  for(i in 1:nsite) {
    
    N.A[i,1] ~ dpois(lambda0.a)
    N.J[i,1] ~ dpois(lambda0.j)
    
    for(t in 2:nyear) {
      

      SA[i,t-1]~dbin(omega.a, N.A[i,t-1])
      SJ[i,t-1]~dbin(omega.j, N.J[i,t-1])
      
      R[i,t-1] ~ dpois(gamma * N.A[i,t-1])
      M[i,t-1] ~ dbin(p.mat, SJ[i,t-1])
      
      N.A[i,t] <- SA[i,t-1]+M[i,t-1]          # +I[i,t-1]-E[i,t-1]
      N.J[i,t] <- SJ[i,t-1]+(R[i,t-1])-M[i,t-1] # +I[i,t-1]-E[i,t-1]
      
    } # t
    
  } # i
  
  # observation model
  
  for (i in 1:nsite) {
    
    for (t in 1:nyear) {
      
      for (j in 1:nreps) {

        yA[i,t,j] ~ dnorm(N.A[i,t], 1/(0.1*N.A[i,t]*0.1*N.A[i,t]))
        yJ[i,t,j] ~ dnorm(N.J[i,t], 1/(0.1*N.J[i,t]*0.1*N.J[i,t])) 
        
      } # j
      
    } # t
    
  } # i

} # model
    ")





###############################################################################################
# run model and store results
###############################################################################################

# bundle data
bugs.data <- list(nsite = nsite, nyear = nyear, nreps = nreps, yA = yA, yJ = yJ, p.mat = p.mat)
parameters = c("lambda0.a", "lambda0.j", "omega.a", "omega.j", "gamma")


# Initial values
inits <- function(){list(lambda0.a = 100, lambda0.j = 100, omega.a = 0.5, omega.j = 0.5, gamma = 1)}

# sampling parameters
nc = 3
nt = 40
ni = 5000
nb = 1000

# Call JAGS from R (jagsUI)
require(jagsUI)
m.Nmix <- jags(data = bugs.data, inits = inits, parameters, "Age_Structured_Nmix.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)


#################
# posterior plots

par(mfrow = c(1,3))

#hist(m.Nmix$sims.list$lambda0.a, breaks = 50, xlim = c(0, 20)) # real value is 10
#abline(v = 10, col = "red", lwd = 2)

#hist(m.Nmix$sims.list$lambda0.j, breaks = 50, xlim = c(0, 20)) # real value is 10
#abline(v = 10, col = "red", lwd = 2)

hist(m.Nmix$sims.list$omega.a, breaks = 10, 
     main = "Adult Survival", xlim = c(0.5, 1)) # real value is 0.9
abline(v = 0.9, col = "red", lwd = 2)

hist(m.Nmix$sims.list$omega.j, breaks = 10, 
     main = "Juvenile Survival", xlim = c(0, 0.5)) # real value is 0.3
abline(v = 0.3, col = "red", lwd = 2)

hist(m.Nmix$sims.list$gamma, breaks = 10, 
     main = "Fecundity", xlim = c(1, 2)) # real value is 1.5
abline(v = 1.5, col = "red", lwd = 2)




#######################
# just 1st year


cat(file = "pois_mix.jags", "
model {
  
  # priors
  lambda0.a ~ dgamma(0.01, 0.01)
  lambda0.j ~ dgamma(0.01, 0.01)  

  pobs ~ dunif(0.8,1)
  
  # process model
  
  for(i in 1:nsite) {
    
    N.A[i,1] ~ dpois(lambda0.a)
    N.J[i,1] ~ dpois(lambda0.j)
    
  } # i
  
  # observation model
  
  for (i in 1:nsite) {
    
      for (j in 1:nreps) {
        
        yA[i,1,j] ~ dbin(pobs, N.A[i,1])
        yJ[i,1,j] ~ dbin(pobs, N.J[i,1])
        
      } # j
      
  } # i
  
} # model
    ")



# bundle data
bugs.data <- list(nsite = nsite, nreps = nreps, yA = yA, yJ = yJ)
parameters = c("lambda0.a", "lambda0.j", "pobs")


# Initial values
inits <- function(){list(lambda0.a = 20, lambda0.j = 20, pobs = 0.95)}

# sampling parameters
nc = 3
nt = 5
ni = 1000
nb = 500

# Call JAGS from R (jagsUI)
require(jagsUI)
m.Nmix <- jags(data = bugs.data, inits = inits, parameters, "pois_mix.jags", n.chains = nc, n.thin = nt, n.iter = ni, n.burnin = nb)

hist(m.Nmix$sims.list$lambda0.a, breaks = 50)
hist(m.Nmix$sims.list$lambda0.j, breaks = 50)
hist(m.Nmix$sims.list$pobs, breaks = 50)

