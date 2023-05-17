############################################################################
# modified from Chapter 17 from Schaub and Kery 2021
# and Zhao, Q., Royle, J.A. and Boomer, G.S., 2017. Spatially explicit dynamic N-mixture models. Population Ecology, 59(4), pp.293-300.
# age structured pop using telemetry & count data 
# estimates proportion of unreported harvest data from catch and survey data
############################################################################

library(jagsUI)

######################
# simulated count data
######################

n.sites = 100 # scalar, number of habitat patches
n.years = 30 # scalar, number of years
#nreps = 3 # scalar, number of surveys within a primary sampling period (needed to estimate observation error)
#adj =  matrix(c(rep(rep(c(0,0), 15), 15), rep(rep(c(1,0), 15), 15)), nrow = nsite) # a matrix[nsite, nsite] describing the adjacency among habitat patches
#nadj = c(rep(1, 15), rep(14, 15)) # a vector[nsite] containing the number of adjacent patches for each patch

# storage arrays
N.A = array(dim = c(n.sites, n.years)) # an array[snite, nyear, nreps] containing the true counts of adults
N.J = array(dim = c(n.sites, n.years)) # an array[snite, nyear, nreps] containing the true counts of juveniles
SA = array(dim = c(n.sites, n.years)) # an array[snite, nyear, nreps] containing the surviving adults from previous year
SJ = array(dim = c(n.sites, n.years)) # an array[snite, nyear, nreps] containing the surviving juveniles from previous year
H = array(dim = c(n.sites, n.years)) # an array[snite, nyear, nreps] containing the number of fished animals
R = array(dim = c(n.sites, n.years)) # an array[snite, nyear, nreps] containing the number of new recruits
M = array(dim = c(n.sites, n.years)) # an array[snite, nyear, nreps] containing the number of juveniles that mature

lambda0_a = 100       # mean number of adults in year 1 
lambda0_j = 250       # mean number of juveniles in year 1
phi.a = 0.6           # adult survival
phi.j = 0.4           # juvenile survival
p.mat = 0.2           # probability of maturity
h = 0.4               # legal fishing pressure
ih = 0.4              # unreported fishing pressure 

#assuming Beverton Holt recruitment, these are general parameters I adjusted by eye
# increase beta to reduce carrying capacity (stochastic beta)
alpha = 800
beta = 100
#adults = seq(0, 400, by = 25)
#R = (alpha*adults)/(beta+adults)
#plot(adults, R)

# create simulated data
for(s in 1:n.sites){
  
  # year 1
  N.A[s,1] = rpois(1, lambda0_a) # if pop sizes across sites are overdispersed, can't use poisson
  N.J[s,1] = rpois(1, lambda0_j)
  
  # year 2
  SA[s,1] = N.A[s,1] * phi.a 
  SJ[s,1] = N.J[s,1] * max(rnorm(1, phi.j, 0.05),0)
  
  M[s,1] = SJ[s,1] * p.mat
  
  N.A[s,2] = as.integer((SA[s,1] + M[s,1]))
  N.J[s,2] = as.integer((SJ[s,1] - M[s,1]))
  
  for(t in 3:n.years){
    
    SA[s,t-1] = N.A[s,t-1] * phi.a 
    H[s,t-1] = as.integer((SA[s,t-1] * (h + ih)))
    SJ[s,t-1] = N.J[s,t-1] * max(rnorm(1, phi.j, 0.05),0)  # improve stochastic functions
    
    R[s,t-2] = (alpha*N.A[s,t-2])/(beta+N.A[s,t-2])
    M[s,t-1] = SJ[s,t-1] * p.mat
    
    N.A[s,t] = max(as.integer((SA[s,t-1] + M[s,t-1] - H[s,t-1])),1)
    N.J[s,t] = max(as.integer((SJ[s,t-1] + R[s,t-2] - M[s,t-1])),1)
    
  }
}


# simulate observation error of counts (but to accurately estimate, you need repeated counts within a year)
yA = matrix(round(rpois(length(N.A),N.A)), nrow = 100)     # observed counts of adults
yJ = matrix(round(rpois(length(N.J),N.J)), nrow = 100)     # observed counts of juveniles

plot(N.A, yA) # observed vs true counts of adults
plot(N.J, yJ) # observed vs true counts of juveniles

# simulate observation error of catch (but to accurately estimate, you need surveys that tell you compliance vs abundance)
C.true = H                                              # observed catch with unreported fishing
C = round(H * (h/(ih+h)))                            # observed catch with no unreported fishing
#C[,,2] = round(H * (1 - 0.2))                      # observed catch with constant unreported fishing (r = 0.2)
#C[,,3] = round(H * (1 - N.A/max(N.A, na.rm = T)))  # observed catch with dynamic unreported fishing (r = (h/N.A)/(max(h/N.A)))


############
# Model Runs
############

#Bundle data
jags.data <- list(lambda0.a = lambda0_a, lambda0.j = lambda0_j,
                  yA = yA, yJ = yJ,                             # counts of adults and juveniles
                  p.mat = p.mat,                                # probability of maturity (might be possible to estimate)
                  C=C, h = h,                                   # catch data
                  n.years=ncol(yA), n.sites = nrow(yA),         # number of years and sites
                  alpha = alpha, beta = beta)                       # BH parameters 


#MCMC settings
ni <- 10000; nb <- 5000; nc <- 3; nt <-50; na <- 5000


########################################
# JAGS model assuming no illegal harvest

#Initial values
inits <- function() {list(phi.a = 0.6, mu.phi.j = 0.4)}  # guess of natural mortality

#parameters monitored
parameters <- c("phi.a", "mu.phi.j", "sd.j",
                "N.A", "N.J")

#call JAGS from R
out1= jags(data = jags.data, inits = inits, parameters.to.save = parameters, model.file = "model2_no_ih.txt", 
           n.iter=ni, n.burnin = nb, n.chains= nc,
           n.thin= nt, n.adapt = na, parallel = TRUE)

save(out1, file = "no_ih_2.RData")



#####################################
# JAGS model file WITHOUT survey data, but assuming some unreported fishing

#Initial values
inits <- function() {list(phi.a = 0.6, mu.phi.j = 0.4, ih = 0.2)}  # guess of natural mortality


#parameters monitored
parameters <- c("phi.a", "mu.phi.j", "sd.j",
                "N.A", "N.J",
                "ih")

#call JAGS from R
out2= jags(data = jags.data, inits = inits, parameters.to.save = parameters, model.file = "model2_ih_uninformed.txt", 
           n.iter=ni, n.burnin = nb, n.chains= nc,
           n.thin= nt, n.adapt = na, parallel = TRUE)

save(out2, file = "ih_uninformed_2.RData")


###########################################################################
# JAGS model file WITH survey data that informs rates of unreported fishing

#Initial values
inits <- function() {list(phi.a = 0.6, mu.phi.j = 0.4, mu.ih = 5)}  # guess of natural mortality

#parameters monitored
parameters <- c("phi.a", "mu.phi.j", "sd.j",
                "N.A", "N.J",
                "mu.ih")

#call JAGS from R
out3= jags(data = jags.data, inits = inits, parameters.to.save = parameters, model.file = "model2_ih_informed.txt", 
           n.iter=ni, n.burnin = nb, n.chains= nc,
           n.thin= nt, n.adapt = na, parallel = TRUE)

save(out3, file = "ih_informed_2.RData")

