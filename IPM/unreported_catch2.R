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

lambda0_a = 500       # mean number of adults in year 1 
lambda0_j = 500       # mean number of juveniles in year 1
phi.a = 0.6           # adult survival
phi.j = 0.6           # juvenile survival
p.mat = 0.2           # probability of maturity
h = 0.2               # legal fishing pressure
ih = 55               # unreported catch 

#assuming Beverton Holt recruitment, these are general parameters I adjusted by eye
# increase beta to reduce carrying capacity (stochastic beta)
alpha = 400 
beta = 10
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
  SJ[s,1] = N.J[s,1] * phi.j
  
  M[s,1] = SJ[s,1] * p.mat
  
  N.A[s,2] = as.integer((SA[s,1] + M[s,1]))
  N.J[s,2] = as.integer((SJ[s,1] - M[s,1]))
  
  for(t in 3:n.years){
    
    SA[s,t-1] = N.A[s,t-1] * phi.a 
    H[s,t-1] = as.integer((N.A[s,t-1] * h) + ih)
    SJ[s,t-1] = N.J[s,t-1] * max(runif(1, phi.j-0.1, phi.j+0.1),0)  # improve stochastic functions
    
    R[s,t-2] = (alpha*N.A[s,t-2])/(beta+N.A[s,t-2])
    M[s,t-1] = SJ[s,t-1] * p.mat
    
    N.A[s,t] = max(as.integer((SA[s,t-1] + M[s,t-1] - H[s,t-1])),1)
    N.J[s,t] = max(as.integer((SJ[s,t-1] + R[s,t-2] - M[s,t-1])),1)
    
  }
}

# simulate observation error of counts (but to accurately estimate, you need repeated counts within a year)
e = 0.1                                          # observation error in counts
yA = round(rnorm(length(N.A),0,N.A*e) + N.A)     # observed counts of adults
yJ = round(rnorm(length(N.J),0,N.J*e) + N.J)     # observed counts of juveniles


# simulate observation error of catch (but to accurately estimate, you need surveys that tell you compliance vs abundance)
C = array(dim = c(n.sites, n.years, 3)) # an array[snite, nyear, nscenarios] to hold observed catch (true catch - unreported catch)
C[,,1] = H                                              # observed catch with unreported fishing
C[,,2] = H - ih                                         # observed catch with no unreported fishing
#C[,,2] = round(H * (1 - 0.2))                      # observed catch with constant unreported fishing (r = 0.2)
#C[,,3] = round(H * (1 - N.A/max(N.A, na.rm = T)))  # observed catch with dynamic unreported fishing (r = (h/N.A)/(max(h/N.A)))



# make fake survey data (to estimate rate of illegal fishing)
surveys = matrix(c(sample(40:80, size = 10, replace = T),
                   sample(20:40, size = 10, replace = T)),
                 ncol = 10, byrow = T)

colnames(surveys) = 2011:2020                                 # years
rownames(surveys) = c("Fishers surveyes", "Registered catch") # number of fishers vs number that report catch


# make fake telemetry data (to help separate natural mortality and fishing mortality)
#telemetry = matrix(c(sample(2:3, size = 10, replace = T),
#                     sample(2:3, size = 10, replace = T),
#                     sample(20:30, size = 10, replace = T)),
#                   ncol = 10, byrow = T)

#colnames(telemetry) = 2011:2020
#rownames(telemetry) = c("fished", "died naturally", "survived")



#Bundle data and produce overview 
jags.data <- list(yA = yA, yJ = yJ,                             # counts of adults and juveniles
                  p.mat = p.mat,                                # probability of maturity (might be possible to estimate)
                  C=C,                                          # catch data
                  a=surveys[1,], b=surveys[2,],                 # survey data
                  n.years=ncol(yA), n.sites = nrow(yA),         # number of years and sites
                  alpha = 400, beta = 10)                       # BH parameters 



#####################################
# JAGS model file WITHOUT survey data


cat(file="model2.txt", " 
    model {
    #priors
    lambda0.a ~ dgamma(0.001, 0.001)    # initial pop size of adults
    lambda0.j ~ dgamma(0.001, 0.001)    # initial pop size of juveniles
    phi.j ~ dunif(0, 1)             # natural juvenile mortality
    phi.a ~ dunif(0, 1)             # natural adult mortality  
    h ~ dunif(0,0.5)                  # fishing mortality
    e ~ dgamma(1, 1)                  # observation error
    # model for the initial population size: poisson priors

for (s in 1:n.sites){  
    # year 1
  N.A[s,1] ~ dpois(lambda0.a) 
  N.J[s,1] ~ dpois(lambda0.j)
  
  # year 2
  SA[s,1] ~ dbin(phi.a, N.A[s,1])
  SJ[s,1] ~ dbin(phi.j, N.J[s,1])
  
  M[s,1] ~ dbin(p.mat, SJ[s,1])
  
  N.A[s,2] <- (SA[s,1] + M[s,1])
  N.J[s,2] <- (SJ[s,1] - M[s,1])
  
  for(t in 3:n.years){
    
    SA[s,t-1] ~ dbin((1-h)*phi.a, N.A[s,t-1]) 
    SJ[s,t-1] ~ dbin(phi.j, N.J[s,t-1])  # improve stochastic functions
    
    R[s,t-2] ~ dpois((alpha*N.A[s,t-2])/(beta+N.A[s,t-2]))
    M[s,t-1] ~ dbin(p.mat, SJ[s,t-1])
    
    N.A[s,t] <- max((SA[s,t-1] + M[s,t-1]),1)
    N.J[s,t] <- max((SJ[s,t-1] + R[s,t-2] - M[s,t-1]),1)
    } # t
} #s
    
    #observation model 
    for (t in 3:n.years){
    for (s in 1:n.sites){
    C[s,t,2] ~ dbin(h, N.A[s,t])
    } #s
    } #t
    
    # observation model
    for (s in 1:n.sites) {
    for (t in 1:n.years) {
    yA[s,t] ~ dnorm(N.A[s,t], 1/pow(e*N.A[s,t],2))
    yJ[s,t] ~ dnorm(N.J[s,t], 1/pow(e*N.J[s,t],2))
    } # t
    } # s
}  
    ")

#Initial values
inits <- function() {list(h = 0.2,  # guess of fishing mortality
                          lambda0.a = 500, lambda0.j = 500,          # guess of mean pop size
                          phi.a = 0.7, phi.j = 0.6,                  # guess of natural mortality
                          e = 0.2)}                                # guess of recruitment

#parameters monitored
parameters <- c("phi.a", "phi.j", 
                "lambda0.a", "lambda0.j",
                "h",
                "e",
                "N.A", "N.J")

#MCMC settings
ni <- 30000; nb <- 5000; nc <- 3; nt <-200; na <- 5000

#call JAGS from R
out1= jags(data = jags.data, inits = inits, parameters.to.save = parameters, model.file = "model2.txt", 
           n.iter=ni, n.burnin = nb, n.chains= nc,
           n.thin= nt, n.adapt = na, parallel = TRUE)

save(out1, file = "model2_C2.RData")

##########
# OUTPUT
##########
#load("no_survey_C1.RData")

traceplot(out1)

# posterior means
la = out1$sims.list$lambda0.a # true value 500
lj = out1$sims.list$lambda0.j # true value 500
phij = out1$sims.list$phi.j # true value 0.6
phia = out1$sims.list$phi.a # true value 0.6
hest = out1$sims.list$h # true value 0.2
eest = out1$sims.list$e # true value 0.1

# storage vectors
NApred = matrix(NA, 30,375)
NJpred = matrix(NA, 30,375)
SApred = matrix(NA, 30,375)
SJpred = matrix(NA, 30,375)
Mpred = matrix(NA, 30,375)
Hpred = matrix(NA, 30,375)
Rpred = matrix(NA, 30,375)

# estimate population sizes
# year 1
NApred[1,] = rpois(375, la) 
NJpred[1,] = rpois(375, lj)

# year 2
SApred[1,] = NApred[1] * phia 
SJpred[1,] = NJpred[1] * phij

Mpred[1,] = SJpred[1,] * p.mat

NApred[2,] = (SApred[1,] + Mpred[1,])
NJpred[2,] = (SJpred[1,] - Mpred[1,])

for(t in 3:n.years){
  
  SApred[t-1,] = NApred[t-1,] * phia 
  Hpred[t-1,] = (NApred[t-1,] * hest)
  SJpred[t-1,] = NJpred[t-1,] * runif(375, phij-0.1, phij+0.1)
  
  Rpred[t-2,] = (alpha*NApred[t-2,])/(beta+NApred[t-2,])
  Mpred[t-1,] = SJpred[t-1,] * p.mat
  
  NApred[t,] = (SApred[t-1,] + Mpred[t-1,] - Hpred[t-1,])
  NJpred[t,] = (SJpred[t-1,] + Rpred[t-2,] - Mpred[t-1,])
  
}

# plot real vs expected adults
par(mfrow = c(2,1))
plot(colSums(N.A[,-c(1:5)])/100, type = "l", ylim = c(0,50), axes = F, main = "Adults", xlab = "", ylab = "")
lines(apply(N.A[,-c(1:5)], 2, quantile, probs = 0.15, na.rm = TRUE), lty = "dashed")
lines(apply(N.A[,-c(1:5)], 2, quantile, probs = 0.85, na.rm = TRUE), lty = "dashed")

lines(rowSums(NApred[-c(1:5),])/375, col = "red")
lines(apply(NApred[-c(1:5),], 1, quantile, probs = 0.025, na.rm = TRUE), lty = "dashed", col = "red")
lines(apply(NApred[-c(1:5),], 1, quantile, probs = 0.975, na.rm = TRUE), lty = "dashed", col = "red")

axis(1)
axis(2)

legend(20, 50, legend=c("True", "Predicted"),
       col=c("black", "red"), lty=1, cex=0.8)

# plot real vs expected juveniles
plot(colSums(N.J[,-c(1:5)])/100, type = "l", ylim = c(0,700), axes = F, main = "Juveniles", xlab = "Year", ylab = "")
lines(apply(N.J[,-c(1:5)], 2, quantile, probs = 0.15, na.rm = TRUE), lty = "dashed")
lines(apply(N.J[,-c(1:5)], 2, quantile, probs = 0.85, na.rm = TRUE), lty = "dashed")

lines(rowSums(NJpred[-c(1:5),])/375, col = "red")
lines(apply(NJpred[-c(1:5),], 1, quantile, probs = 0.025, na.rm = TRUE), lty = "dashed", col = "red")
lines(apply(NJpred[-c(1:5),], 1, quantile, probs = 0.975, na.rm = TRUE), lty = "dashed", col = "red")

axis(1)
axis(2)


###################
# demographic model WITH SURVEY DATA (untested)
###################

#write JAGS model file with survey data / unreported catch
cat(file="model1.txt", " 
    model {
    #priors
    lambda0.a ~ dgamma(0.01, 0.01)  # initial pop size of adults
    lambda0.j ~ dgamma(0.01, 0.01)  # initial pop size of juveniles
    phi.j ~ dunif(0, 1)             # natural juvenile mortality
    phi.a ~ dunif(0, 1)             # natural adult mortality  
    gamma ~ dgamma(1, 1)      # recruitment 
    h ~ dunif(0,1)                  # fishing mortality
    e ~ dgamma(1, 1)
    for(s in 1:n.sites) {
    for (t in 1:n.years) {
    r[s,t] ~ dunif(0,1)      # proportion of legal fishing
    }
    }
    #model for the initial population size: poisson priors
    for(s in 1:nsite) {
    
    N.A[s,1] ~ dpois(lambda0.a)
    N.J[s,1] ~ dpois(lambda0.j)
    
    for(t in 2:nyear) {
    
    SA[s,t-1]~dbin((1-h) * phi.a, N.A[s,t-1])
    SJ[s,t-1]~dbin(phi.j, N.J[s,t-1])
    
    REC[s,t-1] ~ dpois(gamma * N.A[s,t-1])
    M[s,t-1] ~ dbin(p.mat, SJ[s,t-1])
    
    N.A[s,t] <- SA[s,t-1] + M[s,t-1]                  # +I[i,t-1]-E[i,t-1]
    N.J[s,t] <- SJ[s,t-1] + (REC[i,t-1]) - M[s,t-1]   # +I[i,t-1]-E[i,t-1]
    
    } # t
    
    } # s

    #observation model 
    for (t in 1:n.years){
    for (s in 1:n.sites){
    C[s,t] ~ dbin(h * r[s,t], N.A[s,t])
    } #a
    } #t
    #hunter survey data (logistic regression)
    for (t in 1:n.years) {
    for (s in 1:n.sites){
    b[s,t] ~ dbin(r[s,t], a[s,t])
    }
    }
    #radio tracking data (multinomial)
    R[,t] ~ dmulti(prt, n.tracked)
    prt[1] <- h
    prt[2] <- (1-h) * (1-phi.a)
    prt[3] <- (1-h) * phi.a
    # observation model
    
    for (s in 1:nsite) {
    
    for (t in 1:nyear) {
    
    yA[s,t] ~ dnorm(N.A[s,t], 1/(pow(e*N.A[s,t],2)))
    yJ[s,t] ~ dnorm(N.J[s,t], 1/(pow(e*N.J[s,t],2))) 
    } # t
    
    } # s
    
    }  
    ")

#Initial values 
inits <- function() {list(h = runif(jags.data$n.years, 0.05, 0.15),  # guess of fishing mortality
                          r = runif(jags.data$n.years, 0.45, 0.55),  # guess of reporting rate (proportion of legal fishing)
                          lambda0.a = 100, lambda0.j = 100,          # guess of mean pop size
                          phi.a = 0.8, phi.j = 0.3,                  # guess of natural mortality
                          gamma = 1)}                                # guess of recruitment

#parameters monitored
parameters <- c("phi.a", "phi.j", "h", "r", "lambda0.a", "lambda0.j", "gamma")

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

