########################################
# JAGS model assuming no illegal harvest


cat(file="model2_no_ih.txt", " 
    model {
    #priors
 
    mu.phi.j ~ dunif(0, 1)             # natural juvenile mortality
    tau.j <-pow(sd.j, -2)
    sd.j ~dgamma(1,1)

    phi.a ~ dunif(0,1)               # natural adult mortality  
    

# juvenile survival random effect
for (s in 1:n.sites){  
  for(t in 1:(n.years-1)){
 phi.j[s,t] ~ dnorm(mu.phi.j, tau.j)
  }
}
    
# model for the initial population size: poisson priors

for (s in 1:n.sites){  
    # year 1
  N.A[s,1] ~ dpois(lambda0.a) 
  N.J[s,1] ~ dpois(lambda0.j)
  
  # year 2
  
  SA[s,1] ~ dbin(phi.a, N.A[s,1])
  SJ[s,1] ~ dbin(phi.j[s,1], N.J[s,1])
  
  M[s,1] ~ dbin(p.mat, SJ[s,1])
  
  N.A[s,2] <- (SA[s,1] + M[s,1])
  N.J[s,2] <- (SJ[s,1] - M[s,1])
  
  for(t in 3:n.years){
    
    SA[s,t-1] ~ dbin((1-h)*phi.a, N.A[s,t-1]) 
    SJ[s,t-1] ~ dbin(phi.j[s,t-1], N.J[s,t-1])  # improve stochastic functions
    
    R[s,t-2] ~ dpois((alpha*N.A[s,t-2])/(beta+N.A[s,t-2]))
    M[s,t-1] ~ dbin(p.mat, SJ[s,t-1])
    
    N.A[s,t] <- max((SA[s,t-1] + M[s,t-1]),1)
    N.J[s,t] <- max((SJ[s,t-1] + R[s,t-2] - M[s,t-1]),1)
    } # t
} #s
    
    #observation model 

    for (t in 3:n.years){
      for (s in 1:n.sites){
      
        C[s,t] ~ dbin(h, N.A[s,t])
        
      } #s
    } #t
    
    for (s in 1:n.sites) {
    for (t in 1:n.years) {

    yA[s,t] ~ dpois(N.A[s,t])
    yJ[s,t] ~ dpois(N.J[s,t])

    } # t
    } # s
}  
    ")


#####################################
# JAGS model file WITHOUT survey data, but assuming some unreported fishing

cat(file="model2_ih_uninformed.txt", " 
    model {
    #priors
 
    mu.phi.j ~ dunif(0, 1)             # natural juvenile mortality
    tau.j <-pow(sd.j, -2)
    sd.j ~dgamma(1,1)

    phi.a ~ dunif(0,1)               # natural adult mortality  
    
    ih ~ dunif(0,0.5)                # illegal fishing mortality


# juvenile survival random effect
for (s in 1:n.sites){  
  for(t in 1:(n.years-1)){
 phi.j[s,t] ~ dnorm(mu.phi.j, tau.j)
  }
}
    
# model for the initial population size: poisson priors

for (s in 1:n.sites){  
    # year 1
  N.A[s,1] ~ dpois(lambda0.a) 
  N.J[s,1] ~ dpois(lambda0.j)
  
  # year 2
  
  SA[s,1] ~ dbin(phi.a, N.A[s,1])
  SJ[s,1] ~ dbin(phi.j[s,1], N.J[s,1])
  
  M[s,1] ~ dbin(p.mat, SJ[s,1])
  
  N.A[s,2] <- (SA[s,1] + M[s,1])
  N.J[s,2] <- (SJ[s,1] - M[s,1])
  
  for(t in 3:n.years){
    
    SA[s,t-1] ~ dbin((1-(h+ih))*phi.a, N.A[s,t-1]) 
    SJ[s,t-1] ~ dbin(phi.j[s,t-1], N.J[s,t-1])
    
    R[s,t-2] ~ dpois((alpha*N.A[s,t-2])/(beta+N.A[s,t-2]))
    M[s,t-1] ~ dbin(p.mat, SJ[s,t-1])
    
    N.A[s,t] <- max((SA[s,t-1] + M[s,t-1]),1)
    N.J[s,t] <- max((SJ[s,t-1] + R[s,t-2] - M[s,t-1]),1)
    } # t
} #s
    
    #observation model 

    for (t in 3:n.years){
      for (s in 1:n.sites){
      
        C[s,t] ~ dbin(h, N.A[s,t])
        
      } #s
    } #t
    
    for (s in 1:n.sites) {
    for (t in 1:n.years) {

    yA[s,t] ~ dpois(N.A[s,t])
    yJ[s,t] ~ dpois(N.J[s,t])

    } # t
    } # s
}  
    ")




#####################################
# JAGS model file WITH survey data that informs rates of unreported fishing

cat(file="model2_ih_informed.txt", " 
    model {
    #priors
 
    mu.phi.j ~ dunif(0, 1)             # natural juvenile mortality
    tau.j <-pow(sd.j, -2)
    sd.j ~dgamma(1,1)

    phi.a ~ dunif(0,1)               # natural adult mortality  
    
    ih ~ dunif(0,0.5)I(0.15,0.25)    # illegal fishing mortality

# juvenile survival random effect
for (s in 1:n.sites){  
  for(t in 1:(n.years-1)){
 phi.j[s,t] ~ dnorm(mu.phi.j, tau.j)
  }
}
    
# model for the initial population size: poisson priors

for (s in 1:n.sites){  
    # year 1
  N.A[s,1] ~ dpois(lambda0.a) 
  N.J[s,1] ~ dpois(lambda0.j)
  
  # year 2
  
  SA[s,1] ~ dbin(phi.a, N.A[s,1])
  SJ[s,1] ~ dbin(phi.j[s,1], N.J[s,1])
  
  M[s,1] ~ dbin(p.mat, SJ[s,1])
  
  N.A[s,2] <- (SA[s,1] + M[s,1])
  N.J[s,2] <- (SJ[s,1] - M[s,1])
  
  for(t in 3:n.years){
    
    SA[s,t-1] ~ dbin((1-(h+ih))*phi.a, N.A[s,t-1]) 
    SJ[s,t-1] ~ dbin(phi.j[s,t-1], N.J[s,t-1])  # improve stochastic functions
    
    R[s,t-2] ~ dpois((alpha*N.A[s,t-2])/(beta+N.A[s,t-2]))
    M[s,t-1] ~ dbin(p.mat, SJ[s,t-1])
    
    N.A[s,t] <- max((SA[s,t-1] + M[s,t-1]),1)
    N.J[s,t] <- max((SJ[s,t-1] + R[s,t-2] - M[s,t-1]),1)
    } # t
} #s
    
    #observation model 

    for (t in 3:n.years){
      for (s in 1:n.sites){
      
        C[s,t] ~ dbin(h, N.A[s,t])
        
      } #s
    } #t
    
    for (s in 1:n.sites) {
    for (t in 1:n.years) {

    yA[s,t] ~ dpois(N.A[s,t])
    yJ[s,t] ~ dpois(N.J[s,t])

    } # t
    } # s
}  
    ")
