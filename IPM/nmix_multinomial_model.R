########################################
# multinomial illegal harvest uninformed

cat(file="ih_multi_baseline.txt", " 
    model {
    #priors

##########
# maybe try Dirichlet? - can inform each phi

phi_raw[1] <- psi * (1-F)
phi_raw[2] <- F
phi_raw[3] ~ dgamma(0.1, 0.1)

for (i in 1:3) {
  phi[i] <- phi_raw[i] / sum(phi_raw)
}

    
    psi ~ dnorm(0.8, 10)I(0,1)
    F ~ dnorm(0.6, 10)I(0,1)


    lambda ~ dnorm(20,10)
    gamma ~ dnorm(10,10)
    
    p ~ dnorm(0.9, 10)I(0,1)
    
    
# model for the initial population size: poisson priors

for (s in 1:n.sites){  

    # year 1
  N[s,1] ~ dpois(lambda) 

  for(t in 1:(n.years-1)){
    
    S[s,t+1] ~ dbin(phi[1], N[s,t])
    C[s,t+1] ~ dbin(phi[2], N[s,t]) 
    D[s,t+1] <- N[s,t] - (S[s,t+1] + C[s,t+1]) 
    
    R[s,t+1] ~ dpois(gamma)

    N[s,t+1] <- (S[s,t+1] + R[s,t+1])

    } # t
    

     for (t in 1:n.years) {
      for (j in 1:nreps) {

    y[s,t,j] ~ dbin(p, N[s,t])

      } # j
     } # t
    } # s
}  
    ")



########################################
# multinomial illegal harvest uninformed

cat(file="ih_multi_uninformed.txt", " 
    model {
    #priors

##########
# maybe try Dirichlet? - can inform each phi

phi_raw[1] <- psi * (1-F)
phi_raw[2] <- F * r
phi_raw[3] <- F * (1-r)
phi_raw[4] ~ dgamma(0.1, 0.1)

for (i in 1:4) {
  phi[i] <- phi_raw[i] / sum(phi_raw)
}

    
    psi ~ dgamma(0.1, 0.1)I(0,1)
    r ~ dgamma(0.1, 0.1)I(0,1)
    F ~ dgamma(0.1, 0.1)I(0,1)
    
    lambda ~ dnorm(20,10)
    gamma ~ dnorm(10,10)
    
    p ~ dnorm(0.9, 10)I(0,1)
    
    
# model for the initial population size: poisson priors

for (s in 1:n.sites){  

    # year 1
  N[s,1] ~ dpois(lambda) 

  for(t in 1:(n.years-1)){
    
    S[s,t+1] ~ dbin(phi[1], N[s,t])
    C[s,t+1] ~ dbin(phi[2], N[s,t]) 
    UC[s,t+1] ~ dbin(phi[3], N[s,t]) 
    D[s,t+1] <- N[s,t] - (S[s,t+1] + C[s,t+1] + UC[s, t+1]) 
    
    R[s,t+1] ~ dpois(gamma)

    N[s,t+1] <- (S[s,t+1] + R[s,t+1])

    } # t
    

     for (t in 1:n.years) {
      for (j in 1:nreps) {

    y[s,t,j] ~ dbin(p, N[s,t])

      } # j
     } # t
    } # s
}  
    ")



########################################
# multinomial illegal harvest informed


cat(file="ih_multi_informed.txt", " 
    model {
    #priors

##########
# maybe try Dirichlet? - can inform each phi
# use beta distribution for informative priors

phi_raw[1] <- psi * (1-F)
phi_raw[2] <- F * r
phi_raw[3] <- F * (1-r)
phi_raw[4] ~ dnorm(0.2, 100)I(0,1)

for (i in 1:4) {
  phi[i] <- phi_raw[i] / sum(phi_raw)
}

    
    psi ~ dnorm(0.8, 100)I(0,1)
    r ~ dnorm(0.8, 100)I(0,1)
    F ~ dnorm(0.6, 100)I(0,1)
    
    lambda ~ dnorm(20,10)
    gamma ~ dnorm(10,10)
    
    p ~ dnorm(0.9, 10)I(0,1)
    
    
# model for the initial population size: poisson priors

for (s in 1:n.sites){  

    # year 1
  N[s,1] ~ dpois(lambda) 

  for(t in 1:(n.years-1)){
    
    S[s,t+1] ~ dbin(phi[1], N[s,t])
    C[s,t+1] ~ dbin(phi[2], N[s,t]) 
    UC[s,t+1] ~ dbin(phi[3], N[s,t]) 
    D[s,t+1] <- N[s,t] - (S[s,t+1] + C[s,t+1] + UC[s, t+1]) 
    
    R[s,t+1] ~ dpois(gamma)

    N[s,t+1] <- (S[s,t+1] + R[s,t+1])

    } # t
    

     for (t in 1:n.years) {
      for (j in 1:nreps) {

    y[s,t,j] ~ dbin(p, N[s,t])

      } # j
     } # t
    } # s
}  
    ")

