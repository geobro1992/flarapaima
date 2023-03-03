##################
# Simple SES model 

# classic logistic growth model 
# plus fishing mortality
# similar to predator-prey models, there is negative feedback between fishing pressure and fish abundance

###################

# Define ecological state variables
N <- numeric(200) # Population size
r <- numeric(200) # Intrinsic growth rate

# Define social state variables
F <- numeric(200) # Fishing effort
C <- numeric(200) # Catch rate
P <- numeric(200) # Policy regulations

# Define model parameters
alpha <- 0.05 # Growth rate coefficient
K <- 500 # Carrying capacity
gamma <- 0.02 # Catch rate coefficient
delta <- 0.1 # Policy coefficient

# Define initial values for year 1
N[1] <- 300
F[1] <- 0.2

# Define model equations
for (t in 2:200) {
  # Ecological dynamics
  r[t] <- alpha * N[t-1]
  N[t] <- N[t-1] + r[t] * (1 - N[t-1]/K) - C[t-1] 
  
  # Social dynamics
  C[t] <- gamma * F[t-1] * N[t-1]     # catch is a function of population size, fishing intensity, and catch rate (gamma)
  P[t] <- (N[t-1]/K) - 0.5            # policy enforcement becomes more lenient when population approaches carrying capacity 
  F[t] <- F[t-1] + (P[t-1] * F[t-1])  # fishing pressure changes with yearly changes in enforcement
}

# Plot population size and fishing effort over time
par(mfrow = c(2,1), mar = c(2,5,2,2))
plot(N, type = "l", ylab = "Population size")
plot(C, type = "l", ylab = "Catch / Fishing Effort", col = "red")
