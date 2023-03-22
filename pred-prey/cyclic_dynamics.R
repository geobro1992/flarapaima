##################
# Simple SES model 

# classic logistic growth model 
# plus fishing mortality
# similar to predator-prey models, there is negative feedback between fishing pressure and fish abundance

###################

# number of years
n.years = 200

# scenarios
F.sens = data.frame(scalar = rep(c(0, 0.5, 1, 2, 3, 4),2), 
                    lag = c(1,1,1,1,1,1,2,2,2,2,2,2)) # sensitivity to changes in fish abundance

# Define ecological state variables
N <- matrix(nrow = n.years, ncol = length(F.sens[,1])) # Population size
r <- matrix(nrow = n.years, ncol = length(F.sens[,1])) # Intrinsic growth rate

# Define social state variables
F <- matrix(nrow = n.years, ncol = length(F.sens[,1])) # Fishing effort
C <- matrix(nrow = n.years, ncol = length(F.sens[,1])) # Catch rate
P <- matrix(nrow = n.years, ncol = length(F.sens[,1])) # Policy regulations

# Define model parameters
alpha <- 0.05 # Growth rate coefficient
K <- 500 # Carrying capacity
gamma <- 0.02 # Catch rate coefficient
delta <- 0.1 # Policy coefficient

# Define initial values for year 1
N[1:2,] <- 300
F[1:2,] <- 0.2

# Define model equations
for(j in 1:length(F.sens[,1])){
  for (t in 3:n.years) {

        # Social dynamics
    C[t-1,j] <- gamma * F[t-1,j] * N[t-1,j]             # catch is a function of population size, fishing intensity, and catch rate (gamma)
    P[t-1,j] <- ((N[t-F.sens[j,2],j]/K) -  0.5) * F.sens[j,1]  # policy enforcement becomes more lenient when population approaches carrying capacity (can include lagged policy changes) 
    F[t,j] <- F[t-1,j] + (P[t-1,j] * F[t-1,j])          # fishing pressure changes with yearly changes in enforcement
    
    # Ecological dynamics
    r[t-1,j] <- alpha * N[t-1,j]
    N[t,j] <- N[t-1,j] + r[t-1,j] * (1 - N[t-1,j]/K) - C[t-1,j] 
    
    
  }
}


# Plot population size and fishing effort over time
pdf(file = "cyclic_lag.pdf", width = 6, height = 8)

par(mar=c(0,0,0,0))

layout(matrix(c(1,2,3,4,5,1,6,7,8,9),ncol=2),heights=c(1,3,3,3,3))

plot.new()

text(0.25,0.2,"Lag = 1 year",cex=2,font=2)
text(0.8,0.2,"Lag = 2 years",cex=2,font=2)

par(mar=c(2.5,2.5,1,1))
plot(N[,1], type = "l", ylab = "Population size", ylim = c(0,max(N, na.rm=T)+50))
plot(N[,3], type = "l", ylab = "Population size", ylim = c(0,max(N, na.rm=T)+50))
plot(N[,5], type = "l", ylab = "Population size", ylim = c(0,max(N, na.rm=T)+50))
plot(N[,6], type = "l", ylab = "Population size", ylim = c(0,max(N, na.rm=T)+50))
plot(N[,7], type = "l", ylab = "Population size", ylim = c(0,max(N, na.rm=T)+50))
plot(N[,9], type = "l", ylab = "Population size", ylim = c(0,max(N, na.rm=T)+50))
plot(N[,11], type = "l", ylab = "Population size", ylim = c(0,max(N, na.rm=T)+50))
plot(N[,12], type = "l", ylab = "Population size", ylim = c(0,max(N, na.rm=T)+50))

dev.off()

# Plot population size and fishing effort over time
pdf(file = "cyclic_lag_catch.pdf", width = 6, height = 8)

par(mar=c(0,0,0,0))

layout(matrix(c(1,2,3,4,5,1,6,7,8,9),ncol=2),heights=c(1,3,3,3,3))

plot.new()

text(0.25,0.2,"Lag = 1 year",cex=2,font=2)
text(0.8,0.2,"Lag = 2 years",cex=2,font=2)

par(mar=c(2.5,2.5,1,1))
plot(C[,1], type = "l", ylab = "Catch")
plot(C[,3], type = "l", ylab = "Catch")
plot(C[,5], type = "l", ylab = "Catch")
plot(C[,6], type = "l", ylab = "Catch")
plot(C[,7], type = "l", ylab = "Catch")
plot(C[,9], type = "l", ylab = "Catch")
plot(C[,11], type = "l", ylab = "Catch")
plot(C[,12], type = "l", ylab = "Catch")

dev.off()
