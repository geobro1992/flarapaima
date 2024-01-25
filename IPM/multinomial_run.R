library(jagsUI)

######################
# simulated count data
######################

n.sites = 100 # scalar, number of habitat patches
n.years = 20 # scalar, number of years
nreps = 3 # scalar, number of surveys within a primary sampling period (needed to estimate observation error)


#################################
# compare with unreported harvest (proportional to abundance)

# storage arrays
N = array(dim = c(n.sites, n.years)) # an array[snite, nyear, nreps] containing the true counts of adults
S = array(dim = c(n.sites, n.years)) # an array[snite, nyear, nreps] containing the surviving adults from previous year
H = array(dim = c(n.sites, n.years)) # an array[snite, nyear, nreps] containing the number of fished animals
IH = array(dim = c(n.sites, n.years))
R = array(dim = c(n.sites, n.years)) # an array[snite, nyear, nreps] containing the number of new recruits

lambda = 20       # mean number of adults in year 1 

psi = 0.8           # adult survival
F = 0.6               # legal fishing pressure
gamma = 10
r = 1
p = 0.9

ps = c(psi - (psi*F), psi*F*r, psi*F*(1-r), (1-psi))
ps = ps/sum(ps)

# create simulated data

# year 1
N[,1] = rpois(n.sites, lambda) # if pop sizes across sites are overdispersed, can't use poisson

for(t in 1:(n.years-1)){
  
  R[,t+1] = rpois(n.sites, gamma)
  
  for(s in 1:n.sites){
  
    m = rmultinom(1, N[s,t], prob = ps)
    
    S[s,t+1] = m[1] 
    H[s,t+1] = m[2]
    IH[s,t+1] = m[3]
    
    N[s,t+1] = (S[s,t+1] + R[s,t+1])
    
  }
}

# true declines
IHprop.TRUE =  N[,1] - N[,10]
IHprop.TRUE = ifelse(IHprop.TRUE < 0, 0, 1)
sum(IHprop.TRUE)


y = array(dim = c(n.sites, n.years, nreps))

# simulate observation error of counts (but to accurately estimate, you need repeated counts within a year)
for(j in 1:nreps){
  y[,,j] = rbinom(n.sites*n.years, N, p)     # observed counts of adults
}

plot(N, y[,,1]) # observed vs true counts of adults
abline(0,1)

# observed catch
C = round(H)


op <- par(mfrow = c(2, 2), mar = c(5, 5, 4, 3), cex.lab = 1.5, 
          cex.axis = 1.5)
on.exit(par(op))
matplot(t(N), type = "l", main = paste("Population trajectories"), 
        lty = 1, lwd = 3, las = 1, frame = FALSE, xlab = "Year", 
        ylab = "N")
matplot(t(S), type = "l", main = "Number of apparent survivors", 
        lty = 1, lwd = 3, las = 1, frame = FALSE, xlab = "Year", 
        ylab = "S")
hist(N[, 1], main = "Distribution of N in first year", 
     breaks = 50, col = "grey")
hist(N[, n.years], main = "Distribution of N in last year", 
     breaks = 50, col = "grey")

############
# Model Runs
############

#Bundle data
jags.data <- list(y = y,                                       # counts of adults and juveniles
                  C=C,                                        # catch data
                  n.years=ncol(y), n.sites = nrow(y), nreps = nreps) 


#Initial values
Rst <- apply(y, c(1,2), max, na.rm = TRUE)
Rst[Rst == '-Inf'] <- 1
Rst[,1] <- NA

Nst <- array(NA, dim = dim(Rst))
Nst[,1] <- N[,1]

inits <- function(){list(lambda = runif(1, 15, 25), F = runif(1, 0.5, 0.6), 
                         gamma = runif(1, 5, 15), p = runif(1, 0.9, 1), R = Rst+1, N = Nst+2)}


#parameters monitored
parameters <- c("phi", "p", "gamma", "lambda", "F", "psi", "phi_raw")

#MCMC settings
na <- 1000 ; ni <- 30000 ; nt <- 10 ; nb <- 20000 ; nc <- 3


#################
# baseline
###################

#call JAGS from R
(out= jags(data = jags.data, inits = inits, parameters.to.save = parameters, model.file = "ih_multi_baseline.txt", 
           n.iter=ni, n.burnin = nb, n.chains= nc,
           n.thin= nt, n.adapt = na, parallel = TRUE))

traceplot(out)
save(out, file = "ih_multi_baseline.RData")
load("ih_multi_baseline.RData")


#parameters monitored
parameters <- c("phi", "p", "gamma", "lambda", "F", "r", "psi")

#MCMC settings
na <- 1000 ; ni <- 30000 ; nt <- 10 ; nb <- 20000 ; nc <- 3


#################
# uninformed
###################

#call JAGS from R
(out= jags(data = jags.data, inits = inits, parameters.to.save = parameters, model.file = "ih_multi_uninformed.txt", 
            n.iter=ni, n.burnin = nb, n.chains= nc,
            n.thin= nt, n.adapt = na, parallel = TRUE))

traceplot(out)
save(out, file = "ih_multi_uninformed.RData")
load("ih_multi_uninformed.RData")

# posteriors
lest = out$sims.list$lambda # true value 20
phi = out$sims.list$phi # true values are 0.32, 0.24, 0.24, 0.20
pest = out$sims.list$p # true value 0.9
gest = out$sims.list$gamma # true value 0.3
fest = out$sims.list$F # true value 0.6
rest = out$sims.list$r # true value 0.5
psiest = out$sims.list$psi # true value 0.8


hist(lest, breaks = 10, 
     main = "Initial Pop Size", xlim = c(10, 30))
abline(v = 20, col = "red", lwd = 2)



par(mfrow = c(1,2))

hist(gest, breaks = 10, 
     main = "Recruitment", xlim = c(5, 15))
abline(v = 10, col = "red", lwd = 2)

hist(pest, breaks = 10, 
     main = "Detection", xlim = c(0.5, 1))
abline(v = 0.9, col = "red", lwd = 2)


png("ih_prop_r.png", width = 6, height = 6, units = "in", res = 600)

par(mfrow = c(2,2))

hist(psiest, breaks = 10, 
     main = "Survival", xlim = c(0, 1))
abline(v = 0.8, col = "red", lwd = 2)

hist(gest, breaks = 10, 
     main = "Recruitment", xlim = c(5, 15))
abline(v = 10, col = "red", lwd = 2)

hist(fest, breaks = 10, 
     main = "Total Fishing", xlim = c(0, 1))
abline(v = 0.6, col = "red", lwd = 2)

hist(rest, breaks = 10, 
     main = "Reporting Rate", xlim = c(0, 1))
abline(v = 0.8, col = "red", lwd = 2)

dev.off()

png("ih_prop_r2.png", width = 6, height = 4, units = "in", res = 600)

par(mfrow = c(1,2))

hist(fest * rest, breaks = 10, 
     main = "Reported Catch", xlim = c(0, 1))
abline(v = F*r, col = "red", lwd = 2)

hist(fest * (1-rest), breaks = 10, 
     main = "Unreported Catch", xlim = c(0, 1))
abline(v = F*(1-r), col = "red", lwd = 2)

dev.off()


png("ih_prop_r3.png", width = 8, height = 2.5, units = "in", res = 600)

par(mfrow = c(1,3))

hist(phi[,1], breaks = 10, 
     main = "Adult Survival", xlim = c(0, 0.6))
abline(v = ps[1], col = "red", lwd = 2)

hist(gest, breaks = 10, 
     main = "Recruitment", xlim = c(5, 15))
abline(v = 10, col = "red", lwd = 2)

hist(phi[,2], breaks = 10, 
     main = "Legal Fishing", xlim = c(0.4, 0.6))
abline(v = ps[2], col = "red", lwd = 2)

dev.off()

#################
# informed
###################

#call JAGS from R
(out2= jags(data = jags.data, inits = inits, parameters.to.save = parameters, model.file = "ih_multi_informed.txt", 
            n.iter=ni, n.burnin = nb, n.chains= nc,
            n.thin= nt, n.adapt = na, parallel = TRUE))

traceplot(out2)
save(out2, file = "ih_multi_informed.RData")
load("ih_multi_informed.RData")

# posteriors
lest = out2$sims.list$lambda # true value 20
phi = out2$sims.list$phi # true values are 0.32, 0.38, 0.1, 0.20
pest = out2$sims.list$p # true value 0.9
gest = out2$sims.list$gamma # true value 0.3
fest = out2$sims.list$F # true value 0.6
rest = out2$sims.list$r # true value 0.5
psiest = out2$sims.list$psi # true value 0.8


par(mfrow = c(1,2))

hist(lest, breaks = 10, 
     main = "Initial Pop Size", xlim = c(10, 30))
abline(v = lambda, col = "red", lwd = 2)

hist(pest, breaks = 10, 
     main = "Detection", xlim = c(0.5, 1))
abline(v = p, col = "red", lwd = 2)


png("ih_prop_rin.png", width = 6, height = 6, units = "in", res = 600)

par(mfrow = c(2,2))

hist(psiest, breaks = 10, 
     main = "Survival", xlim = c(0, 1))
abline(v = 0.8, col = "red", lwd = 2)

hist(gest, breaks = 10, 
     main = "Recruitment", xlim = c(5, 15))
abline(v = 10, col = "red", lwd = 2)

hist(fest, breaks = 10, 
     main = "Total Fishing", xlim = c(0, 1))
abline(v = 0.6, col = "red", lwd = 2)

hist(rest, breaks = 10, 
     main = "Reporting Rate", xlim = c(0, 1))
abline(v = 0.8, col = "red", lwd = 2)

dev.off()

png("ih_prop_rin2.png", width = 6, height = 4, units = "in", res = 600)

par(mfrow = c(1,2))

hist(fest * rest, breaks = 10, 
     main = "Reported Catch", xlim = c(0, 1))
abline(v = F*r, col = "red", lwd = 2)

hist(fest * (1-rest), breaks = 10, 
     main = "Unreported Catch", xlim = c(0, 1))
abline(v = F*(1-r), col = "red", lwd = 2)

dev.off()



load("ih_multi_informed.RData")
load("ih_multi_uninformed.RData")

# posteriors
lest = out$sims.list$lambda # true value 20
phi = out$sims.list$phi # true values are 0.32, 0.24, 0.24, 0.20
pest = out$sims.list$p # true value 0.9
gest = out$sims.list$gamma # true value 0.3
fest = out$sims.list$F # true value 0.6
rest = out$sims.list$r # true value 0.5
psiest = out$sims.list$psi # true value 0.8

# posteriors
lest2 = out2$sims.list$lambda # true value 20
phi2 = out2$sims.list$phi # true values are 0.32, 0.38, 0.1, 0.20
pest2 = out2$sims.list$p # true value 0.9
gest2 = out2$sims.list$gamma # true value 0.3
fest2 = out2$sims.list$F # true value 0.6
rest2 = out2$sims.list$r # true value 0.5
psiest2 = out2$sims.list$psi # true value 0.8

png("ih_prop_r_combined.png", width = 8, height = 6, units = "in", res = 600)

par(mfrow = c(2,2))
op <- par(cex.main = 1.5, mar = c(3, 6, 2, 2) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5 , font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)

hist(psiest, breaks = 10, col=rgb(1,0,0,0.5), 
     main = "Adult Survival", xlim = c(0, 1), ylim = c(0, 700), xlab = "")

hist(psiest2, breaks = 10, col=rgb(0,0,1,0.5),
     main = "", xlim = c(0, 1), ylim = c(0, 700), xlab = "", add = T)
abline(v = 0.8, col = "red", lwd = 2)


hist(gest, breaks = 10, col=rgb(1,0,0,0.5), 
     main = "Recruitment", xlim = c(8, 12), xlab = "", ylab = "")

hist(gest2, breaks = 10, col=rgb(0,0,1,0.5), 
     main = "", xlim = c(8, 12), xlab = "", ylab = "", add = T)
abline(v = 10, col = "red", lwd = 2)

# Add legend
legend("topright", legend=c("Uninformed","Informed"), col=c(rgb(1,0,0,0.5), rgb(0,0,1,0.5)), pt.cex=2, pch=15 )


hist(fest, breaks = 10, col=rgb(1,0,0,0.5), 
     main = "Total Fishing", xlim = c(0, 1), xlab = "")

hist(fest2, breaks = 10, col=rgb(0,0,1,0.5), 
     main = "", xlim = c(0, 1), xlab = "", ylab = "", add = T)
abline(v = 0.6, col = "red", lwd = 2)


hist(rest, breaks = 10, col=rgb(1,0,0,0.5), 
     main = "Reporting Rate", xlim = c(0, 1), xlab = "", ylab = "")

hist(rest2, breaks = 10, col=rgb(0,0,1,0.5), 
     main = "", xlim = c(0, 1), xlab = "", ylab = "", add = T)
abline(v = 0.8, col = "red", lwd = 2)

dev.off()



png("ih_prop_r2_combined.png", width = 6, height = 6, units = "in", res = 600)

par(mfrow = c(2,2))

hist(fest * rest, breaks = 10, 
     main = "Reported Catch", xlim = c(0, 1))
abline(v = F*r, col = "red", lwd = 2)

hist(fest * (1-rest), breaks = 10, 
     main = "Unreported Catch", xlim = c(0, 1))
abline(v = F*(1-r), col = "red", lwd = 2)


hist(fest2 * rest2, breaks = 10, 
     main = "Reported Catch", xlim = c(0, 1))
abline(v = F*r, col = "red", lwd = 2)

hist(fest2 * (1-rest2), breaks = 10, 
     main = "Unreported Catch", xlim = c(0, 1))
abline(v = F*(1-r), col = "red", lwd = 2)

dev.off()
