##########
# OUTPUT
##########


################################################################
# model 2.1 - assumes fixed legal harvest and no illegal harvest
load("no_ih_2.RData")

traceplot(out1)


################
# posterior plots

# posteriors
phij = out1$sims.list$mu.phi.j # true value 0.4
sdj = out1$sims.list$sd.j # true value 0.05
phia = out1$sims.list$phi.a # true value 0.6
hest = 0.3 # true value 0.2

n.draws = length(phij)

par(mfrow = c(1,2))
hist(phij, breaks = 10, 
     main = "Juvenile Survival", xlim = c(0, 1))
abline(v = 0.4, col = "red", lwd = 2)

hist(phia, breaks = 10, 
     main = "Adult Survival", xlim = c(0, 1))
abline(v = 0.6, col = "red", lwd = 2)

n.sims = 10000
# storage vectors
NApred = matrix(NA, n.years, n.sims)
NJpred = matrix(NA, n.years, n.sims)
SApred = matrix(NA, n.years, n.sims)
SJpred = matrix(NA, n.years, n.sims)
Mpred = matrix(NA, n.years, n.sims)
Hpred = matrix(NA, n.years, n.sims)
Rpred = matrix(NA, n.years, n.sims)

NApredU = matrix(NA, n.years, n.sims)
NJpredU = matrix(NA, n.years, n.sims)
SApredU = matrix(NA, n.years, n.sims)
SJpredU = matrix(NA, n.years, n.sims)
MpredU = matrix(NA, n.years, n.sims)
HpredU = matrix(NA, n.years, n.sims)
RpredU = matrix(NA, n.years, n.sims)

NApredL = matrix(NA, n.years, n.sims)
NJpredL = matrix(NA, n.years, n.sims)
SApredL = matrix(NA, n.years, n.sims)
SJpredL = matrix(NA, n.years, n.sims)
MpredL = matrix(NA, n.years, n.sims)
HpredL = matrix(NA, n.years, n.sims)
RpredL = matrix(NA, n.years, n.sims)

# estimate population sizes
# year 1
NApred[1,] = rpois(n.sims, lambda0_a) 
NApredU[1,] = rpois(n.sims, lambda0_a) 
NApredL[1,] = rpois(n.sims, lambda0_a) 

NJpred[1,] = rpois(n.sims, lambda0_j)
NJpredU[1,] = rpois(n.sims, lambda0_j)
NJpredL[1,] = rpois(n.sims, lambda0_j)

# year 2
SApred[1,] = NApred[1] * quantile(phia, 0.5) 
SApredU[1,] = NApredU[1] * quantile(phia, 0.975) 
SApredL[1,] = NApredL[1] * quantile(phia, 0.025) 

SJpred[1,] = NJpred[1] * rnorm(n.sims, quantile(phij, 0.5), sample(sdj, n.sims, replace = T))
SJpredU[1,] = NJpredU[1] * rnorm(n.sims, quantile(phij, 0.975), sample(sdj, n.sims, replace = T))
SJpredL[1,] = NJpredL[1] * rnorm(n.sims, quantile(phij, 0.025), sample(sdj, n.sims, replace = T))

Mpred[1,] = SJpred[1,] * p.mat
MpredU[1,] = SJpredU[1,] * p.mat
MpredL[1,] = SJpredL[1,] * p.mat

NApred[2,] = (SApred[1,] + Mpred[1,])
NApredU[2,] = (SApredU[1,] + MpredU[1,])
NApredL[2,] = (SApredL[1,] + MpredL[1,])

NJpred[2,] = (SJpred[1,] - Mpred[1,])
NJpredU[2,] = (SJpredU[1,] - MpredU[1,])
NJpredL[2,] = (SJpredL[1,] - MpredL[1,])

for(t in 3:n.years){
  
  SApred[t-1,] = NApred[t-1,] * quantile(phia, 0.5) 
  SApredU[t-1,] = NApredU[t-1,] * quantile(phia, 0.975) 
  SApredL[t-1,] = NApredL[t-1,] * quantile(phia, 0.025) 
  
  Hpred[t-1,] = NApred[t-1,] * quantile(hest, 0.5)
  HpredU[t-1,] = NApredU[t-1,] * quantile(hest, 0.975)
  HpredL[t-1,] = NApredL[t-1,] * quantile(hest, 0.025)
  
  SJpred[t-1,] = NJpred[t-1,] * rnorm(n.sims, quantile(phij, 0.5), sample(sdj, n.sims, replace = T))
  SJpredU[t-1,] = NJpredU[t-1,] * rnorm(n.sims, quantile(phij, 0.975), sample(sdj, n.sims, replace = T))
  SJpredL[t-1,] = NJpredL[t-1,] * rnorm(n.sims, quantile(phij, 0.025), sample(sdj, n.sims, replace = T))
  
  Rpred[t-2,] = (alpha*NApred[t-2,])/(beta+NApred[t-2,])
  RpredU[t-2,] = (alpha*NApredU[t-2,])/(beta+NApredU[t-2,])
  RpredL[t-2,] = (alpha*NApredL[t-2,])/(beta+NApredL[t-2,])
  
  Mpred[t-1,] = SJpred[t-1,] * p.mat
  MpredU[t-1,] = SJpredU[t-1,] * p.mat
  MpredL[t-1,] = SJpredL[t-1,] * p.mat
  
  NApred[t,] = (SApred[t-1,] + Mpred[t-1,] - Hpred[t-1,])
  NApredU[t,] = (SApredU[t-1,] + MpredU[t-1,] - HpredU[t-1,])
  NApredL[t,] = (SApredL[t-1,] + MpredL[t-1,] - HpredL[t-1,])
  
  NJpred[t,] = (SJpred[t-1,] + Rpred[t-2,] - Mpred[t-1,])
  NJpredU[t,] = (SJpredU[t-1,] + RpredU[t-2,] - MpredU[t-1,])
  NJpredL[t,] = (SJpredL[t-1,] + RpredL[t-2,] - MpredL[t-1,])
  
  NApred[which(NApred < 1 )] <- 1 
  NApredU[which(NApredU < 1 )] <- 1 
  NApredL[which(NApredL < 1 )] <- 1 
  
  NJpred[which(NApred < 1 )] <- 1 
  NJpredU[which(NApredU < 1 )] <- 1 
  NJpredL[which(NApredL < 1 )] <- 1 
  
}


# plot real vs expected juveniles
par(mfrow = c(2,1))
plot(apply(N.J, 2, quantile, probs = 0.5, na.rm = TRUE), type = "l", ylim = c(0,650), axes = F, main = "Juveniles", xlab = "Year", ylab = "", col = "red")
lines(apply(N.J, 2, quantile, probs = 0.025, na.rm = TRUE), lty = "dashed", col = "red")
lines(apply(N.J, 2, quantile, probs = 0.975, na.rm = TRUE), lty = "dashed", col = "red")

lines(apply(NJpred, 1, quantile, probs = 0.5, na.rm = TRUE))
lines(apply(NJpred, 1, quantile, probs = 0.025, na.rm = TRUE), lty = "dashed")
lines(apply(NJpred, 1, quantile, probs = 0.975, na.rm = TRUE), lty = "dashed")

axis(1)
axis(2)

legend(20, 600, legend=c("True", "Predicted"),
       col=c("red", "black"), lty=1, cex=0.8)

# plot real vs expected adults
plot(apply(N.A, 2, quantile, probs = 0.5, na.rm = TRUE), type = "l", ylim = c(0,120), axes = F, main = "Adults", xlab = "", ylab = "", col = "red")
lines(apply(N.A, 2, quantile, probs = 0.025, na.rm = TRUE), lty = "dashed", col = "red")
lines(apply(N.A, 2, quantile, probs = 0.975, na.rm = TRUE), lty = "dashed", col = "red")

lines(apply(NApred, 1, quantile, probs = 0.5, na.rm = TRUE))
lines(apply(NApred, 1, quantile, probs = 0.025, na.rm = TRUE), lty = "dashed")
lines(apply(NApred, 1, quantile, probs = 0.975, na.rm = TRUE), lty = "dashed")

axis(1)
axis(2)


ex.p = vector()
est.ex = matrix(NA, n.years,3)

for(i in 1:n.years){
  
  # real extinction probability
  ex.p[i] = length(which(N.A[,i] < 5)) / n.sites
  
  # predicted extinction probability
  est.ex[i,1] = length(which(round(NApred[i,]) < 5)) / n.sims
  est.ex[i,2] = length(which(round(NApredU[i,]) < 5)) / n.sims
  est.ex[i,3] = length(which(round(NApredL[i,]) < 5)) / n.sims
  
}

par(mfrow = c(1,1), cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5, 
    font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)

plot(ex.p[5:n.years], type = "l", ylim = c(0,1), axes = F, xlab = "", ylab = "", col = "red", lwd = 2)
axis(1)
axis(2)
mtext("Year", side = 1, line = 3, cex = 2)
par(las = 0)
mtext("Proportion of Sites Extirpated", side = 2, line = 3.5, cex = 2)

legend(1, 1, legend=c("True", "Predicted"),
       col=c("red", "black"), lty=1, cex=1.5, lwd = 2)

lines(est.ex[5:n.years,1], lwd = 2)
lines(est.ex[5:n.years,2], lwd = 2, lty = "dashed")
lines(est.ex[5:n.years,3], lwd = 2, lty = "dashed")



################################################################
# model 2.2 - assumes fixed legal harvest and uninformed illegal harvest

load("ih_uninformed_2.RData")

traceplot(out2)


################
# posterior plots

# posteriors
phij = out2$sims.list$mu.phi.j # true value 0.4
sdj = out2$sims.list$sd.j # true value 0.05
phia = out2$sims.list$phi.a # true value 0.6
hest = 0.2 # true value 0.2
ihest = out2$sims.list$ih # true value = 0.2

n.draws = length(phij)

par(mfrow = c(2,2))
hist(phij, breaks = 10, 
     main = "Juvenile Survival", xlim = c(0, 1))
abline(v = 0.4, col = "red", lwd = 2)

hist(phia, breaks = 10, 
     main = "Adult Survival", xlim = c(0, 1))
abline(v = 0.6, col = "red", lwd = 2)

hist(ihest, breaks = 10, 
     main = "Illegal Harvest", xlim = c(0, 1))
abline(v = 0.2, col = "red", lwd = 2)

n.sims = 10000
# storage vectors
NApred = matrix(NA, n.years, n.sims)
NJpred = matrix(NA, n.years, n.sims)
SApred = matrix(NA, n.years, n.sims)
SJpred = matrix(NA, n.years, n.sims)
Mpred = matrix(NA, n.years, n.sims)
Hpred = matrix(NA, n.years, n.sims)
Rpred = matrix(NA, n.years, n.sims)

NApredU = matrix(NA, n.years, n.sims)
NJpredU = matrix(NA, n.years, n.sims)
SApredU = matrix(NA, n.years, n.sims)
SJpredU = matrix(NA, n.years, n.sims)
MpredU = matrix(NA, n.years, n.sims)
HpredU = matrix(NA, n.years, n.sims)
RpredU = matrix(NA, n.years, n.sims)

NApredL = matrix(NA, n.years, n.sims)
NJpredL = matrix(NA, n.years, n.sims)
SApredL = matrix(NA, n.years, n.sims)
SJpredL = matrix(NA, n.years, n.sims)
MpredL = matrix(NA, n.years, n.sims)
HpredL = matrix(NA, n.years, n.sims)
RpredL = matrix(NA, n.years, n.sims)

# estimate population sizes
# year 1
NApred[1,] = rpois(n.sims, lambda0_a) 
NApredU[1,] = rpois(n.sims, lambda0_a) 
NApredL[1,] = rpois(n.sims, lambda0_a) 

NJpred[1,] = rpois(n.sims, lambda0_j)
NJpredU[1,] = rpois(n.sims, lambda0_j)
NJpredL[1,] = rpois(n.sims, lambda0_j)

# year 2
SApred[1,] = NApred[1] * quantile(phia, 0.5) 
SApredU[1,] = NApredU[1] * quantile(phia, 0.975) 
SApredL[1,] = NApredL[1] * quantile(phia, 0.025) 

SJpred[1,] = NJpred[1] * rnorm(n.sites, quantile(phij, 0.5), sample(sdj, n.sims, replace = T))
SJpredU[1,] = NJpredU[1] * rnorm(n.sites, quantile(phij, 0.975), sample(sdj, n.sims, replace = T))
SJpredL[1,] = NJpredL[1] * rnorm(n.sites, quantile(phij, 0.025), sample(sdj, n.sims, replace = T))

Mpred[1,] = SJpred[1,] * p.mat
MpredU[1,] = SJpredU[1,] * p.mat
MpredL[1,] = SJpredL[1,] * p.mat

NApred[2,] = (SApred[1,] + Mpred[1,])
NApredU[2,] = (SApredU[1,] + MpredU[1,])
NApredL[2,] = (SApredL[1,] + MpredL[1,])

NJpred[2,] = (SJpred[1,] - Mpred[1,])
NJpredU[2,] = (SJpredU[1,] - MpredU[1,])
NJpredL[2,] = (SJpredL[1,] - MpredL[1,])

for(t in 3:n.years){
  
  SApred[t-1,] = NApred[t-1,] * quantile(phia, 0.5) 
  SApredU[t-1,] = NApredU[t-1,] * quantile(phia, 0.975) 
  SApredL[t-1,] = NApredL[t-1,] * quantile(phia, 0.025) 
  
  Hpred[t-1,] = NApred[t-1,] * (quantile(hest, 0.5) + quantile(ihest, 0.5))
  HpredU[t-1,] = NApredU[t-1,] * (quantile(hest, 0.975) + quantile(ihest, 0.975))
  HpredL[t-1,] = NApredL[t-1,] * (quantile(hest, 0.025) + quantile(ihest, 0.025))
  
  SJpred[t-1,] = NJpred[t-1,] * rnorm(n.sims, quantile(phij, 0.5), sample(sdj, n.sims, replace = T))
  SJpredU[t-1,] = NJpredU[t-1,] * rnorm(n.sims, quantile(phij, 0.975), sample(sdj, n.sims, replace = T))
  SJpredL[t-1,] = NJpredL[t-1,] * rnorm(n.sims, quantile(phij, 0.025), sample(sdj, n.sims, replace = T))
  
  Rpred[t-2,] = (alpha*NApred[t-2,])/(beta+NApred[t-2,])
  RpredU[t-2,] = (alpha*NApredU[t-2,])/(beta+NApredU[t-2,])
  RpredL[t-2,] = (alpha*NApredL[t-2,])/(beta+NApredL[t-2,])
  
  Mpred[t-1,] = SJpred[t-1,] * p.mat
  MpredU[t-1,] = SJpredU[t-1,] * p.mat
  MpredL[t-1,] = SJpredL[t-1,] * p.mat
  
  NApred[t,] = (SApred[t-1,] + Mpred[t-1,] - Hpred[t-1,])
  NApredU[t,] = (SApredU[t-1,] + MpredU[t-1,] - HpredU[t-1,])
  NApredL[t,] = (SApredL[t-1,] + MpredL[t-1,] - HpredL[t-1,])
  
  NJpred[t,] = (SJpred[t-1,] + Rpred[t-2,] - Mpred[t-1,])
  NJpredU[t,] = (SJpredU[t-1,] + RpredU[t-2,] - MpredU[t-1,])
  NJpredL[t,] = (SJpredL[t-1,] + RpredL[t-2,] - MpredL[t-1,])
  
  NApred[which(NApred < 1 )] <- 1 
  NApredU[which(NApredU < 1 )] <- 1 
  NApredL[which(NApredL < 1 )] <- 1 
  
  NJpred[which(NApred < 1 )] <- 1 
  NJpredU[which(NApredU < 1 )] <- 1 
  NJpredL[which(NApredL < 1 )] <- 1 
  
}


# plot real vs expected juveniles
par(mfrow = c(2,1))
plot(apply(N.J, 2, quantile, probs = 0.5, na.rm = TRUE), type = "l", ylim = c(0,650), axes = F, main = "Juveniles", xlab = "Year", ylab = "", col = "red")
lines(apply(N.J, 2, quantile, probs = 0.025, na.rm = TRUE), lty = "dashed", col = "red")
lines(apply(N.J, 2, quantile, probs = 0.975, na.rm = TRUE), lty = "dashed", col = "red")

lines(apply(NJpred, 1, quantile, probs = 0.5, na.rm = TRUE))
lines(apply(NJpred, 1, quantile, probs = 0.025, na.rm = TRUE), lty = "dashed")
lines(apply(NJpred, 1, quantile, probs = 0.975, na.rm = TRUE), lty = "dashed")

axis(1)
axis(2)

legend(20, 600, legend=c("True", "Predicted"),
       col=c("red", "black"), lty=1, cex=0.8)

# plot real vs expected adults
plot(apply(N.A, 2, quantile, probs = 0.5, na.rm = TRUE), type = "l", ylim = c(0,120), axes = F, main = "Adults", xlab = "", ylab = "", col = "red")
lines(apply(N.A, 2, quantile, probs = 0.025, na.rm = TRUE), lty = "dashed", col = "red")
lines(apply(N.A, 2, quantile, probs = 0.975, na.rm = TRUE), lty = "dashed", col = "red")

lines(apply(NApred, 1, quantile, probs = 0.5, na.rm = TRUE))
lines(apply(NApred, 1, quantile, probs = 0.025, na.rm = TRUE), lty = "dashed")
lines(apply(NApred, 1, quantile, probs = 0.975, na.rm = TRUE), lty = "dashed")

axis(1)
axis(2)


ex.p = vector()
est.ex = matrix(NA, n.years,3)

for(i in 1:n.years){
  
  # real extinction probability
  ex.p[i] = length(which(N.A[,i] < 2)) / n.sites
  
  # predicted extinction probability
  est.ex[i,1] = length(which(NApred[i,] < 2)) / n.sims
  est.ex[i,2] = length(which(NApredU[i,] < 2)) / n.sims
  est.ex[i,3] = length(which(NApredL[i,] < 2)) / n.sims
  
}

par(mfrow = c(1,1), cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5, 
    font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)

plot(ex.p, type = "l", ylim = c(0,1), axes = F, xlab = "", ylab = "", col = "red", lwd = 2)
axis(1)
axis(2)
mtext("Year", side = 1, line = 3, cex = 2)
par(las = 0)
mtext("Proportion of Sites Extirpated", side = 2, line = 3.5, cex = 2)

legend(1, 1, legend=c("True", "Predicted"),
       col=c("red", "black"), lty=1, cex=1.5, lwd = 2)

lines(est.ex[,1], lwd = 2)
lines(est.ex[,2], lwd = 2, lty = "dashed")
lines(est.ex[,3], lwd = 2, lty = "dashed")



###########################################################################
# model 2.3 - WITH survey data that informs rates of unreported fishing

load("ih_informed_2.RData")

traceplot(out3)

################
# posterior plots

# posteriors
phij = out3$sims.list$mu.phi.j # true value 0.4
sdj = out3$sims.list$sd.j # true value 0.05
phia = out3$sims.list$phi.a # true value 0.6
hest = 0.2 # true value 0.2
ihest = out3$sims.list$ih # true value = 0.2

n.draws = length(phij)

par(mfrow = c(2,2))
hist(phij, breaks = 10, 
     main = "Juvenile Survival", xlim = c(0, 1))
abline(v = 0.4, col = "red", lwd = 2)

hist(phia, breaks = 10, 
     main = "Adult Survival", xlim = c(0, 1))
abline(v = 0.6, col = "red", lwd = 2)

hist(ihest, breaks = 10, 
     main = "Illegal Harvest", xlim = c(0, 1))
abline(v = 0.2, col = "red", lwd = 2)

n.sims = 10000
# storage vectors
NApred = matrix(NA, n.years, n.sims)
NJpred = matrix(NA, n.years, n.sims)
SApred = matrix(NA, n.years, n.sims)
SJpred = matrix(NA, n.years, n.sims)
Mpred = matrix(NA, n.years, n.sims)
Hpred = matrix(NA, n.years, n.sims)
Rpred = matrix(NA, n.years, n.sims)

NApredU = matrix(NA, n.years, n.sims)
NJpredU = matrix(NA, n.years, n.sims)
SApredU = matrix(NA, n.years, n.sims)
SJpredU = matrix(NA, n.years, n.sims)
MpredU = matrix(NA, n.years, n.sims)
HpredU = matrix(NA, n.years, n.sims)
RpredU = matrix(NA, n.years, n.sims)

NApredL = matrix(NA, n.years, n.sims)
NJpredL = matrix(NA, n.years, n.sims)
SApredL = matrix(NA, n.years, n.sims)
SJpredL = matrix(NA, n.years, n.sims)
MpredL = matrix(NA, n.years, n.sims)
HpredL = matrix(NA, n.years, n.sims)
RpredL = matrix(NA, n.years, n.sims)

# estimate population sizes
# year 1
NApred[1,] = rpois(n.sims, lambda0_a) 
NApredU[1,] = rpois(n.sims, lambda0_a) 
NApredL[1,] = rpois(n.sims, lambda0_a) 

NJpred[1,] = rpois(n.sims, lambda0_j)
NJpredU[1,] = rpois(n.sims, lambda0_j)
NJpredL[1,] = rpois(n.sims, lambda0_j)

# year 2
SApred[1,] = NApred[1] * quantile(phia, 0.5) 
SApredU[1,] = NApredU[1] * quantile(phia, 0.975) 
SApredL[1,] = NApredL[1] * quantile(phia, 0.025) 

SJpred[1,] = NJpred[1] * rnorm(n.sites, quantile(phij, 0.5), sample(sdj, n.sims, replace = T))
SJpredU[1,] = NJpredU[1] * rnorm(n.sites, quantile(phij, 0.975), sample(sdj, n.sims, replace = T))
SJpredL[1,] = NJpredL[1] * rnorm(n.sites, quantile(phij, 0.025), sample(sdj, n.sims, replace = T))

Mpred[1,] = SJpred[1,] * p.mat
MpredU[1,] = SJpredU[1,] * p.mat
MpredL[1,] = SJpredL[1,] * p.mat

NApred[2,] = (SApred[1,] + Mpred[1,])
NApredU[2,] = (SApredU[1,] + MpredU[1,])
NApredL[2,] = (SApredL[1,] + MpredL[1,])

NJpred[2,] = (SJpred[1,] - Mpred[1,])
NJpredU[2,] = (SJpredU[1,] - MpredU[1,])
NJpredL[2,] = (SJpredL[1,] - MpredL[1,])

for(t in 3:n.years){
  
  SApred[t-1,] = NApred[t-1,] * quantile(phia, 0.5) 
  SApredU[t-1,] = NApredU[t-1,] * quantile(phia, 0.975) 
  SApredL[t-1,] = NApredL[t-1,] * quantile(phia, 0.025) 
  
  Hpred[t-1,] = NApred[t-1,] * (quantile(hest, 0.5) + quantile(ihest, 0.5))
  HpredU[t-1,] = NApredU[t-1,] * (quantile(hest, 0.975) + quantile(ihest, 0.975))
  HpredL[t-1,] = NApredL[t-1,] * (quantile(hest, 0.025) + quantile(ihest, 0.025))
  
  SJpred[t-1,] = NJpred[t-1,] * rnorm(n.sims, quantile(phij, 0.5), sample(sdj, n.sims, replace = T))
  SJpredU[t-1,] = NJpredU[t-1,] * rnorm(n.sims, quantile(phij, 0.975), sample(sdj, n.sims, replace = T))
  SJpredL[t-1,] = NJpredL[t-1,] * rnorm(n.sims, quantile(phij, 0.025), sample(sdj, n.sims, replace = T))
  
  Rpred[t-2,] = (alpha*NApred[t-2,])/(beta+NApred[t-2,])
  RpredU[t-2,] = (alpha*NApredU[t-2,])/(beta+NApredU[t-2,])
  RpredL[t-2,] = (alpha*NApredL[t-2,])/(beta+NApredL[t-2,])
  
  Mpred[t-1,] = SJpred[t-1,] * p.mat
  MpredU[t-1,] = SJpredU[t-1,] * p.mat
  MpredL[t-1,] = SJpredL[t-1,] * p.mat
  
  NApred[t,] = (SApred[t-1,] + Mpred[t-1,] - Hpred[t-1,] - IHpred[t,])
  NApredU[t,] = (SApredU[t-1,] + MpredU[t-1,] - HpredU[t-1,] - IHpredL[t,])
  NApredL[t,] = (SApredL[t-1,] + MpredL[t-1,] - HpredL[t-1,] - IHpredU[t,])
  
  NJpred[t,] = (SJpred[t-1,] + Rpred[t-2,] - Mpred[t-1,])
  NJpredU[t,] = (SJpredU[t-1,] + RpredU[t-2,] - MpredU[t-1,])
  NJpredL[t,] = (SJpredL[t-1,] + RpredL[t-2,] - MpredL[t-1,])
  
  NApred[which(NApred < 1 )] <- 1 
  NApredU[which(NApredU < 1 )] <- 1 
  NApredL[which(NApredL < 1 )] <- 1 
  
  NJpred[which(NApred < 1 )] <- 1 
  NJpredU[which(NApredU < 1 )] <- 1 
  NJpredL[which(NApredL < 1 )] <- 1 
  
}


# plot real vs expected juveniles
par(mfrow = c(2,1))
plot(apply(N.J, 2, quantile, probs = 0.5, na.rm = TRUE), type = "l", ylim = c(0,650), axes = F, main = "Juveniles", xlab = "Year", ylab = "", col = "red")
lines(apply(N.J, 2, quantile, probs = 0.025, na.rm = TRUE), lty = "dashed", col = "red")
lines(apply(N.J, 2, quantile, probs = 0.975, na.rm = TRUE), lty = "dashed", col = "red")

lines(apply(NJpred, 1, quantile, probs = 0.5, na.rm = TRUE))
lines(apply(NJpred, 1, quantile, probs = 0.025, na.rm = TRUE), lty = "dashed")
lines(apply(NJpred, 1, quantile, probs = 0.975, na.rm = TRUE), lty = "dashed")

axis(1)
axis(2)

legend(20, 600, legend=c("True", "Predicted"),
       col=c("red", "black"), lty=1, cex=0.8)

# plot real vs expected adults
plot(apply(N.A, 2, quantile, probs = 0.5, na.rm = TRUE), type = "l", ylim = c(0,120), axes = F, main = "Adults", xlab = "", ylab = "", col = "red")
lines(apply(N.A, 2, quantile, probs = 0.025, na.rm = TRUE), lty = "dashed", col = "red")
lines(apply(N.A, 2, quantile, probs = 0.975, na.rm = TRUE), lty = "dashed", col = "red")

lines(apply(NApred, 1, quantile, probs = 0.5, na.rm = TRUE))
lines(apply(NApred, 1, quantile, probs = 0.025, na.rm = TRUE), lty = "dashed")
lines(apply(NApred, 1, quantile, probs = 0.975, na.rm = TRUE), lty = "dashed")

axis(1)
axis(2)


ex.p = vector()
est.ex = matrix(NA, n.years,3)

for(i in 1:n.years){
  
  # real extinction probability
  ex.p[i] = length(which(N.A[,i] < 2)) / n.sites
  
  # predicted extinction probability
  est.ex[i,1] = length(which(NApred[i,] < 2)) / n.sims
  est.ex[i,2] = length(which(NApredU[i,] < 2)) / n.sims
  est.ex[i,3] = length(which(NApredL[i,] < 2)) / n.sims
  
}

par(mfrow = c(1,1), cex.main = 1.5, mar = c(5, 6, 4, 5) + 0.1, mgp = c(3.5, 1, 0), cex.lab = 1.5, 
    font.lab = 2, cex.axis = 1.3, bty = "n", las = 1)

plot(ex.p, type = "l", ylim = c(0,1), axes = F, xlab = "", ylab = "", col = "red", lwd = 2)
axis(1)
axis(2)
mtext("Year", side = 1, line = 3, cex = 2)
par(las = 0)
mtext("Proportion of Sites Extirpated", side = 2, line = 3.5, cex = 2)

legend(1, 1, legend=c("True", "Predicted"),
       col=c("red", "black"), lty=1, cex=1.5, lwd = 2)

lines(est.ex[,1], lwd = 2)
lines(est.ex[,2], lwd = 2, lty = "dashed")
lines(est.ex[,3], lwd = 2, lty = "dashed")


