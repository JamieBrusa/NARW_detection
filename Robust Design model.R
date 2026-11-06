#Sample robust design data prep and model (set up for Cape Cod and Massachusetts Bays)
library(tidyverse)
library(R2jags)
library(coda)
library(ggmcmc)
library(MCMCvis)
library(here)
library(zoo)
library(doParallel)

#Read in data
ccb.ch <- as.matrix(read.csv("CCBch.csv"))
N.det <- readRDS("Ndet.RDS")
CountedSuper <- readRDS("CountedSuper.RDS")
SurveyYears <- read.csv("SurveyYears.csv")
aux.data <- read.csv("AuxData.csv")


###############################################################
##Bayesian robust design code adapted from Riecke et al. 2018##
###############################################################

##Data prep##

n.primary <- nrow(SurveyYears)
n.secondary <- as.numeric(SurveyYears$n.secondary)
EGNo <- aux.data$EGNo
Year <- SurveyYears$Year
first <- aux.data$first
last <- aux.data$last
sex <- aux.data$sex
n.ind <- nrow(aux.data)

n.states <- length(unique(as.vector(ccb.ch))) + 1 #Add 1 for the dead state (ccb.ch only has obs of live animals)

first.state <- vector()
for(i in 1:nrow(ccb.ch)){
  first.state[i] <- ccb.ch[i, first[i]]
}

z.init <- matrix(NA, nrow(ccb.ch), ncol(ccb.ch))
for (i in 1:nrow(ccb.ch)){
  if(first[i] < ncol(z.init)){
    z.init[i,(first[i] + 1):ncol(z.init)] <- 1
  }
}

##################
##Run JAGS model##
##################

rd.mod.code.c <- function() {
  
  #Set priors
  
  phi ~ dunif(0, 1)     #Apparent survival, not a focal parameter for our question
  gamP ~ dunif(0, 1)    #gamma' (probability of staying out of the study site from t-1 to t)
  gammaP <- 1 - gamP    #Probability of staying out of the study site
  beta.sex ~ dnorm(0, 0.1)
  
  mean.gamPP ~ dunif(0, 1)   #gamma'' (probability of moving out of the study site from t-1 to t)
  mu.gamPP <- log(mean.gamPP/(1-mean.gamPP))
  
  mean.p.cap ~ dunif(0, 0.6)
  mu.p.cap <- log(mean.p.cap/(1-mean.p.cap)) #logit transformation
  
  #Random effects
  
  #Detection
  for(t in 1:(n.primary - 1)){
    #Parameter constraints
    eps.year[t] ~ dnorm(0, tau.eps.year)
  }
  eps.year[n.primary] <- -sum(eps.year[1:(n.primary - 1)])
  
  tau.eps.year ~ dgamma(1, 1)
  
  sigma.eps.year <- 1/sqrt(tau.eps.year)
 
  #Temporary emigration
  for(t in 1:(n.primary - 1)){
    #Parameter constraints
    eps.year.te[t] ~ dnorm(0, tau.eps.year.te)
  }
  eps.year.te[n.primary] <- -sum(eps.year.te[1:(n.primary - 1)])
  
  for(i in 1:(n.ind - 1)){
    eps.ind.te[i] ~ dnorm(0, tau.eps.ind.te)
  }
  eps.ind.te[n.ind] <- -sum(eps.ind.te[1:(n.ind - 1)])
  
  
  tau.eps.year.te ~ dgamma(1, 1)
  tau.eps.ind.te ~ dgamma(1, 1)
  
  sigma.eps.year.te <- 1/sqrt(tau.eps.year.te)
  sigma.eps.ind.te <- 1/sqrt(tau.eps.ind.te)
  
  
  #Secondary occasions detection
  #assume constant detection for all secondaries within a primary
  #Need to loop when using random effect of year because logit cannot handle vectors
  
  for(i in 1:n.ind){
    for(t in (first[i]):n.primary){
      for (j in 1:max(n.secondary[1:n.primary])){
        logit(p.cap[t, j, i]) <- mu.p.cap + eps.year[Year[t]]
      }
    }
  }
  
  
  #Estimate number of detections from probability of detection per primary and secondary occasions and number of
  #individuals in the super-population
  for(i in 1:n.ind){
    for(t in (first[i]):n.primary){
      for(j in 1:n.secondary[t]){
        N.det[t, j, i] ~ dbin(p.cap[t, j, i], CountedSuper[t, j, i])
      }
    }
  }
  
  #Primary occasions detection or pooled detection probability (1 - prob of not detecting on any secondary occasions)
  for(i in 1:n.ind){
    for(t in (first[i]):n.primary){
      p.star[t, i] <- 1 - prod(1 - p.cap[t, 1:n.secondary[t], i])
    }
  }
  
  #Set up temporary emigration
  for(i in 1:n.ind){
    for(t in (first[i]):n.primary){
      logit(gamPP[t, i]) <- mu.gamPP + beta.sex*sex[i] + eps.ind.te[EGNo[i]] + eps.year.te[Year[t]]
    }
  }
  
  #State transition matrices
  #1 is in study area
  #2 out of study area
  #3 is dead
  #gamP and gamPP are opposite of MARK notation (i.e., gamPP is stay in the study site (not emigrate), and 
  #gamP is not stay out of study site (immigrate into study site)
  for(i in 1:n.ind){
    for(t in first[i]:n.primary){
      
      TPM[t, 1, 1, i] <- phi * gamPP[t, i]
      TPM[t, 1, 2, i] <- phi * (1 - gamPP[t, i])
      TPM[t, 1, 3, i] <- 1 - phi
      TPM[t, 2, 1, i] <- phi * gamP
      TPM[t, 2, 2, i] <- phi * (1 - gamP)
      TPM[t, 2, 3, i] <- 1 - phi
      TPM[t, 3, 1, i] <- 0
      TPM[t, 3, 2, i] <- 0
      TPM[t, 3, 3, i] <- 1
      
      
      #Observation matrices
      #1 down (first position) is in study area
      #2 down (first position) is out of study area
      #1 across (second position) is detected
      #2 across (second position) is not detected
      #3 is dead
      
      OPM[t, 1, 1, i] <- p.star[t, i]
      OPM[t, 1, 2, i] <- 1 - p.star[t, i]
      OPM[t, 2, 1, i] <- 0
      OPM[t, 2, 2, i] <- 1
      OPM[t, 3, 1, i] <- 0
      OPM[t, 3, 2, i] <- 1
    }
  }
  
  
  #Likelihood
  for(i in 1:n.ind){
    z[i, first[i]] <- ch[i, first[i]]
    for(t in (first[i] + 1):n.primary){
      z[i, t] ~ dcat(TPM[t, z[i, t-1], , i])   
      ch[i, t] ~ dcat(OPM[t, z[i, t], , i])
    } 
  }
}


rd.dat.c <- list(first = first, 
                 last = last,
                 ch = ccb.ch, 
                 sex = sex,
                 EGNo = EGNo,
                 Year = Year,
                 n.secondary = n.secondary, 
                 n.primary = n.primary, 
                 n.ind = n.ind,
                 N.det = N.det, 
                 CountedSuper = CountedSuper)

inits <- function(){list(z = z.init)}

params <- c("mu.p.cap", "phi", "gammaP", "beta.sex", "sigma.eps.year", "mu.gamPP", "sigma.eps.year.te",
            "sigma.eps.ind.te", "eps.year.te", "eps.ind.te")

cl <- makePSOCKcluster(3)
clusterSetRNGStream(cl)
registerDoParallel(cl)


rd.mod.c <- foreach(i = 1:3, .packages = "R2jags", .combine = 'c') %dopar% {
  load.module("lecuyer")
  jags(data = rd.dat.c,
       inits = inits,
       parameters.to.save = params,
       model.file = rd.mod.code.c,
       n.chains = 1,
       n.iter = 10000,
       n.burnin = 5000)
}



rd.mcmc <- as.mcmc(rd.mod.c)

#Investigate model diagnostics
test <- rd.mod.c[[2]]
sims.mat1 <- test[["sims.matrix"]]
test <- rd.mod.c[[8]]
sims.mat2 <- test[["sims.matrix"]]
test <- rd.mod.c[[14]]
sims.mat3 <- test[["sims.matrix"]]
sims.list <- list(sims.mat1, sims.mat2, sims.mat3)

mod_mcmc <- as.mcmc.list(lapply(sims.list, function(x) as.mcmc(x)))


#Diagnostics
mod_ggs <- ggs(mod_mcmc)

ggs_geweke(mod_ggs)
ggs_Rhat(mod_ggs)

ggs_traceplot(mod_ggs, c("phi"))
ggs_traceplot(mod_ggs, c("gamma"))
ggs_traceplot(mod_ggs, c("beta"))
ggs_traceplot(mod_ggs, c("mu"))
ggs_traceplot(mod_ggs, c("sigma"))


mod.summ <- MCMCsummary(mod_mcmc)









