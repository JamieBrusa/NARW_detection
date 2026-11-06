library(tidyverse)
library(here)
library(R2jags)
library(coda)
library(ggmcmc)
library(MCMCvis)
library(doParallel)

set.seed(123)

#Create dataset for goodness-of-fit testing
#Set up 5 years (primary occasions) with 5 days (secondary occasions) each
n_iter <- 100
n.ind <- 500
n.years <- 5
n.days <- c(5, 5, 5, 5, 5)

phi.gen <- runif(1, 0.95, 0.99) 
beta.sex.gen <- rnorm(1, 0, 0.1)
mean.p.cap.gen <- runif(1, 0.4, 0.6)
mu.p.cap.gen <- log(mean.p.cap.gen/(1-mean.p.cap.gen))
mean.gamPP.gen <- runif(1, 0.01, 0.9)
mu.gamPP.gen <- log(mean.gamPP.gen/(1-mean.gamPP.gen))
gamP.gen <- runif(1, 0, 1)
eps.year.gen <- rnorm((n.years - 1), 0, 0.1)
eps.year.gen[n.years] <- -sum(eps.year.gen[1:(n.years - 1)])
eps.year.te.gen <- rnorm((n.years - 1), 0, 0.1)
eps.year.te.gen[n.years] <- -sum(eps.year.te.gen[1:(n.years - 1)])
eps.ind.te.gen <- rnorm(n.ind, 0, 0.1)
sex <- sample(c(0, 1), n.ind, replace = TRUE)
Year <- 1:5


#Set up a data generator function to create each dataset
#Each generated dataset is one data frame in an array

det.data.fun <- function(phi, mu.p.cap, mu.gamPP, beta.sex, 
                         eps.year, eps.year.te, gamP, sex,
                         eps.ind.te, n.ind, n.years, n.days){
  
  
  p.cap <- array(NA, dim = c(n.years, max(n.days), n.ind))
  gamPP <- matrix(NA, nrow = n.years, ncol = n.ind)
  for(i in 1:n.ind){
    for(t in 1:n.years){
      for (j in 1:max(n.days[1:n.years])){
        p.cap[t, j, i] <- exp(mu.p.cap + eps.year[t])/(1 + exp(mu.p.cap + eps.year[t]))
        gamPP[t, i] <- exp(mu.gamPP + beta.sex*sex[i] + eps.year.te[t] + eps.ind.te[i])/
          (1 + exp(mu.gamPP + beta.sex*sex[i] + eps.year.te[t] + eps.ind.te[i]))
      }
    }
  }
  
  p.star <- matrix(NA, nrow = n.years, ncol = n.ind)
  for(i in 1:n.ind){
    for(t in 1:n.years){
      p.star[t, i] <- 1 - prod(1 - p.cap[t, 1:n.days[t], i])
    }
  }
  
  #Generate transition probabilities
  #State transition matrices
  #1 is in study area
  #2 out of study area
  #3 is dead
  #gamP and gamPP are opposite of MARK notation (i.e., gamPP is stay, not emigrate, and 
  #gamP is immigrate, not stay out)
  TPM <- array(NA, dim = c(n.years, 3, 3, n.ind))
  OPM <- array(NA, dim = c(n.years, 3, 2, n.ind))
  for(i in 1:n.ind){
    for(t in 1:n.years){
      
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
  
  
  #create capture histories from transition and observation probabilities
  #Get initial states (start all whales at observed in the study area)
  init.state <- rep(1, n.ind)
  
  #Simulate capture histories by first simulating ecological process and then simulating observation process
  pre.ch.p <- matrix(NA, nrow = n.ind, ncol = (n.years + 1))
  ch.p <- matrix(NA, nrow = n.ind, ncol = (n.years + 1))
  pre.ch.p[, 1] <- init.state
  for(t in 1:n.years){ #year for iterating through pre.ch.p because need previous year's state
    for(i in 1:n.ind){
      pre.ch.p[i, t + 1] <- sample(1:3, 1, prob = c(TPM[t, pre.ch.p[i, t], 1, i], TPM[t, pre.ch.p[i, t], 2, i],
                                                    TPM[t, pre.ch.p[i, t], 3, i])) #set prob for each of 3 states
    }
  }
  ch.p[, 1] <- init.state
  for(t in 1:n.years){ #year for iterating through pre.ch.p because need previous year's state
    for(i in 1:n.ind){
      if(pre.ch.p[i, t + 1] == 1){
        ch.p[i, t + 1] <- sample(1:2, 1, prob = c(OPM[t, 1, 1, i], OPM[t, 1, 2, i])) #prob of observed given whale was present
      }else{ch.p[i, t + 1] <- 2} #if not in study area, can only be unobserved
    }
  }
  
  # ch.p <- matrix(NA, nrow = n.ind, ncol = (n.years + 1))
  # ch.p[, 1] <- init.state
  # for(t in 2:(n.years + 1)){ #year for iterating through pre.ch.p because need previous year's state
  #   for(tt in 1:n.years){ #year for probability matrices because need current year
  #     for(i in 1:n.ind){
  #       #Probabilities: state 1: alive, in area, detected, state 2: alive, out, not detected | alive, in, not detected,
  #       #state 3: dead
  #       ch.p[i, t] <- sample(1:3, 1, prob = c(TPM[tt, ch.p[i, t - 1], 1, i] * OPM[tt, 1, 1, i], 
  #                                                 TPM[tt, ch.p[i, t - 1], 2, i] * 1 + 
  #                                                   TPM[tt, ch.p[i, t - 1], 1, i] * OPM[tt, 1, 2, i],
  #                                                 TPM[tt, ch.p[i, t - 1], 3, i]))
  #     }
  #   }
  # }
  
  #Set detections for each day
  # ch.p[ch.p == 3] <- 2
  sim.dat <- ch.p
  
  obs.sim <- array(NA, dim = c(n.ind, n.years, max(n.days)))
  for(i in 1:n.ind){
    for(t in 2:(n.years + 1)){
      for(j in 1:max(n.days)){
        if(sim.dat[i, t] == 1){
          obs.sim[i, t - 1, j] <- sample(0:1, 1, prob = c((1 - p.cap[t - 1, j, i]), p.cap[t - 1, j, i]))
        }else{
          obs.sim[i, t - 1, j] <- 0
        } 
      }
    }
  }
  
  #Find any years in which there was a detection but no days have a detection an replace with detection in last day
  sum.test <- rowSums(obs.sim, dims = 2)
  
  for(i in 1:n.ind){
    for(t in 1:n.years){
      if(sim.dat[i, t + 1] == 1){
        obs.sim[i, t, max(n.days)] <- ifelse(sum.test[i, t] == 0, 1, obs.sim[i, t, max(n.days)])
      }
    }
  }
  
  #Calculate simulated counts per year
  #Number of times individual was detected in secondary occasions per primary occasion
  PrimaryN <- matrix(NA, n.ind, n.years)
  for(i in 1:nrow(PrimaryN)){
    for(j in 1:ncol(PrimaryN)){
      PrimaryN[i,j] <- sum(obs.sim[i,j,], na.rm = TRUE)
    }
  }
  
  #For primary occasions when the individual was seen in at least one secondary occasion, indicate if the animal was seen
  #on the specific secondary occasion to get at detection vs. temporary emigration
  seen <- array(NA, c(n.ind, n.years, max(n.days)))
  missed <- array(NA, c(n.ind, n.years, max(n.days)))
  
  for(i in 1:nrow(PrimaryN)){
    for(t in 1:ncol(PrimaryN)){
      for(j in 1:n.days[t]){
        if(PrimaryN[i,t] > 1 & obs.sim[i,t,j] == 1){seen[i,t,j] <- 1}
        if(PrimaryN[i,t] >= 1 & obs.sim[i,t,j] == 0){missed[i,t,j] <- 1}
      }
    }
  }
  
  #Count number of individuals detected (N.det) and not detected (N.notdet) per secondary occasion by primary occasion
  N.det <- matrix(NA, n.years, max(n.days))
  N.notdet <- matrix(NA, n.years, max(n.days))
  
  for(i in 1:nrow(N.det)){
    for(j in 1:ncol(N.det)){
      N.det[i,j] <- sum(seen[,i,j], na.rm = TRUE)
      N.notdet[i,j] <- sum(missed[,i,j], na.rm = TRUE)
    }
  }
  
  #Make an array with a matrix of the number of individuals detected for each individual (each matrix is identical - just
  #for indexing purposes later)
  N.det <- replicate(n.ind, N.det, simplify = "array")
  N.notdet <- replicate(n.ind, N.notdet, simplify = "array")
  
  CountedSuper <- N.det + N.notdet
  
  sim.list <- list("sim.dat" = sim.dat, "N.det" = N.det, "CountedSuper" = CountedSuper)
  return(sim.list)
  
}


Sim.fun <- function(n_iter, mu.p.cap.gen, phi.gen, beta.sex.gen, mu.gamPP.gen, gamP.gen, Year, n.years, n.days){
  
  
  sim_summ <- as.list(1:n_iter)
  sim.mcmc.list <- as.list(1:n_iter)
  computation_time <- as.list(1:n_iter)
  df_list <- as.list(1:n_iter)
  
  #set progress bar
  pb <- txtProgressBar(min = 0, max = n_iter, style = 3)
  
  #Simulate  model
  for(k in 1:n_iter) {
    
    #Generate data
    sim.list <- det.data.fun(phi = phi.gen, mu.p.cap = mu.p.cap.gen, mu.gamPP = mu.gamPP.gen, beta.sex = beta.sex.gen, 
                             eps.year = eps.year.gen, eps.year.te = eps.year.te.gen, gamP = gamP.gen, sex = sex,
                             eps.ind.te = eps.ind.te.gen, n.ind = n.ind, n.years = n.years, n.days = n.days)
    
    sim.dat.mat <- sim.list[["sim.dat"]]
    sim.dat.mat <- sim.dat.mat[, -1] #remove initial state column from data-generating process
    sim.dat.df <- as.data.frame(sim.dat.mat)
    sim.dat.df$sums <- apply(sim.dat.df, 1, sum)
    sim.dat.df$ogEGNo <- 1:nrow(sim.dat.df)
    sim.dat.df <- sim.dat.df %>% filter(sums != 10) #remove individuals never detected
    sim.dat <- as.matrix(subset(sim.dat.df, select = -c(ogEGNo, sums)))
    N.det <- sim.list[["N.det"]]
    N.det <- N.det[, , c(sim.dat.df$ogEGNo)]
    CountedSuper <- sim.list[["CountedSuper"]]
    CountedSuper <- CountedSuper[, , c(sim.dat.df$ogEGNo)] 
    
    EGNo <- 1:nrow(sim.dat)
    n.ind2 <- nrow(sim.dat)
    sex2 <- sample(c(0, 1), n.ind, replace = TRUE)
    
    FIRST <- function(x) {min((1:n.years)[x == 1])}
    first <- apply(sim.dat, 1, FIRST)
    
    LAST <- function(x) {max((1:n.years)[x > 0 & x < 2])}
    last <- apply(sim.dat, 1, LAST)
    
    z.init <- matrix(NA, nrow(sim.dat), ncol(sim.dat))
    for (i in 1:nrow(sim.dat)){
      if(first[i] < ncol(z.init)){
        z.init[i, (first[i] + 1):ncol(z.init)] <- 1
      }
    }
    
    
    
    ##Run model##
    
    rd.sim.mod.code <- function() {
      
      #Set priors
      
      phi ~ dunif(0, 1)     #Not really interested in estimating phi...more of a nuisance param
      gamP ~ dunif(0, 1)    #gamma' 
      gammaP <- 1 - gamP    #MARK parameterization (prob of staying out)
      beta.sex ~ dnorm(0, 0.1)
      mean.gamPP ~ dunif(0, 1)   #gamma'' 
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
      #1 is in study area (in CCB)
      #2 out of study area (out of CCB)
      #3 is dead
      #gamP and gamPP are opposite of MARK notation (i.e., gamPP is stay, not emigrate, and 
      #gamP is immigrate, not stay out)
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
          #1 down (first position) is in study area (in CCB)
          #2 down (first position) is out of study area (out of CCB)
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
    
    
    params = c("mu.p.cap", "phi", "gammaP", "beta.sex", "mu.gamPP")
    
    
    rd.dat.sim <- list(first = first,
                       last = n.years, #All sim whales end in last year
                       ch = sim.dat, 
                       sex = sex2,
                       EGNo = EGNo,
                       Year = Year,
                       n.secondary = n.days,
                       n.primary = n.years,
                       n.ind = n.ind2,
                       N.det = N.det,
                       CountedSuper = CountedSuper)
    
    inits <- function(){list(z = z.init)}
    
    
    
    
    cl <- makePSOCKcluster(3)
    clusterSetRNGStream(cl)
    registerDoParallel(cl)
    
    
    
    rd.sim.mod <- foreach(i = 1:3, .packages = "R2jags", .combine = 'c') %dopar% {
      load.module("lecuyer")
      jags(data = rd.dat.sim,
           inits = inits,
           parameters.to.save = params,
           model.file = rd.sim.mod.code,
           n.chains = 1,
           n.iter = 10000,
           n.burnin = 5000)
    }
    
    #organize output
    rd.mcmc.sim <- as.mcmc(rd.sim.mod)
    
    test <- rd.sim.mod[[2]]
    sims.mat1 <- test[["sims.matrix"]]
    test <- rd.sim.mod[[8]]
    sims.mat2 <- test[["sims.matrix"]]
    test <- rd.sim.mod[[14]]
    sims.mat3 <- test[["sims.matrix"]]
    sims.list <- list(sims.mat1, sims.mat2, sims.mat3)
    
    mod_mcmc.sim <- as.mcmc.list(lapply(sims.list, function(x) as.mcmc(x)))
    
    
    out_summ <- MCMCsummary(mod_mcmc.sim)
    out_summ <- out_summ[!(row.names(out_summ) %in% c("deviance")), ]
    colnames(out_summ)[c(3, 5)] <- c("Lower95", "Upper95")
    
    # compare posterior summaries to data generating values to track 
    #coverage 
    out_summ$iter <- k
    out_summ$par <- c(rownames(out_summ)) 
    
    #Adjust gamPP to gammaPP
    og_values <- c(beta.sex.gen, (1 - gamP.gen), mu.gamPP.gen, mu.p.cap.gen, phi.gen)  
    
    out_summ$capture <- as.numeric(og_values >= out_summ$Lower95 & 
                                     og_values <= out_summ$Upper95)
    out_summ$cri95_width <- out_summ$Upper95 - out_summ$Lower95
    
    
    #save the summary output from this realization in the list
    sim_summ[[k]] <- out_summ
    sim.mcmc.list[[k]] <- mod_mcmc.sim
    df_list[[k]] <- sim.dat
    setTxtProgressBar(pb, k)
    print(paste("Completed iteration", k, "of", n_iter))
    close(pb)
  }
  
  #combine summaries
  rd_summ <- do.call(rbind, sim_summ)
  out <- list(rd_summ = rd_summ, df_list = df_list)
  
}


####Run the simulation####
sim_rd <- Sim.fun(n_iter = n_iter, mu.p.cap.gen = mu.p.cap.gen, phi.gen = phi.gen, beta.sex.gen = beta.sex.gen, 
                  mu.gamPP.gen = mu.gamPP.gen, gamP.gen = gamP.gen, Year = Year, n.years = n.years, n.days = n.days)


#summarize simulation results 
sim_out <- sim_rd$rd_summ %>% 
  group_by(par) %>% 
  summarise(cap_rate = mean(capture), 
            avg_est = mean(mean), 
            ci_width = mean(Upper95 - Lower95),  
            avg_L95 = mean(Lower95), 
            avg_U95 = mean(Upper95),
            avg_rhat = mean(Rhat), 
            avg_neff = mean(n.eff))

#Calculate biases and summarize across first detections
sim_out$og_values <- c(beta.sex.gen, (1 - gamP.gen), mu.gamPP.gen, mu.p.cap.gen, phi.gen)
sim_out$avg_bias <- sim_out$avg_est - sim_out$og_values
avg.bias <- mean(sim_out$avg_bias)
avg.capture <- mean(sim_out$cap_rate) 



