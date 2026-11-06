#Hierarchical logistic regression model
library(tidyverse)
library(here)
library(coda)
library(ggmcmc)
library(MCMCvis)


#Bring in aerial survey data
aerial1 <- read.csv("AerialData1.csv")
aerial2 <- read.csv("AerialData2.csv")
aerial3 <- read.csv("AerialData3.csv")
aerial4 <- read.csv("AerialData4.csv")

aerial.data <- rbind(aerial1, aerial2, aerial3, aerial4)


#############
##Stan code##
#############

library(rstanarm)
library(bayesplot)
library(HDInterval)
library(shinystan)
library(cowplot)

#Assume that non-detections are because the animal was not present in the survey area, in the area but hidden or unavailable
#for detection, or present and available but missed by the observer
#SNE (southern New England) takes on the intercept

#Intercept estimates SNE zone
det.mod <- stan_glmer(Detected ~ w.grp.sz.sq + w.grp.sz + BOF + CCB + GB + GOM + GSC + GSL + JL + MIDA + RB + 
                        SEUS + SAG + (1 | Year) + (1 | survey) + (1 | SightingEGNo),
                      data = aerial.data, binomial(link = "logit"), prior = normal(0, 1),
                      prior_intercept = normal(-1.5, 1), cores = 6)



#Look at results
print(det.mod, digits = 4)


#Diagnostics
det.samps.temp <- as.array(det.mod)

det.mcmc <- do.call(mcmc.list, plyr:::alply(det.samps.temp[, , ], 2, as.mcmc))
det.summ.full <- MCMCsummary(det.mcmc)


mcmc_params <- as.mcmc.list(lapply(det.mcmc, function(x) as.mcmc(x[,!grepl("b", colnames(x))])))
det.summ <- MCMCsummary(mcmc_params)
det.ggs <- ggs(mcmc_params)
ggs_geweke(det.ggs)
ggs_Rhat(det.ggs)


det.ggs$Parameter <- gsub("w.grp.sz.sq", "group size squared", det.ggs$Parameter)
det.ggs$Parameter <- gsub("w.grp.sz", "group size", det.ggs$Parameter)
det.ggs$Parameter <- gsub("SAG", "Surface active group", det.ggs$Parameter)
det.ggs$Parameter[str_sub(det.ggs$Parameter, start = 1, end = 12) == "Sigma[survey"] <- "Sigma survey"
det.ggs$Parameter[str_sub(det.ggs$Parameter, start = 1, end = 11) == "Sigma[Sight"] <- "Sigma individual"
det.ggs$Parameter[str_sub(det.ggs$Parameter, start = 1, end = 10) == "Sigma[Year"] <- "Sigma year"
det.ggs$Parameter[str_sub(det.ggs$Parameter, start = 2, end = 10) == "Intercept"] <- "SNE"
det.ggs$Parameter <- gsub("(Intercept)", "SNE", det.ggs$Parameter)

det.ggs$Parameter <- as.factor(det.ggs$Parameter)


ggs_traceplot(det.ggs, c("(Sigma)"))
ggs_traceplot(det.ggs, c("group"))
ggs_traceplot(det.ggs, c("SNE"))


#Post-predictive checks (need to limit draws because dataset is so huge).
NARW.dat$AreaCode2 <- NARW.dat.full$AreaCode2
pr.det <- posterior_predict(det.mod, draws = 50, type = "response", re.form = NA)
y = NARW.dat$Detected
ppc_stat(y, pr.det, stat = "median")
ppc_bars(y, pr.det, freq = FALSE)

#summarize by zone
ppc_stat_grouped(NARW.dat$Detected,
                 pr.det,
                 group = NARW.dat$AreaCode2, 
                 stat = "median")

pp_check(det.mod, plotfun = "hist", nreps = 1)

#Calculate residuals by geographic zone
obs <- cbind(NARW.dat$Detected, NARW.dat$AreaCode2)
obs <- as.data.frame(obs)
colnames(obs) <- c("Detected", "Zone")

obs.BOF <- obs %>% filter(Zone == "BOF")
obs.CCB <- obs %>% filter(Zone == "CCB")
obs.GB <- obs %>% filter(Zone == "GB")
obs.GOM <- obs %>% filter(Zone == "GOM")
obs.GSC <- obs %>% filter(Zone == "GSC")
obs.GSL <- obs %>% filter(Zone == "GSL")
obs.JL <- obs %>% filter(Zone == "JL")
obs.MIDA <- obs %>% filter(Zone == "MIDA")
obs.RB <- obs %>% filter(Zone == "RB")
obs.SEUS <- obs %>% filter(Zone == "SEUS")
obs.SNE <- obs %>% filter(Zone == "SNE")

b.det = data.frame(as.matrix(det.mcmc[, c(1:15)]))
set.seed(123)
sim.zone <- runif(4000, 1, 11)
sim.zone <- round(sim.zone, digits = 0)
sim.zone <- as.data.frame(sim.zone)
sim.zone$sim.zone <- ifelse(sim.zone$sim.zone == 1, "BOF", ifelse(sim.zone$sim.zone == 2, "CCB", 
                     ifelse(sim.zone$sim.zone == 3, "GB", ifelse(sim.zone$sim.zone == 4, "GOM", 
                     ifelse(sim.zone$sim.zone == 5, "GSC", ifelse(sim.zone$sim.zone == 6, "GSL", 
                     ifelse(sim.zone$sim.zone == 7, "JL", ifelse(sim.zone$sim.zone == 8, "MIDA", 
                     ifelse(sim.zone$sim.zone == 9, "RB", ifelse(sim.zone$sim.zone == 10, "SEUS", "SNE"))))))))))

sim.zone$sim.zone <- factor(sim.zone$sim.zone)
mm = with(sim.zone, model.matrix(~sim.zone, contrasts.arg = list(sim.zone = "contr.sum")))
dm.z <- mm[ , 2:ncol(mm)]
colnames(dm.z) <- c("BOF", "CCB", "GB", "GOM", "GSC", "GSL", "JL", "MIDA", "RB", "SEUS")
sim.zone <- cbind(sim.zone, dm.z)

post_sim <- as.data.frame(matrix(nrow = 4000, ncol = 2))
for(i in 1:nrow(sim.zone)){
  post_sim$Zone[i] <- sim.zone$sim.zone[i]
  post_sim$Detected[i] <- plogis((b.det$X.Intercept.[i] + b.det$w.grp.sz[i] + b.det$w.grp.sz.sq[i] +
                                    sim.zone$BOF[i]* b.det$BOF[i] +
                                    sim.zone$CCB[i]* b.det$CCB[i] +
                                    sim.zone$GB[i]* b.det$GB[i] +
                                    sim.zone$GOM[i]* b.det$GOM[i] +
                                    sim.zone$GSC[i]* b.det$GSC[i] +
                                    sim.zone$GSL[i]* b.det$GSL[i] +
                                    sim.zone$JL[i]* b.det$JL[i] +
                                    sim.zone$MIDA[i]* b.det$MIDA[i] +
                                    sim.zone$RB[i]* b.det$RB[i] +
                                    sim.zone$SEUS[i]* b.det$SEUS[i]))
}

post_sim <- post_sim[, -c(1:2)]

post_sim.BOF <- post_sim %>% filter(Zone == 1)
post_sim.CCB <- post_sim %>% filter(Zone == 2)
post_sim.GB <- post_sim %>% filter(Zone == 3)
post_sim.GOM <- post_sim %>% filter(Zone == 4)
post_sim.GSC <- post_sim %>% filter(Zone == 5)
post_sim.GSL <- post_sim %>% filter(Zone == 6)
post_sim.JL <- post_sim %>% filter(Zone == 7)
post_sim.MIDA <- post_sim %>% filter(Zone == 8)
post_sim.RB <- post_sim %>% filter(Zone == 9)
post_sim.SEUS <- post_sim %>% filter(Zone == 10)
post_sim.SNE <- post_sim %>% filter(Zone == 11)

#Get residuals
#BOF
obs.BOF <- mean(as.numeric(obs.BOF$Detected))
residuals.BOF <- data.frame(matrix(nrow = 500, ncol = 1))
for(i in 1:nrow(post_sim.BOF)){
  residuals.BOF[i,] <- obs.BOF - post_sim.BOF$Detected[i]
}
colnames(residuals.BOF) <- c("BOF")
residuals.BOF <- residuals.BOF[complete.cases(residuals.BOF),]
residuals.BOF <- as.data.frame(residuals.BOF)

#CCB
obs.CCB <- mean(as.numeric(obs.CCB$Detected))
residuals.CCB <- data.frame(matrix(nrow = 500, ncol = 1))
for(i in 1:nrow(post_sim.CCB)){
  residuals.CCB[i,] <- obs.CCB - post_sim.CCB$Detected[i]
}
colnames(residuals.CCB) <- c("CCB")
residuals.CCB <- residuals.CCB[complete.cases(residuals.CCB),]
residuals.CCB <- as.data.frame(residuals.CCB)

#GB
obs.GB <- mean(as.numeric(obs.GB$Detected))
residuals.GB <- data.frame(matrix(nrow = 500, ncol = 1))
for(i in 1:nrow(post_sim.GB)){
  residuals.GB[i,] <- obs.GB - post_sim.GB$Detected[i]
}
colnames(residuals.GB) <- c("GB")
residuals.GB <- residuals.GB[complete.cases(residuals.GB),]
residuals.GB <- as.data.frame(residuals.GB)

#GOM
obs.GOM <- mean(as.numeric(obs.GOM$Detected))
residuals.GOM <- data.frame(matrix(nrow = 500, ncol = 1))
for(i in 1:nrow(post_sim.GOM)){
  residuals.GOM[i,] <- obs.GOM - post_sim.GOM$Detected[i]
}
colnames(residuals.GOM) <- c("GOM")
residuals.GOM <- residuals.GOM[complete.cases(residuals.GOM),]
residuals.GOM <- as.data.frame(residuals.GOM)

#GSC
obs.GSC <- mean(as.numeric(obs.GSC$Detected))
residuals.GSC <- data.frame(matrix(nrow = 500, ncol = 1))
for(i in 1:nrow(post_sim.GSC)){
  residuals.GSC[i,] <- obs.GSC - post_sim.GSC$Detected[i]
}
colnames(residuals.GSC) <- c("GSC")
residuals.GSC <- residuals.GSC[complete.cases(residuals.GSC),]
residuals.GSC <- as.data.frame(residuals.GSC)

#GSL
obs.GSL <- mean(as.numeric(obs.GSL$Detected))
residuals.GSL <- data.frame(matrix(nrow = 500, ncol = 1))
for(i in 1:nrow(post_sim.GSL)){
  residuals.GSL[i,] <- obs.GSL - post_sim.GSL$Detected[i]
}
colnames(residuals.GSL) <- c("GSL")
residuals.GSL <- residuals.GSL[complete.cases(residuals.GSL),]
residuals.GSL <- as.data.frame(residuals.GSL)

#JL
obs.JL <- mean(as.numeric(obs.JL$Detected))
residuals.JL <- data.frame(matrix(nrow = 500, ncol = 1))
for(i in 1:nrow(post_sim.JL)){
  residuals.JL[i,] <- obs.JL - post_sim.JL$Detected[i]
}
colnames(residuals.JL) <- c("JL")
residuals.JL <- residuals.JL[complete.cases(residuals.JL),]
residuals.JL <- as.data.frame(residuals.JL)

#MIDA
obs.MIDA <- mean(as.numeric(obs.MIDA$Detected))
residuals.MIDA <- data.frame(matrix(nrow = 500, ncol = 1))
for(i in 1:nrow(post_sim.MIDA)){
  residuals.MIDA[i,] <- obs.MIDA - post_sim.MIDA$Detected[i]
}
colnames(residuals.MIDA) <- c("MIDA")
residuals.MIDA <- residuals.MIDA[complete.cases(residuals.MIDA),]
residuals.MIDA <- as.data.frame(residuals.MIDA)

#RB
obs.RB <- mean(as.numeric(obs.RB$Detected))
residuals.RB <- data.frame(matrix(nrow = 500, ncol = 1))
for(i in 1:nrow(post_sim.RB)){
  residuals.RB[i,] <- obs.RB - post_sim.RB$Detected[i]
}
colnames(residuals.RB) <- c("RB")
residuals.RB <- residuals.RB[complete.cases(residuals.RB),]
residuals.RB <- as.data.frame(residuals.RB)

#SEUS
obs.SEUS <- mean(as.numeric(obs.SEUS$Detected))
residuals.SEUS <- data.frame(matrix(nrow = 500, ncol = 1))
for(i in 1:nrow(post_sim.SEUS)){
  residuals.SEUS[i,] <- obs.SEUS - post_sim.SEUS$Detected[i]
}
colnames(residuals.SEUS) <- c("SEUS")
residuals.SEUS <- residuals.SEUS[complete.cases(residuals.SEUS),]
residuals.SEUS <- as.data.frame(residuals.SEUS)

#SNE
obs.SNE <- mean(as.numeric(obs.SNE$Detected))
residuals.SNE <- data.frame(matrix(nrow = 500, ncol = 1))
for(i in 1:nrow(post_sim.SNE)){
  residuals.SNE[i,] <- obs.SNE - post_sim.SNE$Detected[i]
}
colnames(residuals.SNE) <- c("SNE")
residuals.SNE <- residuals.SNE[complete.cases(residuals.SNE),]
residuals.SNE <- as.data.frame(residuals.SNE)


#Make plots
BOF <- ggplot(residuals.BOF) + 
  geom_histogram(aes(x = residuals.BOF), color = "black", fill = "white", binwidth = 0.0003) + 
  # scale_x_continuous(limits = c(-0.04, 0.01), breaks = seq(from = -0.04, to = 0.01, by = 0.01)) +
  scale_x_continuous(limits = c(-0.04, 0.04), breaks = seq(from = -0.04, to = 0.04, by = 0.04)) +
  ylab("Frequency") +
  xlab("BOF residuals") + 
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.text.x = element_text(size = 4),
        axis.title.y = element_text(size = 5),
        axis.title.x = element_text(size = 5)) +
  theme(axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

CCB <- ggplot(residuals.CCB) + 
  geom_histogram(aes(x = residuals.CCB), color = "black", fill = "white", binwidth = 0.0003) + 
  # scale_x_continuous(limits = c(-0.005, 0.005), breaks = seq(from = -0.005, to = 0.005, by = 0.005)) +
  scale_x_continuous(limits = c(-0.01, 0.01), breaks = seq(from = -0.01, to = 0.01, by = 0.01)) +
  ylab("Frequency") +
  xlab("CCB residuals") + 
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.text.x = element_text(size = 4),
        axis.title.y = element_text(size = 5),
        axis.title.x = element_text(size = 5)) +
  theme(axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

GB <- ggplot(residuals.GB) + 
  geom_histogram(aes(x = residuals.GB), color = "black", fill = "white", binwidth = 0.0003) + 
  # scale_x_continuous(limits = c(-0.005, 0.005), breaks = seq(from = -0.005, to = 0.005, by = 0.002)) +
  scale_x_continuous(limits = c(-0.01, 0.01), breaks = seq(from = -0.01, to = 0.01, by = 0.01)) +
  ylab("Frequency") +
  xlab("GB residuals") + 
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.text.x = element_text(size = 4),
        axis.title.y = element_text(size = 5),
        axis.title.x = element_text(size = 5)) +
  theme(axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

GOM <- ggplot(residuals.GOM) + 
  geom_histogram(aes(x = residuals.GOM), color = "black", fill = "white", binwidth = 0.0003) + 
  # scale_x_continuous(limits = c(-0.005, 0.005), breaks = seq(from = -0.005, to = 0.005, by = 0.005)) +
  scale_x_continuous(limits = c(-0.01, 0.01), breaks = seq(from = -0.01, to = 0.01, by = 0.01)) +
  ylab("Frequency") +
  xlab("GOM residuals") + 
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.text.x = element_text(size = 4),
        axis.title.y = element_text(size = 5),
        axis.title.x = element_text(size = 5)) +
  theme(axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

GSC <- ggplot(residuals.GSC) + 
  geom_histogram(aes(x = residuals.GSC), color = "black", fill = "white", binwidth = 0.0003) + 
  # scale_x_continuous(limits = c(-0.005, 0.005), breaks = seq(from = -0.005, to = 0.005, by = 0.005)) +
  scale_x_continuous(limits = c(-0.01, 0.01), breaks = seq(from = -0.01, to = 0.01, by = 0.01)) +
  ylab("Frequency") +
  xlab("GSC residuals") + 
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.text.x = element_text(size = 4),
        axis.title.y = element_text(size = 5),
        axis.title.x = element_text(size = 5)) +
  theme(axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

GSL <- ggplot(residuals.GSL) + 
  geom_histogram(aes(x = residuals.GSL), color = "black", fill = "white", binwidth = 0.0003) + 
  # scale_x_continuous(limits = c(-0.012, 0.012), breaks = seq(from = -0.012, to = 0.012, by = 0.01)) +
  scale_x_continuous(limits = c(-0.01, 0.01), breaks = seq(from = -0.01, to = 0.01, by = 0.01)) +
  ylab("Frequency") +
  xlab("GSL residuals") + 
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.text.x = element_text(size = 4),
        axis.title.y = element_text(size = 5),
        axis.title.x = element_text(size = 5)) +
  theme(axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

JL <- ggplot(residuals.JL) + 
  geom_histogram(aes(x = residuals.JL), color = "black", fill = "white", binwidth = 0.0003) + 
  # scale_x_continuous(limits = c(-0.005, 0.008), breaks = seq(from = -0.005, to = 0.008, by = 0.005)) +
  scale_x_continuous(limits = c(-0.01, 0.01), breaks = seq(from = -0.01, to = 0.01, by = 0.01)) +
  ylab("Frequency") +
  xlab("JL residuals") + 
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.text.x = element_text(size = 4),
        axis.title.y = element_text(size = 5),
        axis.title.x = element_text(size = 5)) +
  theme(axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

MIDA <- ggplot(residuals.MIDA) + 
  geom_histogram(aes(x = residuals.MIDA), color = "black", fill = "white", binwidth = 0.0003) + 
  # scale_x_continuous(limits = c(-0.006, 0.004), breaks = seq(from = -0.006, to = 0.004, by = 0.005)) +
  scale_x_continuous(limits = c(-0.01, 0.01), breaks = seq(from = -0.01, to = 0.01, by = 0.01)) +
  ylab("Frequency") +
  xlab("MIDA residuals") + 
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.text.x = element_text(size = 4),
        axis.title.y = element_text(size = 5),
        axis.title.x = element_text(size = 5)) +
  theme(axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

RB <- ggplot(residuals.RB) + 
  geom_histogram(aes(x = residuals.RB), color = "black", fill = "white", binwidth = 0.0003) + 
  # scale_x_continuous(limits = c(-0.007, 0.003), breaks = seq(from = -0.007, to = 0.003, by = 0.005)) +
  scale_x_continuous(limits = c(-0.01, 0.01), breaks = seq(from = -0.01, to = 0.01, by = 0.01)) +
  ylab("Frequency") +
  xlab("RB residuals") + 
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.text.x = element_text(size = 4),
        axis.title.y = element_text(size = 5),
        axis.title.x = element_text(size = 5)) +
  theme(axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

SEUS <- ggplot(residuals.SEUS) + 
  geom_histogram(aes(x = residuals.SEUS), color = "black", fill = "white", binwidth = 0.0003) + 
  # scale_x_continuous(limits = c(-0.006, 0.004), breaks = seq(from = -0.006, to = 0.004, by = 0.005)) +
  scale_x_continuous(limits = c(-0.01, 0.01), breaks = seq(from = -0.01, to = 0.01, by = 0.01)) +
  ylab("Frequency") +
  xlab("SEUS residuals") + 
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.text.x = element_text(size = 4),
        axis.title.y = element_text(size = 5),
        axis.title.x = element_text(size = 5)) +
  theme(axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

SNE <- ggplot(residuals.SNE) + 
  geom_histogram(aes(x = residuals.SNE), color = "black", fill = "white", binwidth = 0.0003) + 
  # scale_x_continuous(limits = c(-0.006, 0.004), breaks = seq(from = -0.006, to = 0.004, by = 0.005)) +
  scale_x_continuous(limits = c(-0.01, 0.01), breaks = seq(from = -0.01, to = 0.01, by = 0.01)) +
  ylab("Frequency") +
  xlab("SNE residuals") + 
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.text.x = element_text(size = 4),
        axis.title.y = element_text(size = 5),
        axis.title.x = element_text(size = 5)) +
  theme(axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

resids <- plot_grid(BOF, CCB, GB, GOM, GSC, GSL, JL, MIDA, RB, SEUS, SNE)


#################
###Vessel data###
#################

#Bring in aerial survey data
vessel <- read.csv("VesselData.csv")

#############
##Stan code##
#############

#Assume that non-detections are because the animal was not present in the survey area, in the area but hidden or unavailable
#for detection, or present and available but missed by the observer
#Intercept estimates SNE (southern New England) zone

det.mod <- stan_glmer(Detected ~ w.grp.sz.sq + w.grp.sz + BOF + CCB + GB + GOM + GSL + RB + 
                        SAG + (1 | Year) + (1 | survey) + (1 | SightingEGNo),
                      data = NARW.dat, binomial(link = "logit"), prior = normal(0, 1),
                      prior_intercept = normal(-1.5, 1), cores = 6)

#Look at results
print(det.mod, digits = 4)

#Diagnostics
det.samps.temp <- as.array(det.mod)

det.mcmc <- do.call(mcmc.list, plyr:::alply(det.samps.temp[, , ], 2, as.mcmc))
det.summ.full <- MCMCsummary(det.mcmc)


mcmc_params <- as.mcmc.list(lapply(det.mcmc, function(x) as.mcmc(x[,!grepl("b", colnames(x))])))
det.summ <- MCMCsummary(mcmc_params)
det.ggs <- ggs(mcmc_params)
ggs_geweke(det.ggs)
ggs_Rhat(det.ggs)


det.ggs$Parameter <- gsub("w.grp.sz.sq", "group size squared", det.ggs$Parameter)
det.ggs$Parameter <- gsub("w.grp.sz", "group size", det.ggs$Parameter)
det.ggs$Parameter <- gsub("SAG", "Surface active group", det.ggs$Parameter)
det.ggs$Parameter[str_sub(det.ggs$Parameter, start = 1, end = 12) == "Sigma[survey"] <- "Sigma survey"
det.ggs$Parameter[str_sub(det.ggs$Parameter, start = 1, end = 11) == "Sigma[Sight"] <- "Sigma individual"
det.ggs$Parameter[str_sub(det.ggs$Parameter, start = 1, end = 10) == "Sigma[Year"] <- "Sigma year"
det.ggs$Parameter[str_sub(det.ggs$Parameter, start = 2, end = 10) == "Intercept"] <- "SNE"
det.ggs$Parameter <- gsub("(Intercept)", "SNE", det.ggs$Parameter)

det.ggs$Parameter <- as.factor(det.ggs$Parameter)


ggs_traceplot(det.ggs, c("(Sigma)"))
ggs_traceplot(det.ggs, c("group"))
ggs_traceplot(det.ggs, c("SNE"))


#Post-predictive checks (need to limit draws because dataset is so huge).
#Also might want to remove large files from data prep
NARW.dat$AreaCode2 <- NARW.dat.full$AreaCode2
pr.det <- posterior_predict(det.mod, draws = 50, type = "response", re.form = NA)
y = NARW.dat$Detected
ppc_stat(y, pr.det, stat = "median")
ppc_bars(y, pr.det, freq = FALSE)

#summarize by zone
ppc_stat_grouped(NARW.dat$Detected,
                 pr.det,
                 group = NARW.dat$AreaCode2, 
                 stat = "median")

pp_check(det.mod, plotfun = "hist", nreps = 1)

#Calculate residuals by geographic zone
#Find which geographic zones were only surveyed by ship
ships <- NARW.dat.full %>% filter(PlatformCode == "S")
unique(ships$Zone)

obs <- cbind(NARW.dat$Detected, NARW.dat$AreaCode2)
obs <- as.data.frame(obs)
colnames(obs) <- c("Detected", "Zone")

obs.BOF <- obs %>% filter(Zone == "BOF")
obs.CCB <- obs %>% filter(Zone == "CCB")
obs.GB <- obs %>% filter(Zone == "GB")
obs.GOM <- obs %>% filter(Zone == "GOM")
obs.GSL <- obs %>% filter(Zone == "GSL")
obs.RB <- obs %>% filter(Zone == "RB")
obs.SNE <- obs %>% filter(Zone == "SNE")

b.det = data.frame(as.matrix(det.mcmc[, c(1:10)]))
set.seed(123)
sim.zone <- runif(4000, 1, 7)
sim.zone <- round(sim.zone, digits = 0)
sim.zone <- as.data.frame(sim.zone)
sim.zone$sim.zone <- ifelse(sim.zone$sim.zone == 1, "BOF", ifelse(sim.zone$sim.zone == 2, "CCB", 
                                                                  ifelse(sim.zone$sim.zone == 3, "GB", ifelse(sim.zone$sim.zone == 4, "GOM", 
                                                                                                              ifelse(sim.zone$sim.zone == 5, "GSL", ifelse(sim.zone$sim.zone == 6, "RB", "SNE"))))))

sim.zone$sim.zone <- factor(sim.zone$sim.zone)
mm = with(sim.zone, model.matrix(~sim.zone, contrasts.arg = list(sim.zone = "contr.sum")))
dm.z <- mm[ , 2:ncol(mm)]
colnames(dm.z) <- c("BOF", "CCB", "GB", "GOM", "GSL", "RB")
sim.zone <- cbind(sim.zone, dm.z)

post_sim <- as.data.frame(matrix(nrow = 4000, ncol = 2))
for(i in 1:nrow(sim.zone)){
  post_sim$Zone[i] <- sim.zone$sim.zone[i]
  post_sim$Detected[i] <- plogis((b.det$X.Intercept.[i] + b.det$w.grp.sz[i] + b.det$w.grp.sz.sq[i] +
                                    sim.zone$BOF[i]* b.det$BOF[i] +
                                    sim.zone$CCB[i]* b.det$CCB[i] +
                                    sim.zone$GB[i]* b.det$GB[i] +
                                    sim.zone$GOM[i]* b.det$GOM[i] +
                                    sim.zone$GSL[i]* b.det$GSL[i] +
                                    sim.zone$RB[i]* b.det$RB[i]))
}

post_sim <- post_sim[, -c(1:2)]

post_sim.BOF <- post_sim %>% filter(Zone == 1)
post_sim.CCB <- post_sim %>% filter(Zone == 2)
post_sim.GB <- post_sim %>% filter(Zone == 3)
post_sim.GOM <- post_sim %>% filter(Zone == 4)
post_sim.GSL <- post_sim %>% filter(Zone == 5)
post_sim.RB <- post_sim %>% filter(Zone == 6)
post_sim.SNE <- post_sim %>% filter(Zone == 7)

#Get residuals
#BOF
obs.BOF <- mean(as.numeric(obs.BOF$Detected))
residuals.BOF <- data.frame(matrix(nrow = 500, ncol = 1))
for(i in 1:nrow(post_sim.BOF)){
  residuals.BOF[i,] <- obs.BOF - post_sim.BOF$Detected[i]
}
colnames(residuals.BOF) <- c("BOF")
residuals.BOF <- residuals.BOF[complete.cases(residuals.BOF),]
residuals.BOF <- as.data.frame(residuals.BOF)

#CCB
obs.CCB <- mean(as.numeric(obs.CCB$Detected))
residuals.CCB <- data.frame(matrix(nrow = 500, ncol = 1))
for(i in 1:nrow(post_sim.CCB)){
  residuals.CCB[i,] <- obs.CCB - post_sim.CCB$Detected[i]
}
colnames(residuals.CCB) <- c("CCB")
residuals.CCB <- residuals.CCB[complete.cases(residuals.CCB),]
residuals.CCB <- as.data.frame(residuals.CCB)

#GB
obs.GB <- mean(as.numeric(obs.GB$Detected))
residuals.GB <- data.frame(matrix(nrow = 500, ncol = 1))
for(i in 1:nrow(post_sim.GB)){
  residuals.GB[i,] <- obs.GB - post_sim.GB$Detected[i]
}
colnames(residuals.GB) <- c("GB")
residuals.GB <- residuals.GB[complete.cases(residuals.GB),]
residuals.GB <- as.data.frame(residuals.GB)

#GOM
obs.GOM <- mean(as.numeric(obs.GOM$Detected))
residuals.GOM <- data.frame(matrix(nrow = 500, ncol = 1))
for(i in 1:nrow(post_sim.GOM)){
  residuals.GOM[i,] <- obs.GOM - post_sim.GOM$Detected[i]
}
colnames(residuals.GOM) <- c("GOM")
residuals.GOM <- residuals.GOM[complete.cases(residuals.GOM),]
residuals.GOM <- as.data.frame(residuals.GOM)


#GSL
obs.GSL <- mean(as.numeric(obs.GSL$Detected))
residuals.GSL <- data.frame(matrix(nrow = 500, ncol = 1))
for(i in 1:nrow(post_sim.GSL)){
  residuals.GSL[i,] <- obs.GSL - post_sim.GSL$Detected[i]
}
colnames(residuals.GSL) <- c("GSL")
residuals.GSL <- residuals.GSL[complete.cases(residuals.GSL),]
residuals.GSL <- as.data.frame(residuals.GSL)

#RB
obs.RB <- mean(as.numeric(obs.RB$Detected))
residuals.RB <- data.frame(matrix(nrow = 500, ncol = 1))
for(i in 1:nrow(post_sim.RB)){
  residuals.RB[i,] <- obs.RB - post_sim.RB$Detected[i]
}
colnames(residuals.RB) <- c("RB")
residuals.RB <- residuals.RB[complete.cases(residuals.RB),]
residuals.RB <- as.data.frame(residuals.RB)

#SNE
obs.SNE <- mean(as.numeric(obs.SNE$Detected))
residuals.SNE <- data.frame(matrix(nrow = 500, ncol = 1))
for(i in 1:nrow(post_sim.SNE)){
  residuals.SNE[i,] <- obs.SNE - post_sim.SNE$Detected[i]
}
colnames(residuals.SNE) <- c("SNE")
residuals.SNE <- residuals.SNE[complete.cases(residuals.SNE),]
residuals.SNE <- as.data.frame(residuals.SNE)


#Make plots
BOF <- ggplot(residuals.BOF) + 
  geom_histogram(aes(x = residuals.BOF), color = "black", fill = "white", binwidth = 0.0003) + 
  scale_x_continuous(limits = c(-0.1, 0.1), breaks = seq(from = -0.1, to = 0.1, by = 0.1)) +
  ylab("Frequency") +
  xlab("BOF residuals") + 
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.text.x = element_text(size = 4),
        axis.title.y = element_text(size = 5),
        axis.title.x = element_text(size = 5)) +
  theme(axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

CCB <- ggplot(residuals.CCB) + 
  geom_histogram(aes(x = residuals.CCB), color = "black", fill = "white", binwidth = 0.0003) + 
  # scale_x_continuous(limits = c(-0.05, 0.05), breaks = seq(from = -0.05, to = 0.05, by = 0.05)) +
  scale_x_continuous(limits = c(-0.1, 0.1), breaks = seq(from = -0.1, to = 0.1, by = 0.1)) +
  ylab("Frequency") +
  xlab("CCB residuals") + 
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.text.x = element_text(size = 4),
        axis.title.y = element_text(size = 5),
        axis.title.x = element_text(size = 5)) +
  theme(axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

GB <- ggplot(residuals.GB) + 
  geom_histogram(aes(x = residuals.GB), color = "black", fill = "white", binwidth = 0.0003) + 
  # scale_x_continuous(limits = c(-0.03, 0.07), breaks = seq(from = -0.03, to = 0.07, by = 0.05)) +
  scale_x_continuous(limits = c(-0.1, 0.1), breaks = seq(from = -0.1, to = 0.1, by = 0.1)) +
  ylab("Frequency") +
  xlab("GB residuals") + 
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.text.x = element_text(size = 4),
        axis.title.y = element_text(size = 5),
        axis.title.x = element_text(size = 5)) +
  theme(axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

GOM <- ggplot(residuals.GOM) + 
  geom_histogram(aes(x = residuals.GOM), color = "black", fill = "white", binwidth = 0.0003) + 
  # scale_x_continuous(limits = c(-0.02, 0.02), breaks = seq(from = -0.02, to = 0.02, by = 0.02)) +
  scale_x_continuous(limits = c(-0.1, 0.1), breaks = seq(from = -0.1, to = 0.1, by = 0.1)) +
  ylab("Frequency") +
  xlab("GOM residuals") + 
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.text.x = element_text(size = 4),
        axis.title.y = element_text(size = 5),
        axis.title.x = element_text(size = 5)) +
  theme(axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

GSL <- ggplot(residuals.GSL) + 
  geom_histogram(aes(x = residuals.GSL), color = "black", fill = "white", binwidth = 0.0003) + 
  # scale_x_continuous(limits = c(-0.03, 0.03), breaks = seq(from = -0.03, to = 0.03, by = 0.03)) +
  scale_x_continuous(limits = c(-0.1, 0.1), breaks = seq(from = -0.1, to = 0.1, by = 0.1)) +
  ylab("Frequency") +
  xlab("GSL residuals") + 
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.text.x = element_text(size = 4),
        axis.title.y = element_text(size = 5),
        axis.title.x = element_text(size = 5)) +
  theme(axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

RB <- ggplot(residuals.RB) + 
  geom_histogram(aes(x = residuals.RB), color = "black", fill = "white", binwidth = 0.0003) + 
  # scale_x_continuous(limits = c(-0.008, 0.004), breaks = seq(from = -0.008, to = 0.004, by = 0.004)) +
  scale_x_continuous(limits = c(-0.1, 0.1), breaks = seq(from = -0.1, to = 0.1, by = 0.1)) +
  ylab("Frequency") +
  xlab("RB residuals") + 
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.text.x = element_text(size = 4),
        axis.title.y = element_text(size = 5),
        axis.title.x = element_text(size = 5)) +
  theme(axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

SNE <- ggplot(residuals.SNE) + 
  geom_histogram(aes(x = residuals.SNE), color = "black", fill = "white", binwidth = 0.0003) + 
  # scale_x_continuous(limits = c(-0.02, 0.02), breaks = seq(from = -0.02, to = 0.02, by = 0.02)) +
  scale_x_continuous(limits = c(-0.1, 0.1), breaks = seq(from = -0.1, to = 0.1, by = 0.1)) +
  ylab("Frequency") +
  xlab("SNE residuals") + 
  theme_bw() +
  theme(plot.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(legend.position = "none") +
  theme(panel.border = element_blank(),
        axis.text.y = element_text(size = 4),
        axis.text.x = element_text(size = 4),
        axis.title.y = element_text(size = 5),
        axis.title.x = element_text(size = 5)) +
  theme(axis.line.x = element_line(color = "black", linewidth = 0.5),
        axis.line.y = element_line(color = "black", linewidth = 0.5))

resids <- plot_grid(BOF, CCB, GB, GOM, GSL, RB, SNE)





