###### Load packages
library(glmmTMB)
library(tidyverse)
library(DHARMa)
library(cowplot)
library(AICcmodavg)
library(DT)

###### Read in data
dater <- read.csv("All_Oysters_2025-05 csv.csv") %>% mutate(FixedLocationID = factor(FixedLocationID), StationName = factor(StationName))

###### Take a quick look at the data
glimpse(dater)

###### Summarize data: count up number of samples, total number of legal oysters, mean, SD, lower 95% CL, upper 95% CL, and CV of the number of legal oysters, and finally the numbe rof 0 observations by FixedLocationID
datSumm <- dater %>% group_by(FixedLocationID) %>% summarize(nSamples = n(), TotalLegal = sum(NumLegal), MeanLegal = round(mean(NumLegal),2), sdLegal = round(sd(NumLegal),2), lwr = round(quantile(NumLegal, 0.025),2), upr = round(quantile(NumLegal, 0.975),2), CVLegal = round(sdLegal/MeanLegal*100,2), prop0 = round(mean(NumLegal==0),2))

###### Look at the summaries
datatable(datSumm)

###### Model-fitting by FixedLocationID: there are 4 models, 2 nbinom and 2 poisson, one with and one without zero-inflation
dater391 <- dater %>% filter(FixedLocationID=="391")
m1 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~1, family = nbinom2, data = dater391)
m2 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~0, family = nbinom2, data = dater391)
m3 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~1, family = poisson, data = dater391)
m4 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~0, family = poisson, data = dater391)

###### Assess goodness-of-fit
simulateResiduals(m1, n = 1000, plot = T)
simulateResiduals(m2, n = 1000, plot = T)
simulateResiduals(m3, n = 1000, plot = T)
simulateResiduals(m4, n = 1000, plot = T)

###### Compare via AICc
aictab(list(m1, m2))
aictab(list(m3, m4))
aictab(list(m1, m2, m3, m4))

###### Repeat for other FixedLocationIDs
dater392 <- dater %>% filter(FixedLocationID=="392")
m5 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~1, family = nbinom2, data = dater392)
m6 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~0, family = nbinom2, data = dater392)
m7 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~1, family = poisson, data = dater392)
m8 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~0, family = poisson, data = dater392)
simulateResiduals(m5, n = 1000, plot = T)
simulateResiduals(m6, n = 1000, plot = T)
simulateResiduals(m7, n = 1000, plot = T)
simulateResiduals(m8, n = 1000, plot = T)
aictab(list(m5, m6))
aictab(list(m7, m8))
aictab(list(m5, m6, m7, m8))

dater393 <- dater %>% filter(FixedLocationID=="393")
m9 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~1, family = nbinom2, data = dater393)
m10 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~0, family = nbinom2, data = dater393)
m11 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~1, family = poisson, data = dater393)
m12 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~0, family = poisson, data = dater393)
simulateResiduals(m9, n = 1000, plot = T)
simulateResiduals(m10, n = 1000, plot = T)
simulateResiduals(m11, n = 1000, plot = T)
simulateResiduals(m12, n = 1000, plot = T)
aictab(list(m9, m10))
aictab(list(m11, m12))
aictab(list(m9, m10, m11, m12))

dater394 <- dater %>% filter(FixedLocationID=="394")
m13 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~1, family = nbinom2, data = dater394)
m14 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~0, family = nbinom2, data = dater394)
m15 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~1, family = poisson, data = dater394)
m16 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~0, family = poisson, data = dater394)
simulateResiduals(m13, n = 1000, plot = T)
simulateResiduals(m14, n = 1000, plot = T)
simulateResiduals(m15, n = 1000, plot = T)
simulateResiduals(m16, n = 1000, plot = T)
aictab(list(m13, m14))
aictab(list(m15, m16))
aictab(list(m13, m14, m15, m16))

dater395 <- dater %>% filter(FixedLocationID=="395")
m17 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~1, family = nbinom2, data = dater395)
m18 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~0, family = nbinom2, data = dater395)
m19 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~1, family = poisson, data = dater395)
m20 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~0, family = poisson, data = dater395)
simulateResiduals(m17, n = 1000, plot = T)
simulateResiduals(m18, n = 1000, plot = T)
simulateResiduals(m19, n = 1000, plot = T)
simulateResiduals(m20, n = 1000, plot = T)
aictab(list(m17, m18))
aictab(list(m19, m20))
aictab(list(m17, m18, m19, m20))

dater396 <- dater %>% filter(FixedLocationID=="396")
m21 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~1, family = nbinom2, data = dater396)
m22 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~0, family = nbinom2, data = dater396)
m23 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~1, family = poisson, data = dater396)
m24 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~0, family = poisson, data = dater396)
simulateResiduals(m21, n = 1000, plot = T)
simulateResiduals(m22, n = 1000, plot = T)
simulateResiduals(m23, n = 1000, plot = T)
simulateResiduals(m24, n = 1000, plot = T)
aictab(list(m21, m22))
aictab(list(m23, m24))
aictab(list(m21, m22, m23, m24))

dater397 <- dater %>% filter(FixedLocationID=="397")
m25 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~1, family = nbinom2, data = dater397)
m26 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~0, family = nbinom2, data = dater397)
m27 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~1, family = poisson, data = dater397)
m28 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~0, family = poisson, data = dater397)
simulateResiduals(m25, n = 1000, plot = T)
simulateResiduals(m26, n = 1000, plot = T)
simulateResiduals(m27, n = 1000, plot = T)
simulateResiduals(m28, n = 1000, plot = T)
aictab(list(m25, m26))
aictab(list(m27, m28))
aictab(list(m25, m26, m27, m28))

dater398 <- dater %>% filter(FixedLocationID=="398")
m29 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~1, family = nbinom2, data = dater398)
m30 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~0, family = nbinom2, data = dater398)
m31 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~1, family = poisson, data = dater398)
m32 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~0, family = poisson, data = dater398)
simulateResiduals(m29, n = 1000, plot = T)
simulateResiduals(m30, n = 1000, plot = T)
simulateResiduals(m31, n = 1000, plot = T)
simulateResiduals(m32, n = 1000, plot = T)
aictab(list(m29, m30))
aictab(list(m31, m32))
aictab(list(m29, m30, m31, m32))

dater399 <- dater %>% filter(FixedLocationID=="399")
m33 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~1, family = nbinom2, data = dater399)
m34 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~0, family = nbinom2, data = dater399)
m35 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~1, family = poisson, data = dater399)
m36 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~0, family = poisson, data = dater399)
simulateResiduals(m33, n = 1000, plot = T)
simulateResiduals(m34, n = 1000, plot = T)
simulateResiduals(m35, n = 1000, plot = T)
simulateResiduals(m36, n = 1000, plot = T)
aictab(list(m33, m34))
aictab(list(m35, m36))
aictab(list(m33, m34, m35, m36))

dater400 <- dater %>% filter(FixedLocationID=="400")
m37 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~1, family = nbinom2, data = dater400)
m38 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~0, family = nbinom2, data = dater400)
m39 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~1, family = poisson, data = dater400)
m40 <- glmmTMB(NumLegal ~ 1, dispformula = ~1, ziformula = ~0, family = poisson, data = dater400)
simulateResiduals(m37, n = 1000, plot = T)
simulateResiduals(m38, n = 1000, plot = T)
simulateResiduals(m39, n = 1000, plot = T)
simulateResiduals(m40, n = 1000, plot = T)
aictab(list(m37, m38))
aictab(list(m39, m40))
aictab(list(m37, m38, m39, m40))

###### Final list of best-fitting model for each FixedLocationID
## FOR 2026 WE MIGHT NEED TO MODIFY THESE DEPENDING ON UPDATED MODEL SELECTION RESULTS
## ONLY 1 SITE (FixedLocationID=="395" = m19) WAS A ZERO INFLATED (POISSON) MODEL; THE REST WERE REGULAR OLD NB MODELS
modList <- list(m2, m6, m10, m14, m19, m22, m26, m30, m34, m38)

###### Create an empty data frame for storing simulation results
simDater <- data.frame(FixedLocationID = sort(unique(dater$FixedLocationID)), nSamples = NA, meanLog = NA, meanLogSE = NA, ziLogit = NA, ziLogitSE = NA, disp = NA, nbinom = NA, zi = NA) %>% mutate(nSamples = if_else(FixedLocationID %in% c("391", "397"), 50, 30))

###### Populate simDate with parameter values for each model (to be used in the simulation)
for (i in 1:length(modList)){
  meanLog <- unique(predict(modList[[i]], type = "link", se.fit = T)$fit)
  meanLogSE <- unique(predict(modList[[i]], type = "link", se.fit = T)$se.fit)
  ziLogit <- unique(predict(modList[[i]], type = "zlink", se.fit = T)$fit)
  ziLogitSE <- unique(predict(modList[[i]], type = "zlink", se.fit = T)$se.fit)
  disp <- if_else(is.na(sigma(modList[[i]]))==TRUE, unique(predict(modList[[i]], type = "disp", se.fit = T)$fit), sigma(modList[[i]]))
  nbinom <- if_else("nbinom2" %in% family(modList[[i]]), 1, 0)
  zi <- if_else(modList[[i]]$modelInfo$allForm$ziformula=="~Year"|modList[[i]]$modelInfo$allForm$ziformula=="~1", 1, 0)
  simDater[i,3:ncol(simDater)] <- c(meanLog, meanLogSE, ziLogit, ziLogitSE, disp, nbinom, zi)
}

###### Look at the file
datatable(simDater)

###### Desired number of simulation replicates per FixedLocationID
nSims <- 1000

###### Expand simDater base file by the number of simulation replicates (each row is replicated nSims times)
simDater <- simDater[rep(seq_len(nrow(simDater)), each = nSims), ]

###### Number the simulation replicates for each FixedLocationID
simDater$simRep <- rep(1:nSims, times = length(unique(simDater$FixedLocationID)))

###### Create a file for storing results (one row per simRep x FixedLocationID)
resultsFile <- simDater %>% select(simRep, FixedLocationID, nSamples, meanLog, meanLogSE, ziLogit, ziLogitSE, disp, nbinom, zi) %>% mutate(trueMean = meanLog, estMean = NA, lwr = NA, upr = NA, bias = NA)

###### Optional: set seed for reproducibility
set.seed(42)
for (i in 1:nrow(resultsFile)){
  # Build simDat on the fly: one row per quadrat for this simRep x FixedLocationID
  n_samp <- resultsFile$nSamples[i]
  simDat <- data.frame(
    meanLog   = rep(resultsFile$meanLog[i],   n_samp),
    meanLogSE = rep(resultsFile$meanLogSE[i], n_samp),
    ziLogit   = rep(resultsFile$ziLogit[i],   n_samp),
    ziLogitSE = rep(resultsFile$ziLogitSE[i], n_samp),
    disp      = rep(resultsFile$disp[i],      n_samp),
    nbinom    = rep(resultsFile$nbinom[i],    n_samp),
    zi        = rep(resultsFile$zi[i],        n_samp),
    count     = NA_real_
  )
  fit_result <- NULL
  
  # Draw system-level parameter values ONCE for this replicate
  sim_mu <- exp(rnorm(1, mean = simDat$meanLog[1], sd = simDat$meanLogSE[1]))
  sim_disp <- simDat$disp[1]
  
  # Only draw zi when the model uses it
  if (simDat$zi[1] > 0) {
    sim_zi <- plogis(rnorm(1, simDat$ziLogit[1], simDat$ziLogitSE[1]))
  } else {
    sim_zi <- 0
  }
  
  if (simDat$nbinom[1] == 1){
    if(simDat$zi[1] == 0){
      simDat$count <- rnbinom(nrow(simDat), mu = sim_mu, size = sim_disp)
      fit_result <- try(sampFit <- glmmTMB(count ~ 1, family = nbinom2, ziformula = ~0, data = simDat), silent = TRUE)
    }
    if(simDat$zi[1] > 0){
      simDat$count <- ifelse(rbinom(nrow(simDat), 1, sim_zi) > 0, 0,
                             rnbinom(nrow(simDat), mu = sim_mu, size = sim_disp))
      fit_result <- try(sampFit <- glmmTMB(count ~ 1, family = nbinom2, ziformula = ~1, data = simDat), silent = TRUE)
    }
  }
  if (simDat$nbinom[1] == 0){
    if(simDat$zi[1] == 0){
      simDat$count <- rpois(nrow(simDat), lambda = sim_mu)
      fit_result <- try(sampFit <- glmmTMB(count ~ 1, family = poisson, ziformula = ~0, data = simDat), silent = TRUE)
    }
    if(simDat$zi[1] > 0){
      simDat$count <- ifelse(rbinom(nrow(simDat), 1, sim_zi) > 0, 0,
                             rpois(nrow(simDat), lambda = sim_mu))
      fit_result <- try(sampFit <- glmmTMB(count ~ 1, family = poisson, ziformula = ~1, data = simDat), silent = TRUE)
    }
  }
  
  # Guard against failed model fits
  if (is.null(fit_result) || inherits(fit_result, "try-error")) {
    resultsFile$estMean[i] <- NA
    resultsFile$lwr[i] <- NA
    resultsFile$upr[i] <- NA
    resultsFile$bias[i] <- NA
    next
  }
  
  resultsFile$estMean[i] <- fixef(sampFit)$cond[1]
  cl_result <- try(CLs <- confint(sampFit, level = 0.95)[1,], silent = TRUE)
  if (inherits(cl_result, "try-error")) {
    resultsFile$lwr[i] <- NA
    resultsFile$upr[i] <- NA
  } else {
    resultsFile$lwr[i] <- CLs[1]
    resultsFile$upr[i] <- CLs[2]
  }
  resultsFile$bias[i] <- exp(fixef(sampFit)$cond[1]) - exp(resultsFile$trueMean[i])
  if (i %% 100 == 0) cat("Iteration", i, "of", nrow(resultsFile), "\n")
}

###### Summarize the results
EstimatedMeans <- resultsFile %>% group_by(FixedLocationID) %>% summarise(true_mean = mean(trueMean), est_mean = mean(estMean, na.rm = T), RBias = mean(bias, na.rm = T)/mean(exp(trueMean), na.rm = T)*100, RRMSE = sqrt(mean(bias^2, na.rm = T))/mean(exp(trueMean), na.rm = T)*100)

###### Pivot EstimatedMeans and filter for RRMSE
RRMSE_results <- EstimatedMeans %>% pivot_longer(cols = c(RBias, RRMSE), names_to = "Metric", values_to = "Value") %>% filter(Metric =="RRMSE") %>% mutate(Metric = factor(Metric, levels = c("RRMSE"), labels = c("Relative Root Mean Squared Error (% of true mean)")))

###### Plot RRMSE
(RRMSE <- ggplot(RRMSE_results, aes(x = FixedLocationID, y = Value, group = FixedLocationID)) + geom_bar(stat = "identity", color = "black") + scale_y_continuous(limits = c(0,70), breaks = seq(0,70,10), expand = expansion(add = c(0,0))) + theme_bw() + theme(panel.grid.major.y = element_line(color = "grey90", linetype = "solid"), panel.grid.minor.y = element_line(color = "grey90", linetype = "dashed"), axis.text = element_text(color = "black"), legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + labs(x = NULL, y = "RRMSE (%)") + facet_grid(~Metric))

###### Pivot EstimatedMeans and filter for RBIAS
RBIAS_results <- EstimatedMeans %>% pivot_longer(cols = c(RBias, RRMSE), names_to = "Metric", values_to = "Value") %>% filter(Metric =="RBias") %>% mutate(Metric = factor(Metric, levels = c("RBias"), labels = c("Relative Bias (% of true mean)")))

###### Plot RBIAS
(RBIAS <- ggplot(RBIAS_results, aes(x = FixedLocationID, y = Value, group = FixedLocationID)) + geom_bar(stat = "identity", color = "black") + scale_y_continuous(limits = c(-1,10), breaks = seq(-1,10,1), expand = expansion(add = c(0,0))) + theme_bw() + theme(panel.grid.major.y = element_line(color = "grey90", linetype = "solid"), panel.grid.minor.y = element_line(color = "grey90", linetype = "dashed"), axis.text = element_text(color = "black"), legend.position = "bottom", axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) + labs(x = NULL, y = "RBIAS (%)") + facet_grid(~Metric))

###### Import FixedLocationID areas and Station Names for extrapolation
areas <- data.frame(FixedLocationID = sort(unique(dater$FixedLocationID)), StationName = c("RESTORE Cat Point", "RESTORE Monkey's Elbow", "RESTORE Peanut Ridge", "RESTORE Cat Point Spur", "RESTORE Platform", "RESTORE East Bulkhead", "RESTORE Easthole", "SBM Lighthouse-Restoration", "SBM East Lumps-Restoration", "SBM Cat Point-Restoration"), acres = c(50.25, 27.15, 20.53, 12.45, 22.00, 23.95, 42.89, 10.66, 9.38, 18.25)) %>% mutate(m2 = acres*4046.86, qm2 = round(m2*4,0))

###### Look at areas and station names
datatable(areas)

###### Empty data frame for setting up extrapolation bia parametric bootstrapping
extrapDat <- data.frame(areas, meanLog = NA, meanLogSE = NA, ziLogit = NA, ziLogitSE = NA, disp = NA, nbinom = NA, zi = NA)

###### Loop through models and populate extrapDat with relevant parameter values
for (i in 1:length(modList)){
  meanLog <- unique(predict(modList[[i]], type = "link", se.fit = T)$fit)
  meanLogSE <- unique(predict(modList[[i]], type = "link", se.fit = T)$se.fit)
  ziLogit <- unique(predict(modList[[i]], type = "zlink", se.fit = T)$fit)
  ziLogitSE <- unique(predict(modList[[i]], type = "zlink", se.fit = T)$se.fit)
  disp <- if_else(is.na(sigma(modList[[i]]))==TRUE, unique(predict(modList[[i]], type = "disp", se.fit = T)$fit), sigma(modList[[i]]))
  nbinom <- if_else("nbinom2" %in% family(modList[[i]]), 1, 0)
  zi <- if_else(modList[[i]]$modelInfo$allForm$ziformula=="~Year"|modList[[i]]$modelInfo$allForm$ziformula=="~1", 1, 0)
  extrapDat[i,6:ncol(extrapDat)] <- c(meanLog, meanLogSE, ziLogit, ziLogitSE, disp, nbinom, zi)
}

###### Look at population extrapDat
datatable(extrapDat)

###### Create another data frame for storing bootstrapped extrapolated counts
extrapDat2 <- extrapDat %>% mutate(`Mean Total` = NA, `L95 Total` = NA, `U95 Total` = NA, `3% Total` = NA, `5% Total` = NA, `3% L95` = NA, `5% L95` = NA)

####### Extrapolation using scaled distributions instead of summing individual quadrat draws
####### Sum of n iid Poisson(mu) = Poisson(n*mu); Sum of n iid NB(mu, theta) = NB(n*mu, n*theta)
####### For ZI models: draw n_active ~ Binom(n_quad, 1 - zi_prob), then scale by n_active
for (j in 1:nrow(extrapDat2)){
  totalCount <- numeric(nSims)
  n_quad <- extrapDat2$qm2[j]
  for (i in 1:nSims){
    # Draw system-level parameter values ONCE for this replicate
    sim_mu  <- exp(rnorm(1, mean = extrapDat2$meanLog[j], sd = extrapDat2$meanLogSE[j]))
    sim_zi  <- plogis(rnorm(1, extrapDat2$ziLogit[j], extrapDat2$ziLogitSE[j]))
    
    # Determine number of active (non-structural-zero) quadrats
    if (extrapDat2$zi[j] == 1) {
      n_active <- rbinom(1, n_quad, 1 - sim_zi)
    } else {
      n_active <- n_quad
    }
    # Draw total count from scaled distribution
    if (n_active == 0) {
      totalCount[i] <- 0
    } else if (extrapDat2$nbinom[j] == 1) {
      totalCount[i] <- rnbinom(1, mu = n_active * sim_mu, size = n_active * extrapDat2$disp[j])
    } else {
      totalCount[i] <- rpois(1, lambda = n_active * sim_mu)
    }
  }
  extrapDat2$`Mean Total`[j] <- mean(totalCount)
  extrapDat2$`L95 Total`[j]  <- quantile(totalCount, 0.025)
  extrapDat2$`U95 Total`[j]  <- quantile(totalCount, 0.975)
  extrapDat2$`3% Total`[j]   <- 0.03 * mean(totalCount)
  extrapDat2$`5% Total`[j]   <- 0.05 * mean(totalCount)
  extrapDat2$`3% L95`[j]     <- 0.03 * quantile(totalCount, 0.025)
  extrapDat2$`5% L95`[j]     <- 0.05 * quantile(totalCount, 0.025)
}

###### Simplify results
extrapFIN <- extrapDat2 %>% dplyr::select(-c(qm2, meanLog, meanLogSE, ziLogit, ziLogitSE, disp, nbinom, zi)) %>% mutate_at(4:11, round, 0)

warnings()

###### Extrapolated total number of oysters based on parametric bootstrap
extrapOYSTERS_EST <- extrapFIN %>% select(FixedLocationID, StationName, acres, m2, `Mean Total`, `L95 Total`, `U95 Total`) %>% mutate(MeanTotal = `Mean Total`, L95 = `L95 Total`, U95 = `U95 Total`) %>% select(FixedLocationID, StationName, acres, m2, MeanTotal, L95, U95) %>% mutate(`3% Total` = 0.03*MeanTotal, `5% Total` = 0.05*MeanTotal, `3% LCL` = 0.03*L95, `5% LCL` = 0.05*L95) %>% mutate_at(4:11, round, 0)

###### Extrapolated total number of oysters based on raw data
extrapOYSTERS_RAW <- datSumm %>% select(FixedLocationID, MeanLegal, lwr, upr) %>% left_join(areas, by = "FixedLocationID") %>% mutate(MeanTotal = MeanLegal*qm2, L95 = lwr*qm2, U95 = upr*qm2) %>% select(FixedLocationID, StationName, acres, m2, MeanTotal, L95, U95) %>% mutate(L95 = 0, U95 = 0) %>% mutate(`3% Total` = 0.03*MeanTotal, `5% Total` = 0.05*MeanTotal, `3% LCL` = 0.03*L95, `5% LCL` = 0.05*L95) %>% mutate_at(4:11, round, 0)

###### Extrapolated total number of oysters based on parametric bootstrap
datatable(extrapOYSTERS_EST)

###### Extrapolated total number of oysters based on raw data
datatable(extrapOYSTERS_RAW)

###### Extrapolated total number of bags (100 oysters per box) based on parametric bootstrap
extrapBOXES_EST <- extrapFIN %>% select(FixedLocationID, StationName, acres, m2, `Mean Total`, `L95 Total`, `U95 Total`) %>% mutate(MeanTotal = `Mean Total`/100, L95 = `L95 Total`/100, U95 = `U95 Total`/100) %>% select(FixedLocationID, StationName, acres, m2, MeanTotal, L95, U95) %>% mutate(`3% Total` = 0.03*MeanTotal, `5% Total` = 0.05*MeanTotal, `3% LCL` = 0.03*L95, `5% LCL` = 0.05*L95) %>% mutate_at(4:11, round, 0)

###### Extrapolated total number of bags (100 oysters per box) based on raw data
extrapBOXES_RAW <- datSumm %>% select(FixedLocationID, MeanLegal, lwr, upr) %>% left_join(areas, by = "FixedLocationID") %>% mutate(MeanTotal = MeanLegal*qm2/100, L95 = 0, U95 = 0) %>% select(FixedLocationID, StationName, acres, m2, MeanTotal, L95, U95) %>% mutate(`3% Total` = 0.03*MeanTotal, `5% Total` = 0.05*MeanTotal, `3% LCL` = 0.03*L95, `5% LCL` = 0.05*L95) %>% mutate_at(4:11, round, 0)

###### Extrapolated total number of bags (100 oysters per box) based on parametric bootstrap
datatable(extrapBOXES_EST)

###### Extrapolated total number of bags (100 oysters per box) based on raw data
datatable(extrapBOXES_RAW)

###### Extrapolated total number of bags (225 oysters per bag) based on parametric bootstrap
extrapBAGS_EST <- extrapFIN %>% select(FixedLocationID, StationName, acres, m2, `Mean Total`, `L95 Total`, `U95 Total`) %>% mutate(MeanTotal = `Mean Total`/225, L95 = `L95 Total`/225, U95 = `U95 Total`/225) %>% select(FixedLocationID, StationName, acres, m2, MeanTotal, L95, U95) %>% mutate(`3% Total` = 0.03*MeanTotal, `5% Total` = 0.05*MeanTotal, `3% LCL` = 0.03*L95, `5% LCL` = 0.05*L95) %>% mutate_at(4:11, round, 0)

###### Extrapolated total number of bags (225 oysters per bag) based on raw data
extrapBAGS_RAW <- datSumm %>% select(FixedLocationID, MeanLegal, lwr, upr) %>% left_join(areas, by = "FixedLocationID") %>% mutate(MeanTotal = MeanLegal*qm2/225, L95 = 0, U95 = 0) %>% select(FixedLocationID, StationName, acres, m2, MeanTotal, L95, U95) %>% mutate(`3% Total` = 0.03*MeanTotal, `5% Total` = 0.05*MeanTotal, `3% LCL` = 0.03*L95, `5% LCL` = 0.05*L95) %>% mutate_at(4:11, round, 0)

###### Extrapolated total number of bags (225 oysters per bag) based on parametric bootstrap
datatable(extrapBAGS_EST)

###### Extrapolated total number of bags (225 oysters per bag) based on raw data
datatable(extrapBAGS_RAW)

###### Combine everything
extras <- rbind(extrapOYSTERS_EST %>% mutate(group = "Model", variable = "Oysters"),
                extrapOYSTERS_RAW %>% mutate(group = "Raw", variable = "Oysters"),
                extrapBOXES_EST %>% mutate(group = "Model", variable = "Boxes"),
                extrapBOXES_RAW %>% mutate(group = "Raw", variable = "Boxes"),
                extrapBAGS_EST %>% mutate(group = "Model", variable = "Bags"),
                extrapBAGS_RAW %>% mutate(group = "Raw", variable = "Bags"))

###### Plot total number of oysters
ggplot(extras %>% filter(variable == "Oysters") %>% mutate(MeanTotal = MeanTotal/1e6, L95 = L95/1e6, U95 = U95/1e6), aes(x = FixedLocationID, y = MeanTotal, fill = group)) + geom_bar(stat = "identity", position = position_dodge(width = 1, preserve = "single"), color = "black") + geom_errorbar(aes(ymin = L95, ymax = U95), width = 0.5, position = position_dodge(width = 1, preserve = "single")) + theme_bw() + facet_grid(~variable) + scale_fill_brewer(palette = "Dark2") + theme(axis.text = element_text(colour = "black"), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color = "grey90", linetype = "solid"), panel.grid.minor = element_line(color = "grey90", linetype = "dashed"), legend.title = element_blank(), legend.position = "bottom") + scale_y_continuous(limits = c(0,8), breaks = seq(0,8,0.5), expand = expansion(add = c(0,0))) + labs(x = "Fixed Location ID", y = "Total oysters (Mean ± 95% CLs in millions)")

###### Plot total number of boxes
ggplot(extras %>% filter(variable == "Boxes"), aes(x = FixedLocationID, y = MeanTotal, fill = group)) + geom_bar(stat = "identity", position = position_dodge(width = 1, preserve = "single"), color = "black") + geom_errorbar(aes(ymin = L95, ymax = U95), width = 0.5, position = position_dodge(width = 1, preserve = "single")) + theme_bw() + facet_grid(~variable) + scale_fill_brewer(palette = "Dark2") + theme(axis.text = element_text(colour = "black"), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color = "grey90", linetype = "solid"), panel.grid.minor = element_line(color = "grey90", linetype = "dashed"), legend.title = element_blank(), legend.position = "bottom") + scale_y_continuous(limits = c(0,75000), breaks = seq(0,75000, 5000), expand = expansion(add = c(0,0))) + labs(x = "Fixed Location ID", y = "Total boxes (Mean ± 95% CLs; 1 box = 100 oysters)")

###### Plot total number of bags
ggplot(extras %>% filter(variable == "Bags"), aes(x = FixedLocationID, y = MeanTotal, fill = group)) + geom_bar(stat = "identity", position = position_dodge(width = 1, preserve = "single"), color = "black") + geom_errorbar(aes(ymin = L95, ymax = U95), width = 0.5, position = position_dodge(width = 1, preserve = "single")) + theme_bw() + facet_grid(~variable) + scale_fill_brewer(palette = "Dark2") + theme(axis.text = element_text(colour = "black"), panel.grid.major.x = element_blank(), panel.grid.major.y = element_line(color = "grey90", linetype = "solid"), panel.grid.minor = element_line(color = "grey90", linetype = "dashed"), legend.title = element_blank(), legend.position = "bottom") + scale_y_continuous(limits = c(0,35000), breaks = seq(0,35000, 2500), expand = expansion(add = c(0,0))) + labs(x = "Fixed Location ID", y = "Total bags (Mean ± 95% CLs; 1 bag = 225 oysters)")