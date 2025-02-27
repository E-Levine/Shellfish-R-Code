### Load required packages
library(tidyverse)
library(lubridate)

### Read in data: is an observation is -999 we assign an NA value
SH <- read.csv("spatSH2.csv", na.strings = c("-999"))

### Do some data formatting: here, convert date to date format, extract year, month, and day, convert year to a factor, and then omit rows that contain NA values
SH <- SH %>% mutate(Date = as.Date(Date, format = "%m/%d/%Y"), Year = year(Date), Month = month(Date), Day = day(Date), YearF = factor(Year)) %>% filter(complete.cases(.))

### Create a new variable called GRP which represents unique combinations of Date, Site, Station, and Quadrat (so unique quadrats)
SH$GRP <- droplevels(interaction(SH$Date, SH$Site, SH$Station, SH$Quadrat))

### Assign Quarters and then convert to a factor
SH$Quarter <- "First"
SH$Quarter <- ifelse(SH$Month>=4 & SH$Month<=6, "Second", SH$Quarter)
SH$Quarter <- ifelse(SH$Month>=7 & SH$Month<=9, "Third", SH$Quarter)
SH$Quarter <- ifelse(SH$Month>=10 & SH$Month<=12, "Fourth", SH$Quarter)
SH$Quarter <- factor(SH$Quarter, levels = c("First", "Second", "Third", "Fourth"))

# Convert Stations to a factor
SH$Station <- factor(SH$Station)
# Convert Survey to a factor
SH$Survey <- factor(SH$Survey)

### Look at data: sort of a weird distribution. Not easy to fit a simple linear regression here, and a count model (e.g,. Poisson or negative binomial) doesn't make sense. Sometimes a Gamma regression works, but in this case that didn't pan out either. 
hist(SH$SH[SH$YearF=="2015"], breaks = 100)
hist(SH$SH[SH$YearF=="2016"], breaks = 100)
hist(SH$SH[SH$YearF=="2017"], breaks = 100)
hist(SH$SH[SH$YearF=="2018"], breaks = 100)
hist(SH$SH[SH$YearF=="2019"], breaks = 100)

### So, we "resort" to nonparametric bootstrapping. I put "resort" in quotes because this is a perfectly reasonable approach. The procedure is as follows: 

#################################
# BY SURVEY AND STATION N ==100 #
#################################
# Take the original data set, SH, group it by GRP (unique quadrat IDs), determine the number of samples (oysters) measure in each GRP (N), filter for only those with 100 oysters, and then calculate the true mean SH for each Survey
SHResample <- SH %>% group_by(GRP) %>% mutate(N = 1:n(), maxN = max(N)) %>% filter(maxN==100) %>% droplevels(.) %>% ungroup() %>% group_by(Survey, Station) %>% mutate(mnTrue = mean(SH))
# Then set up the simulations with objects we'll need to store results and various simulation conditions
SHSSMN <- NULL
SHSSMNREPS <- NULL
maxSamp <- seq(5, 100, 5)
nSims <- 1000
quadrats <- unique(SHResample$GRP)
# Now start the simulation, looping through maxSamp vector and simulation replicates
for(i in 1:length(maxSamp)){
  for (j in 1:nSims){
    SHTemp <- SHResample %>% group_by(GRP) %>% sample_n(size = maxSamp[i], replace = T)
    SHSSMN <- suppressMessages(SHTemp %>% group_by(Survey, Station, mnTrue) %>% summarise(mnSH = mean(SH)) %>% mutate(bias = mnSH - mnTrue, nSamples = maxSamp[i], replicate = j))
    SHSSMNREPS <- rbind(SHSSMNREPS, SHSSMN)
    if(j %% 1000==0){cat(paste0("iteration: ", j, "; nSamples: ", maxSamp[i], "\n"))}
  }
} 
# Now summarize the full results, calculating Relative Bias (expressed as % of mean) and Relative Root Mean Squared Error (expressed as a % of mean)
finSS100 <- SHSSMNREPS %>% group_by(Survey, Station, nSamples) %>% summarise(RBIAS = mean(bias)/mean(mnTrue)*100, RRMSE = sqrt(mean(bias^2))/mean(mnTrue)*100)
# Then make a plot for RRMSE: RRMSE is equivalent to a CV and values <20% are considered good
SurveysRRMSE100 <- ggplot(finSS100, aes(x = nSamples, y = RRMSE, group = Station, color = Station)) + geom_point(stat = "identity") + geom_line() + facet_grid(~Survey) + theme(legend.position = "bottom") + theme_bw() + theme(legend.title = element_blank(), legend.position = "bottom", panel.grid = element_blank()) + scale_y_continuous(limits = c(0,35), breaks = seq(0,35,5)) + scale_x_continuous(limits = c(0,100), breaks = seq(0,100,10)) + labs(x = "Number of individuals per quadrat", y = "Relative Root Mean Squared Error (% of true mean)") + theme(axis.text = element_text(color = "black")) + geom_hline(yintercept = 20, linetype = "dashed") + geom_vline(xintercept = 50, linetype = "dashed")
ggsave("SurveysRRMSE_N100.png", SurveysRRMSE100, dpi = 600, height = 5, width = 20)

# And one for RBIAS: generally these are unbiased; but these are also expressed as a percentage
SurveysRBIAS100 <- ggplot(finSS100, aes(x = nSamples, y = RBIAS, group = Station, color = Station)) + geom_point(stat = "identity") + geom_line() + facet_grid(~Survey) + theme(legend.position = "bottom") + theme_bw() + theme(legend.title = element_blank(), legend.position = "bottom", panel.grid = element_blank()) + scale_y_continuous(limits = c(-1.8,1.8), breaks = seq(-1.8,1.8,0.2)) + scale_x_continuous(limits = c(0,100), breaks = seq(0,100,10)) + labs(x = "Number of individuals per quadrat", y = "Relative Bias (% of true mean)") + theme(axis.text = element_text(color = "black")) + geom_hline(yintercept = 0, linetype = "solid") + geom_hline(yintercept = 0.2, linetype = "dashed") + geom_hline(yintercept = -0.2, linetype = "dashed") + geom_vline(xintercept = 50, linetype = "dashed")
ggsave("SurveysRBIAS_N100.png", SurveysRBIAS100, dpi = 600, height = 8, width = 20)

##################################
# BY QUARTER AND STATION N ==100 #
##################################
# Take the original data set, SH, group it by GRP (unique quadrat IDs), determine the number of samples (oysters) measure in each GRP (N), filter for only those with 100 oysters, and then calculate the true mean SH for each Quarter
SHResample <- SH %>% group_by(GRP) %>% mutate(N = 1:n(), maxN = max(N)) %>% filter(maxN==100) %>% droplevels(.) %>% ungroup() %>% group_by(Quarter, Station) %>% mutate(mnTrue = mean(SH))
# Then set up the simulations with objects we'll need to store results and various simulation conditions
SHQSMN <- NULL
SHQSMNREPS <- NULL
maxSamp <- seq(5, 100, 5)
nSims <- 1000
quadrats <- unique(SHResample$GRP)
# Now start the simulation, looping through maxSamp vector and simulation replicates
for(i in 1:length(maxSamp)){
  for (j in 1:nSims){
    SHTemp <- SHResample %>% group_by(GRP) %>% sample_n(size = maxSamp[i], replace = T)
    SHQSMN <- suppressMessages(SHTemp %>% group_by(Quarter, Station, mnTrue) %>% summarise(mnSH = mean(SH)) %>% mutate(bias = mnSH - mnTrue, nSamples = maxSamp[i], replicate = j))
    SHQSMNREPS <- rbind(SHQSMNREPS, SHQSMN)
    if(j %% 1000==0){cat(paste0("iteration: ", j, "; nSamples: ", maxSamp[i], "\n"))}
  }
}
# Now summarize the full results, calculating Relative Bias (expressed as % of mean) and Relative Root Mean Squared Error (expressed as a % of mean)
finQS100 <- SHQSMNREPS %>% group_by(Quarter, Station, nSamples) %>% summarise(RBIAS = mean(bias)/mean(mnTrue)*100, RRMSE = sqrt(mean(bias^2))/mean(mnTrue)*100)
# Then make a plot for RRMSE: RRMSE is equivalent to a CV and values <20% are considered good
StationsRRMSE100 <- ggplot(finQS100, aes(x = nSamples, y = RRMSE, group = Station, color = Station)) + geom_point(stat = "identity") + geom_line() + facet_grid(~Quarter) + theme(legend.position = "bottom") + theme_bw() + theme(legend.title = element_blank(), legend.position = "bottom", panel.grid = element_blank()) + scale_y_continuous(limits = c(0,35), breaks = seq(0,35,5)) + scale_x_continuous(limits = c(0,100), breaks = seq(0,100,10)) + labs(x = "Number of individuals per quadrat", y = "Relative Root Mean Squared Error (% of true mean)") + theme(axis.text = element_text(color = "black")) + geom_hline(yintercept = 20, linetype = "dashed") + geom_vline(xintercept = 50, linetype = "dashed")
ggsave("StationRRMSE_N100.png", StationsRRMSE100, dpi = 600, height = 5, width = 10)
# And one for RBIAS: generally these are unbiased; but these are also expressed as a percentage
StationsRBIAS100 <- ggplot(finQS100, aes(x = nSamples, y = RBIAS, group = Station, color = Station)) + geom_point(stat = "identity") + geom_line() + facet_grid(~Quarter) + theme(legend.position = "bottom") + theme_bw() + theme(legend.title = element_blank(), legend.position = "bottom", panel.grid = element_blank()) + scale_y_continuous(limits = c(-1.0,1.2), breaks = seq(-1.0,1.2,0.2)) + scale_x_continuous(limits = c(0,100), breaks = seq(0,100,10)) + labs(x = "Number of individuals per quadrat", y = "Relative Bias (% of true mean)") + theme(axis.text = element_text(color = "black")) + geom_hline(yintercept = 0, linetype = "solid") + geom_hline(yintercept = 0.2, linetype = "dashed") + geom_hline(yintercept = -0.2, linetype = "dashed") + geom_vline(xintercept = 50, linetype = "dashed")
ggsave("StationRBIAS_N100.png", StationsRBIAS100, dpi = 600, height = 5, width = 10)


########################
### BY YEAR ~ N==100 ###
########################
### Take the original data set, SH, group it by GRP (unique quadrat IDs), determine the number of samples (oysters) measure in each GRP (N), filter for only those with 100 oysters, and then calculate the true mean SH for each Year
SHResample <- SH %>% group_by(GRP) %>% mutate(N = 1:n(), maxN = max(N)) %>% filter(maxN==100) %>% droplevels(.) %>% ungroup() %>% group_by(YearF) %>% mutate(mnTrue = mean(SH))

### Then set up the simulations with objects we'll need to store results and various simulation conditions
SHYMN <- NULL # empty object for storing results
SHYMNREPS <- NULL # empty object for storing results
maxSamp <- seq(5, 100, 5) # vector of sample sizes: number of individuals measured per quadrat
nSims <- 1000 # number of simulation reps for each maxSamp (I'd jack this up to 5000 or even 10000)
### Now start the simulation, looping through maxSamp vector and simulation replicates
for(i in 1:length(maxSamp)){
  for (j in 1:nSims){
### First Select a random sample of maxSamp[i] indviduals from each quadrat, where maxSamp ranges from 5 to 100 in increments of 5. 
### Here, we're Sampling with replacement = T, so the same observations in a quadrat can appear more than once.       
    SHTemp <- SHResample %>% group_by(GRP) %>% sample_n(size = maxSamp[i], replace = T)
### Now group by Year and the true mean, calculate the mean SH of the simulated data, and then calculate the bias; also record nSample and the replicate number    
    SHYMN <- suppressMessages(SHTemp %>% group_by(YearF, mnTrue) %>% summarise(mnSH = mean(SH)) %>% mutate(bias = mnSH - mnTrue, nSamples = maxSamp[i], replicate = j))
### Then bind simulation replicate results to previous ones to make a full results data frame
    SHYMNREPS <- rbind(SHYMNREPS, SHYMN)
### This is just a counter to keep track of progress; change 1000 to something else (e.g. 100) if you want more regular updates        
    if(j %% 1000==0){cat(paste0("iteration: ", j, "; nSamples: ", maxSamp[i], "\n"))}
  }
}
### Now summarize the full results, calculating Relative Bias (expressed as % of mean) and Relative Root Mean Squared Error (expressed as a % of mean)
finYear100 <- SHYMNREPS %>% group_by(YearF, nSamples) %>% summarise(RBIAS = mean(bias)/mean(mnTrue)*100, RRMSE = sqrt(mean(bias^2))/mean(mnTrue)*100)

### Then make a plot for RRMSE: RRMSE is equivalent to a CV and values <20% are considered good, in general - i.e., on average, your estimated mean is within 20% of the true mean
YearsRRMSE100 <- ggplot(finYear100, aes(x = nSamples, y = RRMSE, group = YearF, color = YearF)) + geom_point(stat = "identity") + geom_line() + theme_bw() + theme(legend.title = element_blank(), legend.position = "bottom", panel.grid = element_blank()) + scale_y_continuous(limits = c(0,15), breaks = seq(0,15,1)) + scale_x_continuous(limits = c(0,100), breaks = seq(0,100,5)) + labs(x = "Number of individuals per quadrat", y = "Relative Root Mean Squared Error (% of true mean)") + theme(axis.text = element_text(color = "black")) + geom_hline(yintercept = 10, linetype = "dashed")
ggsave("YearsRRMSE_N100.png", YearsRRMSE100, dpi = 600, height = 5, width = 10)

### And one for RBIAS: generally these are unbiased; but these are also expressed as a percentage
YearsRBIAS100 <- ggplot(finYear100, aes(x = nSamples, y = RBIAS, group = YearF, color = YearF)) + geom_point(stat = "identity") + geom_line() + theme_bw() + theme(legend.title = element_blank(), legend.position = "bottom", panel.grid = element_blank()) + scale_y_continuous(limits = c(-1,1), breaks = seq(-1,1,0.1)) + geom_hline(yintercept = 0, linetype = "dashed") + scale_x_continuous(limits = c(0,100), breaks = seq(0,100,5)) + labs(x = "Number of individuals per quadrat", y = "Relative bias (% of true mean)") + labs(x = "Number of individuals per quadrat", y = "Relative bias (% of true mean)") + theme(axis.text = element_text(color = "black")) + geom_hline(yintercept = 0, linetype = "dashed")
ggsave("YearsRBIAS_N100.png", YearsRBIAS100, dpi = 600, height = 5, width = 10)

###################################
### BY YEAR AND QUARTER N ==100 ###
###################################
SHResample <- SH %>% group_by(GRP) %>% mutate(N = 1:n(), maxN = max(N)) %>% filter(maxN==100) %>% droplevels(.) %>% ungroup() %>% group_by(YearF, Quarter) %>% mutate(mnTrue = mean(SH))

SHYQMN <- NULL
SHYQMNREPS <- NULL
maxSamp <- seq(5, 100, 5)
nSims <- 1000
quadrats <- unique(SHResample$GRP)
for(i in 1:length(maxSamp)){
  for (j in 1:nSims){
    SHTemp <- SHResample %>% group_by(GRP) %>% sample_n(size = maxSamp[i], replace = T)
    SHYQMN <- suppressMessages(SHTemp %>% group_by(YearF, Quarter, mnTrue) %>% summarise(mnSH = mean(SH)) %>% mutate(bias = mnSH - mnTrue, nSamples = maxSamp[i], replicate = j))
    SHYQMNREPS <- rbind(SHYQMNREPS, SHYQMN)
    if(j %% 1000==0){cat(paste0("iteration: ", j, "; nSamples: ", maxSamp[i], "\n"))}
  }
}
finYQ100 <- SHYQMNREPS %>% group_by(YearF, Quarter, nSamples) %>% summarise(RBIAS = mean(bias)/mean(mnTrue)*100, RRMSE = sqrt(mean(bias^2))/mean(mnTrue)*100)

QuartersRRMSE100 <- ggplot(finYQ100, aes(x = nSamples, y = RRMSE, group = Quarter, color = Quarter)) + geom_point(stat = "identity") + geom_line() + facet_grid(~YearF) + theme(legend.position = "bottom") + theme_bw() + theme(legend.title = element_blank(), legend.position = "bottom", panel.grid = element_blank()) + scale_y_continuous(limits = c(0,15), breaks = seq(0,15,1)) + scale_x_continuous(limits = c(0,100), breaks = seq(0,100,5)) + labs(x = "Number of individuals per quadrat", y = "Relative Root Mean Squared Error (% of true mean)") + theme(axis.text = element_text(color = "black")) + geom_hline(yintercept = 10, linetype = "dashed")
ggsave("QuartersRRMSE_N100.png", QuartersRRMSE100, dpi = 600, height = 5, width = 10)

QuartersRBIAS100 <- ggplot(finYQ100, aes(x = nSamples, y = RBIAS, group = Quarter, color = Quarter)) + geom_point(stat = "identity") + geom_line() + facet_grid(~YearF) + theme(legend.position = "bottom") + theme_bw() + theme(legend.title = element_blank(), legend.position = "bottom", panel.grid = element_blank()) + scale_y_continuous(limits = c(-1,1), breaks = seq(-1,1,0.1)) + scale_x_continuous(limits = c(0,100), breaks = seq(0,100,5)) + labs(x = "Number of individuals per quadrat", y = "Relative Bias (% of true mean)") + theme(axis.text = element_text(color = "black")) + geom_hline(yintercept = 0, linetype = "dashed")
ggsave("QuartersRBIAS_N100.png", QuartersRBIAS100, dpi = 600, height = 5, width = 10)

