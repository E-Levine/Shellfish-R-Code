## CI and WQ data organization and plots  ##
#
# This code is to prepare the cleaned data sets for gamma analysis and create plots used for reports
#
library(Matrix)
library(glmmTMB)
library(DHARMa)
library(tidyverse)
library(lubridate)
library(emmeans)
library(AICcmodavg)
library(ggplot2)
library(dplyr)
library(cowplot)
library(ggridges)
library(ggh4x)

#### This code is AB specific and removes w.2018, w.2021 and all 2020 to avoid skewing results
#### WQ data is included for all years

AB.CI <-read.csv('Data/Cleaned_Data/AB_filtered_CI_2016_2023.csv', skip = 0, na= "Z")
View(AB.CI)
AB.WQ <-read.csv('Data/Cleaned_Data/AB_filtered_WQ_COLL_2016_2023.csv', skip = 0, na="Z")
View(AB.WQ)

#Removing columns not needed
AB.CI$Location.ID<-NULL
AB.CI$SampleEventID<-NULL
AB.CI$ShellHeight<-NULL
AB.CI$ShellLength<-NULL
AB.CI$ShellWidth<-NULL

AB.WQ$SampleEvent<-NULL
AB.WQ$Depth<-NULL
#
#### Calculate CI 
# CI = (dry tissue weight/dry shell weight)*100
AB.CI$FinalTissueWt <- AB.CI$TissueDryWeight - AB.CI$TarePanWeight
#
# Convert columns to numeric
AB.CI$FinalTissueWt <- as.numeric(AB.CI$FinalTissueWt)
AB.CI$ShellDryWeight <- as.numeric(AB.CI$ShellDryWeight)
# Remove rows with missing values
AB.CI <- na.omit(AB.CI)
#
# Perform calculation and add CI as a column
CI <- (AB.CI$FinalTissueWt / AB.CI$ShellDryWeight) * 100
AB.CI$CI<- CI
#
str(AB.CI)  #checking to see if we need to convert anything to numeric - 
str(AB.WQ)

##There are NAs in the WQ dataset but removing them from one column removes all the rows
#for the others so we're going to just keep them. This will at least make them numeric
AB.WQ$DO[AB.WQ$DO == "NA"] <- NA
AB.WQ$DO <- as.numeric(AB.WQ$DO)

AB.WQ$pH[AB.WQ$pH == "NA"] <- NA
AB.WQ$pH <- as.numeric(AB.WQ$pH)

AB.WQ$Temperature[AB.WQ$Temperature == "NA"] <- NA
AB.WQ$Temperature <- as.numeric(AB.WQ$Temperature)

AB.WQ$Salinity[AB.WQ$Salinity == "NA"] <- NA
AB.WQ$Salinity <- as.numeric(AB.WQ$Salinity)
#
#
##### Dates need to convert to Date format
head(AB.CI$Date) #shows how dates are formatted (- vs /)
head(AB.WQ$Date)
AB.CI$Date <- as.Date(AB.CI$Date, format = "%Y-%m-%d")
AB.WQ$Date <- as.Date(AB.WQ$Date, format = "%Y-%m-%d")
#
## Make WQ match CI observations:
station_counts <- AB.CI %>%
  group_by(Date, Section, Station) %>%
  summarise(Repeat = n(), .groups = 'drop')

WQ_avg <- AB.WQ %>%
  group_by(Date, Section, Station) %>%  
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop") 

AB.CIwq <- AB.CI %>%
  left_join(WQ_avg, by = c("Date", "Section", "Station"))
# 
# Separating into bay Sections
AB.CIwq.West <- subset(AB.CIwq,Section == "W")
AB.CIwq.Central <- subset(AB.CIwq,Section == "C")
AB.CIwq.East <- subset(AB.CIwq,Section == "E")

AB.CIwq <- rbind(
  transform(AB.CIwq.West, Group = "West Section"),
  transform(AB.CIwq.Central, Group = "Central Section"),
  transform(AB.CIwq.East, Group = "East Section")
)
Section.order <- c("West", "Central", "East")

AB.CIwq_year_breaks <- seq(as.Date("2016-01-01"), as.Date("2024-01-01"), by = "years")
year_labels <- format(AB.CIwq_year_breaks, "%Y")

AB.CIwq$Group <- factor(AB.CIwq$Group, levels = c("West Section", "Central Section", "East Section"))

# Assign CB friendly colors to sections and years
Section.colors <- c("West Section" = "#00a884", 
                    "Central Section" = "#FF7979", 
                    "East Section" = "#6D9DFF")
yr.colors <- c("2016"="#FFBBBB",
               "2017"="#56B4E9",
               "2018"="#CC79A7",
               "2019"="#0067A0",
               "2021"="#58FF7E",
               "2022"="#FF6666",
               "2023"="#94DCC8")
strip <- strip_themed(background_x = elem_list_rect(fill = yr.colors))
strip2 <- strip_themed(background_x = elem_list_rect(fill = Section.colors))
desired_order <- c("West Section", "Central Section", "East Section")

#
## Plot
ggplot(AB.CIwq %>% drop_na(CI), aes(x = Date, y = CI, color = Group))+
  geom_point() +
  labs(title = "Condition Index of Eastern Oysters in Apalachicola Bay Over Time", 
       x = "Date", y = "Condition Index",
       color = "Section") +
  scale_color_manual(name = "Section", 
                     values = Section.colors,
                     labels = c("West", "Central", "East")) +  # Assigning specific colors and labels
  theme_minimal() +
  scale_x_date(breaks = AB.CIwq_year_breaks, labels = year_labels)
  
##do one for each WQ parameter
ggplot(AB.CIwq, aes(x = Date, y = Temperature, color = Group))+   
  geom_point() +
  labs(title = "Temperature of Apalachicola Bay Over Time", 
       x = "Date", y = "Temperature",
       color = "Section") +
  scale_color_manual(name = "Section", 
                     values = Section.colors,
                     labels = c("West", "Central", "East")) +  
  theme_minimal() +
  scale_x_date(breaks = AB.CIwq_year_breaks, labels = year_labels)

ggplot(AB.CIwq, aes(x = Date, y = Salinity, color = Group)) +
  geom_point() +
  labs(title = "Salinity of Apalachicola Bay Over Time", 
       x = "Date", y = "Salinity",
       color = "Section") +
  scale_color_manual(name = "Section", 
                     values = Section.colors,
                     labels = c("West", "Central", "East")) + 
  theme_minimal() +
  scale_x_date(breaks = AB.CIwq_year_breaks, labels = year_labels)

#Scatterplots don't seem to work as well for the other 2 so line graphs instead
ggplot(AB.CIwq, aes(x = Date, y = DO, color = Group)) +
  geom_line() +
  labs(title = "Dissolved Oxygen of Apalachicola Bay Over Time", 
       x = "Date", y = "D.O.",
       color = "Section") +
  scale_color_manual(name = "Section", 
                     values = Section.colors,
                     labels = c("West", "Central", "East")) +  
  theme_minimal() +
  scale_x_date(breaks = AB.CIwq_year_breaks, labels = year_labels)

ggplot(AB.CIwq, aes(x = Date, y = pH, color = Group)) +
  geom_line() +
  labs(title = "pH of Apalachicola Bay Over Time", 
       x = "Date", y = "pH",
       color = "Section") +
  scale_color_manual(name = "Section", 
                     values = Section.colors,
                     labels = c("West", "Central", "East")) + 
  theme_minimal() +
  scale_x_date(breaks = AB.CIwq_year_breaks, labels = year_labels)

##### Modify a few variables (Year = year and GroupYear is a combined Group × Year
##### variable as a factor). Including Year as a factor is probably the best way to
##### incorporate temporal change in this analysis
AB.CIwq$Year <- year(AB.CIwq$Date)
AB.CIwq$YearF <- factor(AB.CIwq$Year)

AB.CIwq$GroupYear <- interaction(AB.CIwq$Group, AB.CIwq$Year) %>% droplevels(.)

AB.CIwq$pre.post <- factor(AB.CIwq$Year, levels = c("2016", "2017", "2018", "2019", "2021", "2022", "2023"), 
                         labels = c("Pre-closure", "Pre-closure", "Pre-closure", "Pre-closure", "Post-closure", 
                                    "Post-closure", "Post-closure"))
##### Then combine with Groups to make a new variable, pre_post_group
AB.CIwq$pre.post.group <- interaction(AB.CIwq$pre.post, AB.CIwq$Group)

#### Omit the 2020, w.2018 and w.2021 data to avoid skewing results
AB.CIwq <- AB.CIwq %>% filter(Year !="2020") %>% droplevels(.) 
AB.CIwq <- AB.CIwq %>% filter(GroupYear != "West Section.2018", GroupYear != "West Section.2021") %>% droplevels(.)

###Separate scatterplots then stacked - CI
(plot_west.ci <- ggplot(subset(AB.CIwq, Group == "West Section"), aes(x = Date, y = CI)) +
    geom_point(color = Section.colors["West Section"]) +
    labs(title = "West Section", x = "Date", y = "Condition Index") +
    scale_y_continuous(limits = c(0, 10)) +  
    theme_minimal())

(plot_central.ci <- ggplot(subset(AB.CIwq, Group == "Central Section"), aes(x = Date, y = CI)) +
    geom_point(color = Section.colors["Central Section"]) +
    labs(title = "Central Section", x = "Date", y = "Condition Index") + 
    scale_y_continuous(limits = c(0, 10)) +
    theme_minimal())

(plot_east.ci <- ggplot(subset(AB.CIwq, Group == "East Section"), aes(x = Date, y = CI)) +
    geom_point(color = Section.colors["East Section"]) +
    labs(title = "East Section", x = "Date", y = "Condition Index") + 
    scale_y_continuous(limits = c(0, 10)) +
    theme_minimal())

# Combine the three scatterplots into a single plot using cowplot::plot_grid()
combined_plot.ci <- plot_grid(plot_west.ci, plot_central.ci, plot_east.ci, ncol = 1)
combined_plot.ci

### Gamma plots ###
## Scale WQ parameters
str(AB.CIwq)
AB.CIwq$CI <- as.numeric(AB.CIwq$CI)
AB.CIwq$Temperature <- as.numeric(AB.CIwq$Temperature)
AB.CIwq$Salinity <- as.numeric(AB.CIwq$Salinity)
AB.CIwq$DO <- as.numeric(AB.CIwq$DO)
AB.CIwq$pH <- as.numeric(AB.CIwq$pH)

scDO <- scale(AB.CIwq$DO)
scTemp <- scale(AB.CIwq$Temperature)
scSalinity <- scale(AB.CIwq$Salinity)
scpH <- scale(AB.CIwq$pH)
#
#
table(is.na(AB.CIwq))
## We have some NA values so good to remove them first (most functions will remove them automatically)
AB.CIwq <- AB.CIwq %>% filter(complete.cases(.))

##### Fit gamma regressions with the fixed effect all wq parameters 
AB.CIwq$GYD <- droplevels(interaction(AB.CIwq$GroupYear, AB.CIwq$Date))
mGamma.A <- glmmTMB(CI ~ Temperature + Salinity + DO + pH + GroupYear + (1|GYD), 
                    dispformula = ~ GroupYear,
                    family = Gamma(link = "log"), 
                    data = AB.CIwq)
##### Fit gamma regressions with the fixed effect GroupYear to see year to year differences in CI only
mGamma.B <- glmmTMB(CI ~ GroupYear, family = Gamma(link = "log"), data = AB.CIwq)

### Use AICc to compare the 2 models: we'll use the model with the lowest score for making any inferences moving forward
model.names <- c('A', 'B')
aictab(list(mGamma.A, mGamma.B), modnames = model.names) #will use .A for analysis

##### Assess goodness-of-fit 
simulateResiduals(mGamma.A, n = 1000, plot = T)

##### Look at estimates
summary(mGamma.A)

## Varying one variable at a time; Station and Date = NA to make population-level predictions (i.e., predictions that ignore the random effects subjects)
## Temp
newDat.Temp <- expand.grid(GroupYear = unique(AB.CIwq$GroupYear), Temperature = unique(AB.CIwq$Temperature, na.rm = T), 
                           Salinity = mean(AB.CIwq$Salinity, na.rm = T), DO = mean(AB.CIwq$DO, na.rm = T), 
                           pH = mean(AB.CIwq$pH, na.rm = T), Station = NA, Date = NA, GYD = NA)
predsGamma.At <- predict(mGamma.A, newdata = newDat.Temp, type = "link", se.fit = T)
preddsGamma.At <- data.frame(newDat.Temp, fit = predsGamma.At$fit, se = predsGamma.At$se.fit)
preddsGamma.At$mean <- exp(preddsGamma.At$fit)
preddsGamma.At$lwr <- exp(preddsGamma.At$fit - 1.96*preddsGamma.At$se)
preddsGamma.At$upr <- exp(preddsGamma.At$fit + 1.96*preddsGamma.At$se)
preddsGamma.At <- preddsGamma.At %>% separate(GroupYear, into = c("Region", "Group", "Year"), remove = FALSE) %>%
  mutate(LocationGroup = paste0(Region," ", Group))

preddsGamma.At$LocationGroup <- factor(preddsGamma.At$LocationGroup, levels = desired_order)
preddsGamma.At$pre.post <- ifelse(preddsGamma.At$Year < 2020, "pre", "post")
preddsGamma.At$Yearn<- as.numeric(preddsGamma.At$Year)
# plot affects of temperature on CI over years
(gammaPlot.Temp <- ggplot(preddsGamma.At, aes(x = Temperature, y = mean, fill = LocationGroup)) +
    geom_line() + 
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.25) + 
    facet_wrap(~Year) +
    facet_wrap2(~Year, strip = strip)+
    theme_bw() + 
    theme(axis.text = element_text(color = "black"), 
          panel.grid.major.y = element_line(colour = "grey90"), 
          panel.grid.minor.y = element_line(colour = "grey90", linetype = "dashed"), 
          panel.grid.major.x = element_blank(), legend.title = element_blank(), legend.position = "bottom") + 
    scale_fill_manual(values = Section.colors)+
    scale_y_continuous(limits = c(0, 5), 
                       breaks = seq(0, 5, 0.5), expand = expansion(add = c(0, 0))) + 
    labs(x = "Mean temperature (°C)", y = "Mean condition index (± 95% CL)", 
         title = "Affect of Temperature on the Condition Index of Eastern Oysters in Apalachicola Bay"))
## Salinity
newDat.Sal <- expand.grid(GroupYear = unique(AB.CIwq$GroupYear), Temperature = mean(AB.CIwq$Temperature, na.rm = T), 
                          Salinity = unique(AB.CIwq$Salinity, na.rm = T), DO = mean(AB.CIwq$DO, na.rm = T), 
                          pH = mean(AB.CIwq$pH, na.rm = T), Station = NA, Date = NA, GYD = NA)
predsGamma.As <- predict(mGamma.A, newdata = newDat.Sal, type = "link", se.fit = T)
preddsGamma.As <- data.frame(newDat.Sal, fit = predsGamma.As$fit, se = predsGamma.As$se.fit)
preddsGamma.As$mean <- exp(preddsGamma.As$fit)
preddsGamma.As$lwr <- exp(preddsGamma.As$fit - 1.96*preddsGamma.As$se)
preddsGamma.As$upr <- exp(preddsGamma.As$fit + 1.96*preddsGamma.As$se)
preddsGamma.As <- preddsGamma.As %>% separate(GroupYear, into = c("Region", "Group", "Year"), remove = FALSE) %>%
  mutate(LocationGroup = paste0(Region," ", Group))

preddsGamma.As$LocationGroup <- factor(preddsGamma.As$LocationGroup, levels = desired_order)
preddsGamma.As$pre.post <- ifelse(preddsGamma.As$Year < 2020, "pre", "post")
preddsGamma.As$Yearn<- as.numeric(preddsGamma.As$Year)
#plot affect of salinity on CI over years 
(gammaPlot.Sal <- ggplot(preddsGamma.As, aes(x = Salinity, y = mean, fill = LocationGroup)) +
    geom_line() + 
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.25) + 
    facet_wrap(~Year) +
    facet_wrap2(~Year, strip = strip)+
    theme_bw() + 
    theme(axis.text = element_text(color = "black"), 
          panel.grid.major.y = element_line(colour = "grey90"), 
          panel.grid.minor.y = element_line(colour = "grey90", linetype = "dashed"), 
          panel.grid.major.x = element_blank(), legend.title = element_blank(), legend.position = "bottom") + 
    scale_fill_manual(values = Section.colors)+
    scale_y_continuous(limits = c(0, 5), 
                       breaks = seq(0, 5, 0.5), expand = expansion(add = c(0, 0))) + 
    labs(x = "Mean salinity (ppt)", y = "Mean condition index (± 95% CL)", 
         title = "Affect of Salinity on the Condition Index of Eastern Oysters in Apalachicola Bay"))
## DO
newDat.DO <- expand.grid(GroupYear = unique(AB.CIwq$GroupYear), Temperature = mean(AB.CIwq$Temperature, na.rm = T), 
                         Salinity = mean(AB.CIwq$Salinity, na.rm = T), DO = unique(AB.CIwq$DO, na.rm = T), 
                         pH = mean(AB.CIwq$pH, na.rm = T), Station = NA, Date = NA, GYD = NA)
predsGamma.Ad <- predict(mGamma.A, newdata = newDat.DO, type = "link", se.fit = T)
preddsGamma.Ad <- data.frame(newDat.DO, fit = predsGamma.Ad$fit, se = predsGamma.Ad$se.fit)
preddsGamma.Ad$mean <- exp(preddsGamma.Ad$fit)
preddsGamma.Ad$lwr <- exp(preddsGamma.Ad$fit - 1.96*preddsGamma.Ad$se)
preddsGamma.Ad$upr <- exp(preddsGamma.Ad$fit + 1.96*preddsGamma.Ad$se)
preddsGamma.Ad <- preddsGamma.Ad %>% separate(GroupYear, into = c("Region", "Group", "Year"), remove = FALSE) %>%
  mutate(LocationGroup = paste0(Region," ", Group))

preddsGamma.Ad$LocationGroup <- factor(preddsGamma.Ad$LocationGroup, levels = desired_order)
preddsGamma.Ad$pre.post <- ifelse(preddsGamma.Ad$Year < 2020, "pre", "post")
preddsGamma.Ad$Yearn<- as.numeric(preddsGamma.Ad$Year)
#plot affect of DO on CI over years 
(gammaPlot.DO <- ggplot(preddsGamma.Ad, aes(x = DO, y = mean, fill = LocationGroup)) +
    geom_line() + 
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.25) + 
    facet_wrap(~Year) +
    facet_wrap2(~Year, strip = strip)+
    theme_bw() + 
    theme(axis.text = element_text(color = "black"), 
          panel.grid.major.y = element_line(colour = "grey90"), 
          panel.grid.minor.y = element_line(colour = "grey90", linetype = "dashed"), 
          panel.grid.major.x = element_blank(), legend.title = element_blank(), legend.position = "bottom") + 
    scale_fill_manual(values = Section.colors)+
    scale_y_continuous(limits = c(0, 5), 
                       breaks = seq(0, 5, 0.5), expand = expansion(add = c(0, 0))) + 
    labs(x = "Mean Dissolved Oxygen (mg/L)", y = "Mean condition index (± 95% CL)", 
         title = "Affect of DO on the Condition Index of Eastern Oysters in Apalachicola Bay"))

## pH
newDat.pH <- expand.grid(GroupYear = unique(AB.CIwq$GroupYear), Temperature = mean(AB.CIwq$Temperature, na.rm = T), 
                         Salinity = mean(AB.CIwq$Salinity, na.rm = T), DO = mean(AB.CIwq$DO, na.rm = T), 
                         pH = unique(AB.CIwq$pH, na.rm = T), Station = NA, Date = NA, GYD = NA)
predsGamma.Ap <- predict(mGamma.A, newdata = newDat.pH, type = "link", se.fit = T)
preddsGamma.Ap <- data.frame(newDat.pH, fit = predsGamma.Ap$fit, se = predsGamma.Ap$se.fit)
preddsGamma.Ap$mean <- exp(preddsGamma.Ap$fit)
preddsGamma.Ap$lwr <- exp(preddsGamma.Ap$fit - 1.96*preddsGamma.Ap$se)
preddsGamma.Ap$upr <- exp(preddsGamma.Ap$fit + 1.96*preddsGamma.Ap$se)
preddsGamma.Ap <- preddsGamma.Ap %>% separate(GroupYear, into = c("Region", "Group", "Year"), remove = FALSE) %>%
  mutate(LocationGroup = paste0(Region," ", Group))

preddsGamma.Ap$LocationGroup <- factor(preddsGamma.Ap$LocationGroup, levels = desired_order)
preddsGamma.Ap$pre.post <- ifelse(preddsGamma.Ap$Year < 2020, "pre", "post")
preddsGamma.Ap$Yearn<- as.numeric(preddsGamma.Ap$Year)
#plot affect of pH on CI over years 
(gammaPlot.pH <- ggplot(preddsGamma.Ap, aes(x = pH, y = mean, fill = LocationGroup)) +
    geom_line() + 
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.25) + 
    facet_wrap(~Year) +
    facet_wrap2(~Year, strip = strip)+
    theme_bw() + 
    theme(axis.text = element_text(color = "black"), 
          panel.grid.major.y = element_line(colour = "grey90"), 
          panel.grid.minor.y = element_line(colour = "grey90", linetype = "dashed"), 
          panel.grid.major.x = element_blank(), legend.title = element_blank(), legend.position = "bottom") + 
    scale_fill_manual(values = Section.colors)+
    scale_y_continuous(limits = c(0, 5), 
                       breaks = seq(0, 5, 0.5), expand = expansion(add = c(0, 0))) + 
    labs(x = "Mean pH", y = "Mean condition index (± 95% CL)", 
         title = "Affect of pH on the Condition Index of Eastern Oysters in Apalachicola Bay"))

## NewDat.B to get CI over years
newDat.B <- expand.grid(GroupYear = unique(AB.CIwq$GroupYear), Section = NA)
predsGamma.B <- predict(mGamma.B, newdata = newDat.B, type = "link", se.fit = T)
preddsGamma.B <- data.frame(newDat.B, fit = predsGamma.B$fit, se = predsGamma.B$se.fit)
preddsGamma.B$mean <- exp(preddsGamma.B$fit)
preddsGamma.B$lwr <- exp(preddsGamma.B$fit - 1.96*preddsGamma.B$se)
preddsGamma.B$upr <- exp(preddsGamma.B$fit + 1.96*preddsGamma.B$se)
preddsGamma.B <- preddsGamma.B %>% separate(GroupYear, into = c("Region", "Group", "Year"), remove = FALSE) %>%
  mutate(LocationGroup = paste0(Region," ", Group))

preddsGamma.B$LocationGroup <- factor(preddsGamma.B$LocationGroup, levels = desired_order)
preddsGamma.B$pre.post <- ifelse(preddsGamma.B$Year < 2020, "pre", "post")
preddsGamma.B$Yearn<- as.numeric(preddsGamma.B$Year)
# plot CI change over years
(gammaPlot.B <- ggplot(preddsGamma.B, aes(x = Yearn, y = mean, fill = Year)) + 
    geom_bar(stat = "identity", color = "black") + 
    scale_fill_manual(values = yr.colors) +
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.25) + 
    facet_wrap(~LocationGroup) +
    facet_wrap2(~ LocationGroup, strip = strip2) +
    theme_bw() + 
    theme(axis.text = element_text(color = "black"), 
          panel.grid.major.y = element_line(colour = "grey90"), 
          panel.grid.minor.y = element_line(colour = "grey90", linetype = "dashed"), 
          panel.grid.major.x = element_blank()) + 
    scale_x_continuous(breaks = unique(preddsGamma.B$Yearn)) + 
    scale_y_continuous(limits = c(0,4), 
                       breaks = seq(0,4, 0.25), expand = expansion(add = c(0,0))) + 
    labs(x = NULL, y = "Gamma: Mean CI", title = "Condition Index of Eastern Oysters in Apalachicola Bay 2016-2023") +
    geom_vline(data = preddsGamma.B, aes(xintercept = 2020), linetype = "dashed", color = "#4B0082"))

##### Summarize the raw CI data (useful to compare to model estimates)
c_summ <- AB.CIwq %>% group_by(GroupYear) %>% 
  summarise(N = n(), mn = mean(CI, na.rm = T), se = sd(CI, na.rm = T)/sqrt(N), 
            lwr = mn-1.96*se, upr = mn+1.96*se) %>%
  separate(GroupYear, into = c("Section", "Group", "Year"), remove = FALSE) %>% 
  mutate(LocationGroup = paste0(Section," ", Group))
c_summ$LocationGroup <- factor(c_summ$LocationGroup, levels = desired_order)
c_summ <- c_summ %>% arrange(LocationGroup, Year)

##### Compare estimated/predicted CI means to those observed in the raw data: 
(rawdataPlot.CI <- ggplot(c_summ, aes(x = Year, y = mn, fill = Year)) + 
    geom_bar(stat = "identity", color = "black") + 
    geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.25) + 
    scale_fill_manual(values = yr.colors) +
    facet_wrap(~LocationGroup) + 
    theme_bw() + 
    theme(axis.text = element_text(color = "black"), 
          panel.grid.major.y = element_line(colour = "grey90"),
          panel.grid.minor.y = element_line(colour = "grey90", linetype = "dashed"), 
          panel.grid.major.x = element_blank()) + 
    scale_y_continuous(limits = c(0,4), breaks = seq(0,4, 0.25), expand = expansion(add = c(0,0))) + 
    labs(x = NULL, y = "Raw data: Mean Condition Index", title = "Raw data"))

### export AB.CIwq as new .csv and move on to new script for only analysis
write.csv(AB.CIwq, 
          file = paste0("Data/Cleaned_data/","AB_CIwq_2016_2023.csv"), 
          row.names = FALSE)

