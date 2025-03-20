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

#### This code is AB specific and is run w/ w.2018, w.2021 and all 2020 removed to avoid skewing results
#### WQ data is included for all years

AB.CI <-read.csv('ConditionIndexAB2016-2022.csv', skip = 0, na= "Z")
View(AB.CI)
AB.WQ <-read.csv('ABWQ_Coll16-22.csv', skip = 0, na="Z")
View(AB.WQ)

#Removing columns not needed
AB.WQ$Sample.Depth<-NULL
AB.WQ$Turbidity.YSI<-NULL
AB.CI$Location.ID<-NULL
AB.CI$SampleEventID<-NULL
AB.WQ$Location.ID<-NULL
AB.CI$Shell.Height<-NULL
AB.CI$Shell.Length<-NULL
AB.CI$Shell.Width<-NULL
AB.CI$Total.Weight<-NULL
AB.CI$Tare.Pan.Weight <- NULL
AB.CI$Tissue.Wet.Weight <-NULL
AB.CI$Shell.Wet.Weight<-NULL
AB.CI$Tissue.Dry.Weight<-NULL
AB.CI$Final.tissue <-NULL
AB.CI$Shell.Dry.Weight <-NULL

# May need to convert to numeric - 
str(AB.CI)
AB.CI$CI <- as.numeric(AB.CI$CI)

##There are NAs in the WQ dataset but removing them from one column removes all the rows
#for the others so we're going to just keep them. This will at least make them numeric
str(AB.WQ)
AB.WQ$DO[AB.WQ$DO == "NA"] <- NA
AB.WQ$DO <- as.numeric(AB.WQ$DO)

AB.WQ$pH[AB.WQ$pH == "NA"] <- NA
AB.WQ$pH <- as.numeric(AB.WQ$pH)

AB.WQ$Temp[AB.WQ$Temp == "NA"] <- NA
AB.WQ$Temp <- as.numeric(AB.WQ$Temp)

AB.WQ$Salinity[AB.WQ$Salinity == "NA"] <- NA
AB.WQ$Salinity <- as.numeric(AB.WQ$Salinity)

##################
##### Dates
class(AB.CI$Date)
AB.CI$Date <- as.Date(AB.CI$Date, format = "%m/%d/%Y")
class(AB.WQ$Date) ##since WQ and CI are one df now, should be ok to not do this again
unique(AB.WQ$Date)
AB.WQ$Date <- as.Date(AB.WQ$Date, format = "%m/%d/%Y")


## to make WQ match CI observations:
station_counts <- AB.CI %>%
  group_by(Date, Region, Station) %>%
  summarise(Repeat = n(), .groups = 'drop')

WQ_avg <- AB.WQ %>%
  group_by(Date, Region, Station) %>%  
  summarise(across(where(is.numeric), mean, na.rm = TRUE), .groups = "drop") 

class(WQ_avg$Date)
WQ_avg$Date <- as.Date(WQ_avg$Date, format = "%m/%d/%Y")

AB.CI <- AB.CI %>%
  left_join(WQ_avg, by = c("Date", "Region", "Station"))

# Separating into bay Regions
AB.CI.West <- subset(AB.CI,Region == "W")
AB.CI.Central <- subset(AB.CI,Region == "C")
AB.CI.East <- subset(AB.CI,Region == "E")

AB.CI <- rbind(
  transform(AB.CI.West, Group = "West Region"),
  transform(AB.CI.Central, Group = "Central Region"),
  transform(AB.CI.East, Group = "East Region")
)
Region.order <- c("West", "Central", "East")

AB.CI_year_breaks <- seq(as.Date("2016-01-01"), as.Date("2022-01-01"), by = "years")
year_labels <- format(AB.CI_year_breaks, "%Y")

Region.colors <- c("West Region" = "#6D9DFF", 
                    "Central Region" = "#FFC368", 
                    "East Region" = "#FF7979")

AB.CI$Group <- factor(AB.CI$Group, levels = c("West Region", "Central Region", "East Region"))

ggplot(AB.CI %>% drop_na(CI), aes(x = Date, y = CI, color = Group))+
  geom_point() +
  labs(title = "Condition Index of Eastern Oysters in Apalachicola Bay Over Time", 
       x = "Date", y = "Condition Index",
       color = "Region") +
  scale_color_manual(name = "Region", 
                     values = Region.colors,
                     labels = c("West", "Central", "East")) +  # Assigning specific colors and labels
  theme_minimal() +
  scale_x_date(breaks = AB.CI_year_breaks, labels = year_labels)
  
##do one for each WQ parameter
ggplot(AB.CI, aes(x = Date, y = Temp, color = Group))+   
  geom_point() +
  labs(title = "Temperature of Apalachicola Bay Over Time", 
       x = "Date", y = "Temperature",
       color = "Region") +
  scale_color_manual(name = "Region", 
                     values = Region.colors,
                     labels = c("West", "Central", "East")) +  # Assigning specific colors and labels
  theme_minimal() +
  scale_x_date(breaks = AB.CI_year_breaks, labels = year_labels)

ggplot(AB.CI, aes(x = Date, y = Salinity, color = Group)) +
  geom_point() +
  labs(title = "Salinity of Apalachicola Bay Over Time", 
       x = "Date", y = "Salinity",
       color = "Region") +
  scale_color_manual(name = "Region", 
                     values = Region.colors,
                     labels = c("West", "Central", "East")) +  # Assigning specific colors and labels
  theme_minimal() +
  scale_x_date(breaks = AB.CI_year_breaks, labels = year_labels)

#Scatterplots don't seem to work as well for the other 2 so line graphs instead?
ggplot(AB.CI, aes(x = Date, y = DO, color = Group)) +
  geom_line() +
  labs(title = "Dissolved Oxygen of Apalachicola Bay Over Time", 
       x = "Date", y = "D.O.",
       color = "Region") +
  scale_color_manual(name = "Region", 
                     values = Region.colors,
                     labels = c("West", "Central", "East")) +  # Assigning specific colors and labels
  theme_minimal() +
  scale_x_date(breaks = AB.CI_year_breaks, labels = year_labels)

ggplot(AB.CI, aes(x = Date, y = pH, color = Group)) +
  geom_line() +
  labs(title = "pH of Apalachicola Bay Over Time", 
       x = "Date", y = "pH",
       color = "Region") +
  scale_color_manual(name = "Region", 
                     values = Region.colors,
                     labels = c("West", "Central", "East")) +  # Assigning specific colors and labels
  theme_minimal() +
  scale_x_date(breaks = AB.CI_year_breaks, labels = year_labels)

##### Modify a few variables (Year = year and GroupYear is a combined Group × Year
##### variable as a factor). Including Year as a factor is probably the best way to
##### incorporate temporal change in this analysis
AB.CI$Year <- year(AB.CI$Date)
AB.CI$YearF <- factor(AB.CI$Year)

AB.CI$GroupYear <- interaction(AB.CI$Group, AB.CI$Year) %>% droplevels(.)

AB.CI$pre.post <- factor(AB.CI$Year, levels = c("2016", "2017", "2018", "2019", "2021", "2022", "2023"), 
                         labels = c("Pre-closure", "Pre-closure", "Pre-closure", "Pre-closure", "Post-closure", 
                                    "Post-closure", "Post-closure"))
##### Then combine with Groups to make a new variable, pre_post_group
AB.CI$pre.post.group <- interaction(AB.CI$pre.post, AB.CI$Group)

#### Omit the 2020 data because it's not really all that comparable to the other year, also w.2018 and w.2021
AB.CI <- AB.CI %>% filter(Year !="2020") %>% droplevels(.) 
AB.CI <- AB.CI %>% filter(GroupYear != "West Region.2018", GroupYear != "West Region.2021") %>% droplevels(.)

###Separate scatterplots then stacked - CI
(plot_west.ci <- ggplot(subset(AB.CI, Group == "West Region"), aes(x = Date, y = CI)) +
    geom_point(color = Region.colors["West Region"]) +
    labs(title = "West Region", x = "Date", y = "Condition Index") +
    scale_y_continuous(limits = c(0, 10)) +  
    theme_minimal())

(plot_central.ci <- ggplot(subset(AB.CI, Group == "Central Region"), aes(x = Date, y = CI)) +
    geom_point(color = Region.colors["Central Region"]) +
    labs(title = "Central Region", x = "Date", y = "Condition Index") + 
    scale_y_continuous(limits = c(0, 10)) +
    theme_minimal())

(plot_east.ci <- ggplot(subset(AB.CI, Group == "East Region"), aes(x = Date, y = CI)) +
    geom_point(color = Region.colors["East Region"]) +
    labs(title = "East Region", x = "Date", y = "Condition Index") + 
    scale_y_continuous(limits = c(0, 10)) +
    theme_minimal())

# Combine the three scatterplots into a single plot using cowplot::plot_grid()
combined_plot.ci <- plot_grid(plot_west.ci, plot_central.ci, plot_east.ci, ncol = 1)
combined_plot.ci


### on to the Gamma analysis ###
## Scale WQ parameters? Taking awhile
scDO <- scale(AB.CI$DO)
scTemp <- scale(AB.CI$Temp)
scSalinity <- scale(AB.CI$Salinity)
scpH <- scale(AB.CI$pH)
#
#
table(is.na(AB.CI))
## We have some NA values so good to remove them first (most functions will remove them automatically)
AB.CI <- AB.CI %>% filter(complete.cases(.))

##### Fit gamma regressions with the fixed effect GroupYear 
mGamma.A <- glmmTMB(CI ~ GroupYear, family = Gamma(link = "log"), data = AB.CI)

##### Fit a model that ignores the spatial groupings and just assesses temporal change (need to convert
##### Year to a factor first)
AB.CI$Year <- factor(AB.CI$Year)
mGamma.B <- glmmTMB(CI ~ Year, family = Gamma(link = "log"), data = AB.CI)

##### Fit a gamma regressions with pre.post.group as a predictor 
mGamma.C <- glmmTMB(CI ~ pre.post.group, family = Gamma(link = "log"), data = AB.CI)

##### Fit gamma regressions with the fixed effect Temp 
mGamma.D <- glmmTMB(CI ~ Temp, family = Gamma(link = "log"), data = AB.CI)

##### Fit gamma regressions with the fixed effect Sal 
mGamma.E <- glmmTMB(CI ~ Salinity, family = Gamma(link = "log"), data = AB.CI)

##### Fit gamma regressions with the fixed effect DO 
mGamma.F <- glmmTMB(CI ~ DO, family = Gamma(link = "log"), data = AB.CI)

##### Fit gamma regressions with the fixed effect pH 
mGamma.G <- glmmTMB(CI ~ pH, family = Gamma(link = "log"), data = AB.CI)

##### Fit gamma regressions with the fixed effect all wq parameters 
mGamma.H <- glmmTMB(CI ~ Temp + Salinity + DO + pH + GroupYear+ (1 | Station) +(1 | Date), ##Since we are more section focused, should I change Station to Section?
                    family = Gamma(link = "log"), 
                    data = AB.CI)
mGamma.Hb <- glmmTMB(CI ~ Temp + Salinity + DO + pH + GroupYear+ (1 | Station) +(1 | Date), 
                     dispformula = ~ GroupYear, #Modeling dispersion as a function of fixed and random effects 
                     family = Gamma(link = "log"), 
                     data = AB.CI)

#### Use AICc to compare the 9 models: 
model.names <- c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'Hb')
aictab(list(mGamma.A, mGamma.B, mGamma.C, mGamma.D, mGamma.E, mGamma.F, mGamma.G, mGamma.H, mGamma.Hb), 
       modnames = model.names)
#Model Hb lowest AICc

##### Assess goodness-of-fit 
simulateResiduals(mGamma.Hb, n = 1000, plot = T)

##### Look at estimates
summary(mGamma.Hb)

## Varying one variable at a time; Station and Date = NA to make population-level predictions (i.e., predictions that ignore the random effects subjects)
newDat <- expand.grid(GroupYear = unique(AB.CI$GroupYear), Temp = unique(AB.CI$Temp, na.rm = T), 
                      Salinity = mean(AB.CI$Salinity, na.rm = T), DO = mean(AB.CI$DO, na.rm = T), 
                      pH = mean(AB.CI$pH, na.rm = T), Station = NA, Date = NA) #Again, should Station be Region/Section instead?

##### Gamma: Now make some predictions, calculating the back-transformed mean and 95% confidence intervals by exponentiating the log-scale parameter estimates/predictions
predsGamma <- predict(mGamma.Hb, newdata = newDat, type = "link", se.fit = T)
preddsGamma <- data.frame(newDat, fit = predsGamma$fit, se = predsGamma$se.fit)
preddsGamma$mean <- exp(preddsGamma$fit)
preddsGamma$lwr <- exp(preddsGamma$fit - 1.96*preddsGamma$se)
preddsGamma$upr <- exp(preddsGamma$fit + 1.96*preddsGamma$se)
preddsGamma <- preddsGamma %>% separate(GroupYear, into = c("Region", "Group", "Year"), remove = FALSE) %>%
  mutate(LocationGroup = paste0(Region," ", Group))

desired_order <- c("West Region", "Central Region", "East Region")
preddsGamma$LocationGroup <- factor(preddsGamma$LocationGroup, levels = desired_order)
preddsGamma$pre.post <- ifelse(preddsGamma$Year < 2020, "pre", "post")
preddsGamma$Yearn<- as.numeric(preddsGamma$Year)

## Following graph has color coded facet strips to match the map and scatter plots above
strip <- strip_themed(background_x = elem_list_rect(fill = Region.colors))

bar.colors <- c("2016"= "#FFBBBB",
                "2017"="#56B4E9",
                "2018"="#CC79A7",
                "2019"="#FF7000",
                "2021"= "#0067A0",
                "2022"="#94DCC8")

(gammaPlot.a <- ggplot(preddsGamma, aes(x = Yearn, y = mean, fill = Year)) + 
    geom_bar(stat = "identity", color = "black") + 
    scale_fill_manual(values = bar.colors) +
    geom_errorbar(aes(ymin = lwr, ymax = upr),width = 0.25) + 
    facet_wrap(~LocationGroup) +
    facet_wrap2(~ LocationGroup, strip = strip) +
    theme_bw() + 
    theme(axis.text = element_text(color = "black"), 
          panel.grid.major.y = element_line(colour = "grey90"), 
          panel.grid.minor.y = element_line(colour = "grey90", linetype = "dashed"), 
          panel.grid.major.x = element_blank()) + 
    scale_y_continuous(limits = c(0, 4), 
                       breaks = seq(0, 4, 0.25), expand = expansion(add = c(0, 0))) + 
    labs(x = NULL, y = "Gamma: Mean CI", title = "Condition Index of Eastern Oysters in Apalachicola Bay: Pre- & Post- Fishery Closure") +
    geom_vline(data = preddsGamma, aes(xintercept = 2020), linetype = "dashed", color = "#4B0082"))

### Because of the station numbers being added to the AB.CI df, looks like maybe each stations 
### error margin is included (lots of tick marks)??  
### Also my annual CI average is higher/lower than before?? 

#### Another way to do above but also change font color
# Define custom theme with strip_themed() for background color and strip.text for text color
#custom_theme <- function() {
#  theme_bw() +
#    theme(strip.background = element_rect(fill = Region.colors),
#          strip.text = element_text(color = "white"))
# }

# Apply custom theme to the plot with facet_wrap2()
#(gammaPlot.b <- ggplot(preddsGamma, aes(x = Yearn, y = mean, fill = Year)) + 
#  geom_bar(stat = "identity", color = "black") + 
#  geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.25) + 
#  facet_wrap2(~LocationGroup, strip = strip_themed(background_x = elem_list_rect(fill = Region.colors))) +
#  custom_theme() + 
#  theme(axis.text = element_text(color = "black"), 
#        panel.grid.major.y = element_line(colour = "grey90"), 
#        panel.grid.minor.y = element_line(colour = "grey90", linetype = "dashed"), 
#        panel.grid.major.x = element_blank()) + 
#  scale_y_continuous(limits = c(0,4), 
#                     breaks = seq(0,4, 0.25), expand = expansion(add = c(0,0))) + 
#  labs(x = NULL, y = "Gamma: Mean CI", 
#       title = "Condition Index of Eastern Oysters in Apalachicola Bay: Pre- & Post- Fishery Closure") +
#  geom_vline(data = preddsGamma, aes(xintercept = 2020), linetype = "dashed", color = "#4B0082"))

##### Summarize the raw data (useful to compare to model estimates)
c_summ <- AB.CI %>% group_by(GroupYear) %>% summarise(N = n(), mn = mean(CI, na.rm = T), se = sd(CI, na.rm = T)/sqrt(N), lwr = mn-1.96*se, upr = mn+1.96*se) %>% separate(GroupYear, into = c("Region", "Group", "Year"), remove = FALSE) %>% mutate(LocationGroup = paste0(Region," ", Group))
desired.order <- c("West", "Central", "East")
c_summ$LocationGroup <- factor(c_summ$LocationGroup, levels = desired_order)
c_summ <- c_summ %>% arrange(LocationGroup, Year)

##### Compare estimated/predicted means to those observed in the raw data: 
(rawdataPlot <- ggplot(c_summ, aes(x = Year, y = mn, fill = Year)) + 
    geom_bar(stat = "identity", color = "black") + geom_errorbar(aes(ymin = lwr, ymax = upr), width = 0.25) + 
    facet_wrap(~LocationGroup) + theme_bw() + theme(axis.text = element_text(color = "black"), 
                                                    panel.grid.major.y = element_line(colour = "grey90"), 
                                                    panel.grid.minor.y = element_line(colour = "grey90", linetype = "dashed"), 
                                                    panel.grid.major.x = element_blank()) + 
    scale_y_continuous(limits = c(0,4), breaks = seq(0,4, 0.25), expand = expansion(add = c(0,0))) + 
  labs(x = NULL, y = "Raw data: Mean Condition Index"))
#### does not match

##### Compare means among years for each region
(emmGamma <- data.frame(emmeans(mGamma.Hb, specs = pairwise ~ GroupYear)$emmeans)) # these are marginal or group-specific means
(conGamma <- data.frame(emmeans(mGamma.Hb, specs = pairwise ~ GroupYear)$contrasts)) # these are comparisons among groups

##### Break contrasts down by region and year...the contrasts above (con) are to extensive because 
##### there are contrasts we don't care about, so this code reduced the number of comparisons and then 
##### adjusts p-values accordingly for multiple comparisons. Here, we are only comparing among years for
##### each group (west, central, and east)
(Gamma_contrasts <- conGamma %>% separate(contrast, into = c("contrast1", "contrast2"), sep = "\\-", extra = "merge") %>% 
    separate(contrast1, into = c("Group1", "Year1"), sep = "\\.", extra = "merge") %>% mutate(contrast2 = trimws(contrast2)) %>% 
    separate(contrast2, into = c("Group2", "Year2"), sep = "\\.", extra = "merge") %>% 
    mutate(keep = case_when(Group1==Group2 ~ "Keep", .default = "Drop")) %>% filter(keep == "Keep") %>%
    mutate(p_adj = p.adjust(p.value, method = "fdr")) %>% arrange(Group1, Group2) %>% 
    mutate(p.value = round(p.value, 3), p_adj = round(p_adj, 3)))
# Convert Group1 and Group2 to factors with desired order
(Gamma_contrasts <- Gamma_contrasts %>%
    mutate(Group1 = factor(Group1, levels = desired_order),
           Group2 = factor(Group2, levels = desired_order)) %>%
    arrange(Group1, Group2))

##### Filter for statistically important (not necessarily biologically important) pairwise comparisons
(gamma_important_contrasts <- Gamma_contrasts[Gamma_contrasts$p_adj<0.05,])
## nothing significant so zero rows

# Save contrast outputs as tab-delimited text
#write.table(gamma_important_contrasts, "Gamma_contrastsa.txt", sep = "\t", row.names = FALSE)

##Make a table for CI means
means <- aggregate(CI ~ AB.CI$GroupYear, data = AB.CI, FUN= mean, levels = desired_order)
means
write.table(means, "means.a.txt", sep = "\t", row.names = FALSE)
#### 
## You could maybe have, say, 3 levels of another variable like so 
## (eg. min, mean, and max of salinity) if you want to show the effect of two 
## continuous predictors...
newDat2 <- expand.grid(GroupYear = unique(AB.CI$GroupYear), 
              Temp = unique(AB.CI$Temp, na.rm = T), 
              Salinity = c(min(AB.CI$Salinity, na.rm = T), mean(AB.CI$Salinity, na.rm = T), max(AB.CI$Salinity, na.rm = T)), 
              DO = mean(AB.CI$DO, na.rm = T),
              pH = mean(AB.CI$pH, na.rm = T), 
              Station = NA, 
              Date = NA
              ) 


