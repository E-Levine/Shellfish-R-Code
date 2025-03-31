## CI Gamma Analysis  ##
#
# This code is to run the gamma regression model using the data sets previously prepped
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
AB.CIwq <-read.csv('Data/Cleaned_Data/AB_CIwq_2016_2023.csv', skip = 0, na= "Z")
View(AB.CIwq)


### on to the Gamma analysis ###
## Scale WQ parameters
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

##### Assess goodness-of-fit 
simulateResiduals(mGamma.A, n = 1000, plot = T)

##### Look at estimates
summary(mGamma.A)

## Varying one variable at a time; Station and Date = NA to make population-level predictions (i.e., predictions that ignore the random effects subjects)
newDat <- expand.grid(GroupYear = unique(AB.CIwq$GroupYear), Temperature = unique(AB.CIwq$Temperature, na.rm = T), 
                      Salinity = mean(AB.CIwq$Salinity, na.rm = T), DO = mean(AB.CIwq$DO, na.rm = T), 
                      pH = mean(AB.CIwq$pH, na.rm = T), Station = NA, Date = NA, GYD = NA) #Again, should Station be Region/Section instead?
predsGamma <- predict(mGamma.A, newdata = newDat, type = "link", se.fit = T)
preddsGamma <- data.frame(newDat, fit = predsGamma$fit, se = predsGamma$se.fit)
preddsGamma$mean <- exp(preddsGamma$fit)
preddsGamma$lwr <- exp(preddsGamma$fit - 1.96*preddsGamma$se)
preddsGamma$upr <- exp(preddsGamma$fit + 1.96*preddsGamma$se)
preddsGamma <- preddsGamma %>% separate(GroupYear, into = c("Region", "Group", "Year"), remove = FALSE) %>%
  mutate(LocationGroup = paste0(Region," ", Group))

desired_order <- c("West Section", "Central Section", "East Section")
preddsGamma$LocationGroup <- factor(preddsGamma$LocationGroup, levels = desired_order)
preddsGamma$pre.post <- ifelse(preddsGamma$Year < 2020, "pre", "post")
preddsGamma$Yearn<- as.numeric(preddsGamma$Year)

## Adding customized colors for years and include designated section colors to match the map and scatter plots above
strip <- strip_themed(background_x = elem_list_rect(fill = Section.colors))

bar.colors <- c("2016"= "#FFBBBB",
                "2017"="#56B4E9",
                "2018"="#CC79A7",
                "2019"="#FF7000",
                "2021"= "#0067A0",
                "2022"="#94DCC8",
                "2023"="grey")

(gammaPlot.a <- ggplot(preddsGamma, aes(x = Temperature, y = mean, fill = LocationGroup)) +
    geom_line() + 
    geom_ribbon(aes(ymin = lwr, ymax = upr), alpha = 0.25) + 
    facet_wrap(~Year) + 
    theme_bw() + 
    theme(axis.text = element_text(color = "black"), 
          panel.grid.major.y = element_line(colour = "grey90"), 
          panel.grid.minor.y = element_line(colour = "grey90", linetype = "dashed"), 
          panel.grid.major.x = element_blank(), legend.title = element_blank(), legend.position = "bottom") + 
    scale_fill_manual(values = Section.colors)+
    scale_y_continuous(limits = c(0, 5), 
                       breaks = seq(0, 5, 0.5), expand = expansion(add = c(0, 0))) + 
    labs(x = "Mean temperature (°C)", y = "Mean condition index (± 95% CL)", title = "Condition Index of Eastern Oysters in Apalachicola Bay 2016-2023"))


##### Summarize the raw data (useful to compare to model estimates)
c_summ <- AB.CI %>% group_by(GroupYear) %>% summarise(N = n(), mn = mean(CI, na.rm = T), se = sd(CI, na.rm = T)/sqrt(N), lwr = mn-1.96*se, upr = mn+1.96*se) %>% separate(GroupYear, into = c("Section", "Group", "Year"), remove = FALSE) %>% mutate(LocationGroup = paste0(Section," ", Group))
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

##### Compare means among years for each Section
(emmGamma <- data.frame(emmeans(mGamma.A, specs = pairwise ~ GroupYear)$emmeans)) # these are marginal or group-specific means
(conGamma <- data.frame(emmeans(mGamma.A, specs = pairwise ~ GroupYear)$contrasts)) # these are comparisons among groups

##### Break contrasts down by Section and year...the contrasts above (con) are to extensive because 
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
                       Temperature= unique(AB.CI$Temperature, na.rm = T), 
                       Salinity = c(min(AB.CI$Salinity, na.rm = T), mean(AB.CI$Salinity, na.rm = T), max(AB.CI$Salinity, na.rm = T)), 
                       DO = mean(AB.CI$DO, na.rm = T),
                       pH = mean(AB.CI$pH, na.rm = T), 
                       Station = NA, 
                       Date = NA
) 


