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
#
## Need to rerun the mGamma code again. Because mGamma.A had lowest AIC score, we will only use that one moving forward for comparisons
## Scale WQ parameters
str(AB.CIwq)
## run following code if any are not numeric
#AB.CIwq$CI <- as.numeric(AB.CIwq$CI)
#AB.CIwq$Temperature <- as.numeric(AB.CIwq$Temperature)
#AB.CIwq$Salinity <- as.numeric(AB.CIwq$Salinity)
#AB.CIwq$DO <- as.numeric(AB.CIwq$DO)
#AB.CIwq$pH <- as.numeric(AB.CIwq$pH)

scDO <- scale(AB.CIwq$DO)
scTemp <- scale(AB.CIwq$Temperature)
scSalinity <- scale(AB.CIwq$Salinity)
scpH <- scale(AB.CIwq$pH)
#
#
table(is.na(AB.CIwq))
## Run following code if NAs are present. most functions will remove them automatically but good to check 
#AB.CIwq <- AB.CIwq %>% filter(complete.cases(.))

##### Fit gamma regressions with the fixed effect all wq parameters 
AB.CIwq$GYD <- droplevels(interaction(AB.CIwq$GroupYear, AB.CIwq$Date))
mGamma.A <- glmmTMB(CI ~ Temperature + Salinity + DO + pH + GroupYear + (1|GYD), 
                    dispformula = ~ GroupYear,
                    family = Gamma(link = "log"), 
                    data = AB.CIwq)

##### Compare means among years for each Section
(emmGamma <- data.frame(emmeans(mGamma.A, specs = pairwise ~ GroupYear)$emmeans)) # these are marginal or group-specific means
(conGamma <- data.frame(emmeans(mGamma.A, specs = pairwise ~ GroupYear)$contrasts)) # these are comparisons among groups

##### Break contrasts down by Section and year...the contrasts above (con) are to extensive because 
##### there are contrasts we don't care about, so this code reduces the number of comparisons and then 
##### adjusts p-values accordingly for multiple comparisons. Here, we are only comparing among years for
##### each group (west, central, and east)
(Gamma_contrasts <- conGamma %>% separate(contrast, into = c("contrast1", "contrast2"), sep = "\\-", extra = "merge") %>% 
    separate(contrast1, into = c("Group1", "Year1"), sep = "\\.", extra = "merge") %>% mutate(contrast2 = trimws(contrast2)) %>% 
    separate(contrast2, into = c("Group2", "Year2"), sep = "\\.", extra = "merge") %>% 
    mutate(keep = case_when(Group1==Group2 ~ "Keep", .default = "Drop")) %>% filter(keep == "Keep") %>%
    mutate(p_adj = p.adjust(p.value, method = "fdr")) %>% arrange(Group1, Group2) %>% 
    mutate(p.value = round(p.value, 3), p_adj = round(p_adj, 3)))
# Convert Group1 and Group2 to factors with desired order
desired_order <- c("West Section", "Central Section", "East Section")
(Gamma_contrasts <- Gamma_contrasts %>%
    mutate(Group1 = factor(Group1, levels = desired_order),
           Group2 = factor(Group2, levels = desired_order)) %>%
    arrange(Group1, Group2))

##### Filter for statistically important (not necessarily biologically important) pairwise comparisons
(gamma_important_contrasts <- Gamma_contrasts[Gamma_contrasts$p_adj<0.05,])
## nothing significant so zero rows. 
#@colin, does this mean there is no significant effect on CI by WQ over the years? or that there was no significant change of WQ over the years?

# Save contrast outputs as tab-delimited text
#write.table(gamma_important_contrasts, "Gamma_contrasts.txt", sep = "\t", row.names = FALSE)

# CI contrasts 
mGamma.B <- glmmTMB(CI ~ GroupYear, family = Gamma(link = "log"), data = AB.CIwq)

newDat.CI <- expand.grid(GroupYear = unique(AB.CIwq$GroupYear), Section = NA)
predsGamma.CI <- predict(mGamma.B, newdata = newDat.CI, type = "link", se.fit = T)
preddsGamma.CI <- data.frame(newDat.CI, fit = predsGamma.CI$fit, se = predsGamma.CI$se.fit)
preddsGamma.CI$mean <- exp(preddsGamma.CI$fit)
preddsGamma.CI$lwr <- exp(preddsGamma.CI$fit - 1.96*preddsGamma.CI$se)
preddsGamma.CI$upr <- exp(preddsGamma.CI$fit + 1.96*preddsGamma.CI$se)
preddsGamma.CI <- preddsGamma.CI %>% separate(GroupYear, into = c("Section", "Group", "Year"), remove = FALSE) %>%
  mutate(LocationGroup = paste0(Section," ", Group))
preddsGamma.CI$LocationGroup <- factor(preddsGamma.CI$LocationGroup, levels = desired_order)
summary(preddsGamma.CI)
str(preddsGamma.CI)

preddsGamma.CI$pre.post <- ifelse(preddsGamma.CI$Year < 2020, "pre", "post")
preddsGamma.CI$Yearn<- as.numeric(preddsGamma.CI$Year)

##### Compare means among years for each region
(emmGamma.CI <- data.frame(emmeans(mGamma.B, specs = pairwise ~ GroupYear)$emmeans)) # these are marginal or group-specific means
(conGamma.CI <- data.frame(emmeans(mGamma.B, specs = pairwise ~ GroupYear)$contrasts)) # these are comparisons among groups

##### Break contrasts down by region and year...the contrasts above (con) are to extensive because there are contrasts we don't care about, so this code reduced the number of comparisons and then adjusts p-values accordingly for multiple comparisons. Here, we are only comparing among years for each group (west, central, and east)
(Gamma_contrasts.CI <- conGamma.CI %>% separate(contrast, into = c("contrast1", "contrast2"), sep = "\\-", extra = "merge") %>% 
    separate(contrast1, into = c("Group1", "Year1"), sep = "\\.", extra = "merge") %>% mutate(contrast2 = trimws(contrast2)) %>% 
    separate(contrast2, into = c("Group2", "Year2"), sep = "\\.", extra = "merge") %>% 
    mutate(keep = case_when(Group1==Group2 ~ "Keep", .default = "Drop")) %>% filter(keep == "Keep") %>%
    mutate(p_adj = p.adjust(p.value, method = "fdr")) %>% arrange(Group1, Group2) %>% 
    mutate(p.value = round(p.value, 3), p_adj = round(p_adj, 3)))

(Gamma_contrasts.CI <- Gamma_contrasts.CI %>%
    mutate(Group1 = factor(Group1, levels = desired_order),
           Group2 = factor(Group2, levels = desired_order)) %>%
    arrange(Group1, Group2))

##### Filter for statistically important (not necessarily biologically important) pairwise comparisons
(gamma_important_contrasts_CI <- Gamma_contrasts.CI[Gamma_contrasts.CI$p_adj<0.05,])
write.table(gamma_important_contrasts_CI, "contrasts.CI.txt", sep = "\t", row.names = FALSE)

##Make a table for CI means
means.CI <- aggregate(CI ~ GroupYear, data = AB.CIwq, FUN = mean)

# Create separate Group and Year columns in the output table only
means.CI$Group <- sub("\\..*", "", means.CI$GroupYear)  # Extract before first period
means.CI$Year  <- sub(".*\\.", "", means.CI$GroupYear)  # Extract after last period


# Reorder columns if needed
means.CI <- means.CI[, c("Group", "Year", "CI")]
means.CI$Group <- factor(means.CI$Group, levels = desired_order, ordered = TRUE)

# Arrange the data frame by Year (following desired_order) and then by Group if needed
means.CI <- means.CI[order(means.CI$Group, means.CI$Year), ]
means.CI

write.table(means.CI, "means.CI.txt", sep = "\t", row.names = FALSE)

