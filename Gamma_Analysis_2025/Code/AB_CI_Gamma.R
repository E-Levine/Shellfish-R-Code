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
desired_order <- c("West Section", "Central Section", "East Section")
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
means <- aggregate(CI ~ AB.CIwq$GroupYear, data = AB.CIwq, FUN= mean, levels = desired_order)
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


