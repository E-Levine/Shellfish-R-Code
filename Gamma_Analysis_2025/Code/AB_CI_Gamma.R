## CI Gamma Analysis  ##
#
# This code is to run the gamma regressionmodel using the datasets previously prepped
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