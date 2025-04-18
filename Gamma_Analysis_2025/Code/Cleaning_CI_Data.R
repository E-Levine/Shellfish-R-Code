### Gamma Regression Analysis - Clean CI Data ### 
#
## Compile CI data from ODIN - Select parameters, add Dates, Stations, & Sections
## Output of cleaned data
# 
#Load require packages (install as necessary)
library(dplyr)
library(lubridate)
#
##Set Estuary working within, and start and end years of data. 
Est <- c("AB") #Enter two-letter Estuary code.
Start_year <- c("2016") #Year the dataset will begin (varies by bay)
End_year <- c("2023") #Year the dataset ends 
#
## Make sure datasets are saved in .csv format
##### Load datasets
CI <-read.csv("Gamma_Analysis_2025/Data/Raw_data/ConditionIndex.csv")
View(CI)
#
FixedLoc <- read.csv("Gamma_Analysis_2025/Data/Raw_data/FixedLocations.csv")
#
# Extract year, month, and day, add it to CI df 
CI <- CI %>% 
  mutate(
    year = substr(SampleEventID, 8, 11),
    month = substr(SampleEventID, 12, 13),
    day = substr(SampleEventID, 14, 15),
    Date = as.Date(paste(year, month, day, sep = "-")) 
  )%>%
  select(-year, -month, -day)
#
# Extracting and adding FixedLocationID to CI
substr_from_19 <- substr(CI$SampleEventID, 19, nchar(CI$SampleEventID))
pos <- regexpr("_", substr_from_19)+14
FixedLocationID <- substr(CI$SampleEventID, pos, pos+3)
CI$FixedLocationID <- FixedLocationID
#
# Merge FixedLoc to CI based on FixedLocationID to filter only collection based stations and sections
CI <- CI %>%
  left_join(
    FixedLoc %>% select(FixedLocationID, Estuary, SectionName, StationName, StationNumber), 
    by ="FixedLocationID")
#
#
#Simplify Column Names
CI <- CI %>%
  rename(
    Section = SectionName,
    Station = StationNumber)

# Tidy up the df
CI <- CI %>%  
  filter(Estuary == Est) %>%
  filter(year(Date) >= Start_year & year(Date) <= End_year)%>%
  select(OysterID, Date, Section, StationName, Station, 
         ShellHeight, ShellLength, ShellWidth, TotalWeight, 
         TarePanWeight, TissueWetWeight, ShellWetWeight, TissueDryWeight, ShellDryWeight)

#Save to .csv
write.csv(CI, 
          file = paste0("Gamma_Analysis_2025/Data/Cleaned_data/",Est,"_filtered_CI_", Start_year, "_", End_year,".csv"), 
          row.names = FALSE)

