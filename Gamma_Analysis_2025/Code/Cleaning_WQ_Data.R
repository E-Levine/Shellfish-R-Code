### Clean WQ Data ###
#
## Compile WQ data from ODIN - Select parameters, add Dates, Stations, & Sections
## Output of cleaned data
# 
#Load require packages (install as necessary)
library(dplyr)
library(lubridate)
#
##Set Estuary working within, and start and end years of data. 
Est <- c("TB") #Enter two-letter Estuary code.
Start_year <- c("2005") #Year the dataset will begin (varies by bay)
End_year <- c("2024") #Year the dataset ends 
SampleEv <- c("RCRT") #4 letter Sample Event Code
#
## Make sure datasets are saved in .csv format
##### Load datasets
WQ <-read.csv("Raw_data/SampleEventWQ.csv")
View(WQ)
#
FixedLoc <- read.csv("Raw_data/FixedLocations.csv")
#
# Extract year, month, and day, add it to WQ df 
WQ <- WQ %>% 
  mutate(
    year = substr(SampleEventID, 8, 11),
    month = substr(SampleEventID, 12, 13),
    day = substr(SampleEventID, 14, 15),
    Date = as.Date(paste(year, month, day, sep = "-")) 
  )%>%
  select(-year, -month, -day)
#
# Extracting and adding FixedLocationID to WQ
substr_from_19 <- substr(WQ$SampleEventID, 19, nchar(WQ$SampleEventID))
pos <- regexpr("_", substr_from_19)+14
FixedLocationID <- substr(WQ$SampleEventID, pos, pos+3)
WQ$FixedLocationID <- FixedLocationID
#
# Extracting and adding SampleEvent to WQ
substr_from_3 <- substr(WQ$SampleEventID, 3, nchar(WQ$SampleEventID))
pos2 <- regexpr("AB", substr_from_3)+4
SampleEvent <- substr(WQ$SampleEventID, pos2, pos2+3)
WQ$SampleEvent <- SampleEvent


#
# Merge FixedLoc to WQ based on FixedLocationID to filter only collection based stations and sections
WQ <- WQ %>%
  left_join(
    FixedLoc %>% select(FixedLocationID, Estuary, SectionName, StationName, StationNumber), 
    by ="FixedLocationID")
#
#
#Simplify Column Names
WQ <- WQ %>%
  rename(
    Section = SectionName,
    Station = StationNumber)

# Tidy up the df
WQ <- WQ %>%  
  filter(Estuary == Est) %>%
  filter(year(Date) >= Start_year & year(Date) <= End_year)%>%
  filter(SampleEvent == SampleEv)%>%
  select(Estuary, SampleEvent, Date, Section, StationName, Station, 
         Temperature, Salinity, DissolvedOxygen, pH, 
         Depth, SampleDepth, Secchi, TurbidityYSI, TurbidityHach)

#Save to .csv
write.csv(WQ, 
          file = paste0("Clean_data/",Est,"_filtered_WQ_", SampleEv, "_", Start_year, "_", End_year,".csv"), 
          row.names = FALSE)
