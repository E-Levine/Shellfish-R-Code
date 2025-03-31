## CI and WQ data organization and plots  ##
#
# This code is to prepare the cleaned data sets for gamma analysis and create plots used for reports
#
library(Matrix)
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

Section.colors <- c("West Section" = "#00a884", 
                    "Central Section" = "#FF7979", 
                    "East Section" = "#6D9DFF")

AB.CIwq$Group <- factor(AB.CIwq$Group, levels = c("West Section", "Central Section", "East Section"))
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

### export AB.CIwq as new .csv and move on to new script for only analysis
write.csv(AB.CIwq, 
          file = paste0("Data/Cleaned_data/","AB_CIwq_2016_2023.csv"), 
          row.names = FALSE)

