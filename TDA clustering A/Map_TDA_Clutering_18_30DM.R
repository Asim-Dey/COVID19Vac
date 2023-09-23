



#install.packages("usmap")
library(usmap) #import the package
library(ggplot2) #use ggplot2 to add layer for visualization

plot_usmap(regions = "counties") + 
  labs(title = "U.S. counties",
       subtitle = "This is a blank map of the United States.") + 
  theme(panel.background=element_blank())



#################################################################################################

library("MazamaSpatialPlots")


data  <- read.csv("C:/Users/asimi/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/cluster_member_TDA_18_30DM.csv")
head(data)
dim(data)

# Added required 'stateCode' and 'countyFIPS' variables
data <- 
  data %>%
  dplyr::mutate(
    stateCode = stateCode,
    countyFIPS = MazamaSpatialUtils::US_countyNameToFIPS(stateCode, county_name),
    pUninsured2010 = Cluster_member,
    .keep = "none"
  ) 



# Create map
countyMap(
  data = data,
  parameter = 'pUninsured2010',
  # palette = "Set1",
  palette= c("green","blue4", "orange", "red","gold"),
  legendTitle = 'Cluster',
  title = ""
)








