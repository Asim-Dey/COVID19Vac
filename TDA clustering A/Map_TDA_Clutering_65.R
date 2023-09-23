
#install.packages("usmap")
library(usmap) #import the package
library(ggplot2) #use ggplot2 to add layer for visualization
#################################################################################################

library("MazamaSpatialPlots")


data  <- read.csv("Cluster_65.csv")
head(data)

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




