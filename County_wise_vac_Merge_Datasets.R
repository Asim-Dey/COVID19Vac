


#D1<- read.csv("C:/Users/adey/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/COVID-19_Vaccinations_in_the_United_States_County.csv")

#D1<- read.csv("C:/Users/asimi/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/COVID-19_Vaccinations_in_the_United_States_County.csv")
#head(D1) # 
#dim(D1) # 666004     22

D2<-D1[D1$Date %in% c('12/28/2022'),]
head(D2) # 
dim(D2) # 

#unique(D2$Recip_State)
#length(unique(D2$Recip_State))



names(D2)


D21<-D2[!(D2$Recip_State %in% c("GU","HI","MP","VI","PR")),]
dim(D21)

unique(D21$Recip_State)
#length(unique(D21$Recip_State))

#write.csv(D21, 'C:/Users/adey/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/County_VaccinationsAA_Dec28_22.csv')



######################################################################################################################

D3<-read.csv('C:/Users/adey/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/County_VaccinationsAA_Dec28_22.csv')

#D3<-read.csv('C:/Users/asimi/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/County_VaccinationsAA_Dec28_22.csv')



county <- read.csv("C:/Users/adey/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/county.csv")
county <- read.csv("C:/Users/asimi/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/county.csv")
head(county)
dim(county) #3331    5


length(unique(county$STATE))


########### Merge #############################################################################################

D0 <- merge(x=D3, y=county, by="FIPS")
dim(D0) # 3322   26

length(unique(D0$STATE))


Da<-read.csv('C:/Users/adey/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/United_States_COVID-19_CasesLevels_by_County.csv')
names(Da)
dim(Da)


Da2<-Da[Da$date_updated %in% c('9/29/2022'),] # data 3 months back---sep
head(Da2) # 
dim(Da2) # 


 


DD0 <- merge(x=D0, y=Da2, by="FIPS")
dim(DD0) # 3204   24


#write.csv(DD0, 'C:/Users/adey/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/County_VaccinationsAA_0.csv')
#write.csv(D0, 'C:/Users/asimi/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/County_VaccinationsAA_0.csv')

names(DD0)

DD0a<-read.csv('C:/Users/adey/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/County_VaccinationsAA_B.csv')
head(DD0a)
names(DD0a)



unique(DD0a$state)
length(unique(DD0a$state))


#Region
#1) Northeast
#2) South
#3) Midwest
#4) West


######################## Education and Income county level ###############


D_Inc<- read.csv("C:/Users/adey/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/Unemployment_Income_county_level_US_DEPARTMENT_OF_AGRICULTURE.csv")
names(D_Inc)

D_Ed<- read.csv("C:/Users/adey/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/Education_county_lelvel_US_DEPARTMENT_OF_AGRICULTURE.csv")
names(D_Ed)


D_InEd <- merge(x=D_Ed, y=D_Inc, by="FIPS")
dim(D_InEd)

head(D_InEd)
names(D_InEd)

#########################################################################################################################
D0_AB <- merge(x=DD0a, y=D_InEd, by="FIPS")
dim(D0_AB)

head(D0_AB)
names(D0_AB)




#write.csv(D0_AB, 'C:/Users/adey/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/County_VaccinationsAA_2.csv')
#write.csv(D0, 'C:/Users/asimi/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/County_VaccinationsAA_2.csv')




