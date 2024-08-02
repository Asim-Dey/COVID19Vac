#install.packages("igraph")

library(igraph) 
library(ggplot2)


D0<-read.csv('~/Library/CloudStorage/OneDrive-Personal/COVID-19 Vac Network//Data and Code/Data/COVID-19 Vaccinations in the US/County_VaccinationsAA.csv')
#D0<-read.csv('C:/Users/adey/OneDrive/COVID-19 Vac Network//Data and Code/Data/COVID-19 Vaccinations in the US/County_VaccinationsAA.csv')
#D0<-read.csv('C:/Users/asimi/OneDrive/COVID-19 Vac Network//Data and Code/Data/COVID-19 Vaccinations in the US/County_VaccinationsAA.csv')

names(D0)
head(D0) 

length(unique(D0$STATE))



################ 18 #####################################

A0<-read.csv("~/Library/CloudStorage/OneDrive-Personal/COVID-19 Vac Network/Data and Code/Data/COVID-19 Vaccinations in the US/edlist_65_Dist_V2_aboveMean.csv")


G01<-graph_from_edgelist(as.matrix(A0), directed = FALSE)
#plot(G01)

V01<-length(V(G01));V01 # 3193
E01<-length(E(G01));E01 #  5794

A1B<-as_adjacency_matrix(G01)
dim(A1B) # 3191 3191

########################################################################################
#install.packages('ergm')
library(sand)
library(ergm)
library(statnet)


library(ergm)



A1C<-get.adjacency(G01)


G18<- network::as.network(as.matrix(A1C), directed=FALSE)

#network::set.vertex.attribute(G18, "vac18", D0$Series_Complete_18PlusPop_Pct)
network::set.vertex.attribute(G18, "region", D0$Region)
network::set.vertex.attribute(G18, "Covid_cases", D0$covid_cases_per_100k)

#network::set.vertex.attribute(G18, "Edu_ass_coll", D0$Percent.of.adults.completing.some.college.or.associate.s.degree..2016.20)
network::set.vertex.attribute(G18, "Edu_coll_higher", D0$Percent.of.adults.with.a.bachelor.s.degree.or.higher.2015.19)

network::set.vertex.attribute(G18, "Income", D0$Median_Household_Income_2020)
#network::set.vertex.attribute(G18, "Unemployment", D0$Unemployment_rate_2021)







###### need to change the state in D0 from character to numeric factor ###########

G56_ergm <- formula(G18 ~ edges                     # AIC: 224893  BIC: 224906  (Smaller is better. MC Std. Err. = 0.01524)
                    + nodemain("Income")           # AIC: 224495  BIC: 224522  (Smaller is better. MC Std. Err. = 0.003345)
                    + nodemain("Covid_cases")      # AIC: 224328  BIC: 224368  (Smaller is better. MC Std. Err. = 0.01086)
                    + nodemain("Edu_coll_higher")  # AIC: 224322  BIC: 224375  (Smaller is better. MC Std. Err. = 0.03218)
                    
                    + nodefactor("region"))        # AIC: 219567  BIC: 219661  (Smaller is better. MC Std. Err. = 0.02516)


 
  

set.seed(1234)
G56_ergm <- ergm(G56_ergm,
                     control = control.ergm(seed=135,force.main=TRUE,MCMC.samplesize=1000)
)



summary(G56_ergm)

anova(G56_ergm)





################## mcmc.diagnostics #################################

mcmc.diagnostics(G56_ergm)


gof1<-gof(G56_ergm, GOF=~model)

plot(gof1)










