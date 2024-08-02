#install.packages("igraph")

library(igraph) 
library(ggplot2)
library(sand)
library(ergm)
library(statnet)

 

D0<-read.csv('~/Library/CloudStorage/OneDrive-Personal/COVID-19 Vac Network/Data and Code/Data/COVID-19 Vaccinations in the US/County_VaccinationsAA.csv')
#D0<-read.csv('C:/Users/adey/OneDrive/COVID-19 Vac Network/Data and Code/Data/COVID-19 Vaccinations in the US/County_VaccinationsAA.csv')


#D0<-read.csv('C:/Users/asimi/OneDrive/COVID-19 Vac Network/Data and Code/Data/COVID-19 Vaccinations in the US/County_VaccinationsAA.csv')

names(D0)
head(D0) 

length(unique(D0$STATE))



################ 18 #####################################


#A0<-read.csv('C:/Users/adey/OneDrive/COVID-19 Vac Network/Data and Code/Data/COVID-19 Vaccinations in the US/edlist_18_Dist_V2_aboveMean.csv')

A0<-read.csv('~/Library/CloudStorage/OneDrive-Personal/COVID-19 Vac Network/Data and Code/Data/COVID-19 Vaccinations in the US/edlist_18_Dist_V2_aboveMean.csv')


G01<-graph_from_edgelist(as.matrix(A0), directed = FALSE)
#plot(G01)

V01<-length(V(G01));V01 # 
E01<-length(E(G01));E01 #  

A1B<-as_adjacency_matrix(G01)
dim(A1B) # 3191 3191

########################################################################################
#install.packages('ergm')
library(sand)
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

###############################################################

G18_ergm <- formula(G18 ~ edges                    # AIC: 90129  BIC: 90142  (Smaller is better. MC Std. Err. = 0.05086)
                    + nodemain("Income")           # IC: 86741  BIC: 86768  (Smaller is better. MC Std. Err. = 0.008115)
                    + nodemain("Covid_cases")      # AIC: 86316  BIC: 86356  (Smaller is better. MC Std. Err. = 0.01157)
                    + nodemain("Edu_coll_higher")  # AIC: 85896  BIC: 85950  (Smaller is better. MC Std. Err. = 0.01053)
                    + nodefactor("region"))        # AIC: 84248  BIC: 84342  (Smaller is better. MC Std. Err. = 0.02465)
                                                 


  
set.seed(1234)
G18_ergm.fit <- ergm(G18_ergm,
                     control = control.ergm(seed=135,force.main=TRUE,MCMC.samplesize=1000))

G18_ergm.fit
summary(G18_ergm.fit)

anova(G18_ergm.fit)

################## mcmc.diagnostics #################################

mcmc.diagnostics(G18_ergm.fit)


gof1<-gof(G18_ergm.fit, GOF=~model)
plot(gof1)

gof1<-gof(G18_ergm.fit, plotlogodds=TRUE)
plot(gof1)









########################################################################################
################################ Including endogeneous variable ###############################

G18_ergm_B <- formula(G18 ~ edges
                     + kstar (2)
                     + triangle
)
                     + nodemain("Income")
                     + nodemain("Covid_cases")
                     + nodemain("Edu_coll_higher")
                     + nodefactor("region")
                      
)

summary(G18_ergm_B)


set.seed(1234)
G18_ergm_B.fit <- ergm(G18_ergm_B, 
                       control = control.ergm(seed=132,
                                   main.method = "MCMLE",
                                   force.main=TRUE,
                                   MCMLE.maxit = 20,
                                   MCMC.samplesize = 5000
                                   )
                     
)


                       




G18_ergm_B.fit
 
  

################## mcmc.diagnostics #################################

mcmc.diagnostics(G18_ergm_B.fit)

gof2<-gof(G18_ergm_B.fit, GOF=~model)
plot(gof2)














########################################################################################
######################### Bayesian ERGM #############################

#*** Convergence inssue 
#* does not work 




#mean.priors <- list(c(-1, 0, 0), c(-1, 0, 0), c(-1, 0, 0), c(-1, 0))

#sigma <- 5
#sigma.priors <- list(diag(sigma, 3), diag(sigma, 3), diag(sigma, 3), diag(sigma, 2))

#set.seed(27)
#mod.sel <- bergm(m1, iters = 1000, mean.priors = mean.priors, sigma.priors = sigma.priors)
#summary(mod.sel)
#set.seed(27)
#bgof(mod.sel, aux.iters = 10000, n.deg = 20, n.dist = 10, n.esp = 15) 


s_Covid<-var(D0$covid_cases_per_100k)
s_Income<-var(D0$Median_Household_Income_2020)
s_Education<-var(D0$Percent.of.adults.with.a.bachelor.s.degree.or.higher.2015.19)

dim2=7
mean.priors <- rep(0,7)

S<-c(100,s_Income,s_Covid,s_Education,100,100,100)
sigma.priors <- diag(S, dim2)
dim(sigma.priors)



library("Bergm")


m11 <- formula(G18 ~ edges
               + nodemain("Income")
               + nodemain("Covid_cases")
               + nodemain("Edu_coll_higher")
               
               + nodefactor("region")
)


set.seed(1127)
mod.sel2 <- bergm(m11, iters = 20000, mean.priors = mean.priors, sigma.priors = sigma.priors)

summary(mod.sel2)

plot(mod.sel2)



set.seed(27)
bgof(mod.sel2, aux.iters = 20000, n.deg = 20, n.dist = 10, n.esp = 15) 










