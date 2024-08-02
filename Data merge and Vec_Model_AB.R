
library(bookdown)
library(skimr)
library(mlmRev)
library(lme4)
library(lmtest)
library(rstanarm)
library(ggplot2)
library(dplyr)



##############################################

dat1 <- read.csv("C:/Users/adey/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/VaccDat_B.csv")
dat1 <- read.csv("C:/Users/asimi/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/VaccDat_B.csv")

dim(dat1)
head(dat1)





###Combining levels of some variables################################



dat1$VACC_STATUS[dat1$VACC_STATUS==1] =0
dat1$VACC_STATUS[dat1$VACC_STATUS==2] =1



#-------------------------------------------------------------------

dat1$GENDER[dat1$GENDER==1] =0
dat1$GENDER[dat1$GENDER==2] =1



#-------------------------------------------------------------------

dat1$EDUC[dat1$EDUC==1] =0
dat1$EDUC[dat1$EDUC==2] =0
dat1$EDUC[dat1$EDUC==3] =0


dat1$EDUC[dat1$EDUC==4] =1
dat1$EDUC[dat1$EDUC==5] =1

dat1$EDUC[dat1$EDUC==6] =2
dat1$EDUC[dat1$EDUC==7] =3

###################################################################

dat1$INCOME[dat1$INCOME==1] =0
dat1$INCOME[dat1$INCOME==2] =0
  
  
  
dat1$INCOME[dat1$INCOME==3] =1
dat1$INCOME[dat1$INCOME==4] =1

dat1$INCOME[dat1$INCOME==5] =2
dat1$INCOME[dat1$INCOME==6] =2

dat1$INCOME[dat1$INCOME==7] =3
dat1$INCOME[dat1$INCOME==8] =3



##################################################################################################
## proportion #########

# Making balanced data ###########

vacP=dat1$VACC_STATUS[dat1$VACC_STATUS==0]
vacN=dat1$VACC_STATUS[dat1$VACC_STATUS==1]

length(vacP)/length(dat1$VACC_STATUS) # 0.8850733
length(vacN)/length(dat1$VACC_STATUS) # 0.1149267



dim(dat1) # 57367     9


dat_N<-dat1[dat1$VACC_STATUS==1,]
dim(dat_N) # 6593    9

dat_P<-dat1[dat1$VACC_STATUS==0,]
dim(dat_P) # 50774    9


dat_P1<-dat_P[sample(nrow(dat_P), 10000), ]
dim(dat_P1) # 1000    9




dat2<-rbind(dat_P1,dat_N)
dim(dat2)

#############################################
############################################################################################

#write.csv(dat2,"C:/Users/adey/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/VaccDat_A.csv")




dat22 <- read.csv("C:/Users/adey/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/VaccDat_A.csv")




dat22$VACC_STATUS<-as.factor(dat22$VACC_STATUS) 
dat22$GENDER <- as.factor(dat22$GENDER)
dat22$RACE <- as.factor(dat22$RACE)
dat22$EDUC <- as.factor(dat22$EDUC) 
dat22$INCOME <- as.factor(dat22$INCOME)
dat22$STATE <- as.factor(dat22$STATE)

dat22$REGION <- as.factor(dat22$REGION)


################ Model -- glmer #########################################


Vec_M1_0 <- glmer(VACC_STATUS ~ GENDER + RACE + 
                    + EDUC + INCOME + 
                    (1|STATE), family=binomial, data=dat22)


summary(Vec_M1_0)



############## Bayesian inference using rstanarm #################################

Vacc_M1 <- stan_glmer(VACC_STATUS ~ GENDER + RACE + EDUC + INCOME + 
                        (1|STATE), 
                      family = binomial("logit"), 
                      data = dat1,
                      iter = 100,
                      chains = 1,
                      seed = 349)




# display the median and the median absolute deviation (MAD) of the posterior draws
print(Vacc_M1, digits = 2)


# display the mean and standard deviation of the marginal posterior distribution of each parameter
summary(Vacc_M1, 
        probs = c(0.025, 0.975),
        digits = 2)




######## What are the priors ########

prior_summary(Vacc_M1)





