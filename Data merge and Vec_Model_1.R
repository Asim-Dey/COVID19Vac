
library(bookdown)
library(skimr)
library(mlmRev)
library(lme4)
library(lmtest)
library(rstanarm)
library(ggplot2)
library(dplyr)



##############################################

dat1 <- read.csv("C:/Users/adey/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/VaccDat.csv")
dat1 <- read.csv("C:/Users/asimi/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/VaccDat.csv")

dim(dat1)
head(dat1)





###Combining levels of some variables################################



dat1$VACC_STATUS[dat1$VACC_STATUS==1] =0
dat1$VACC_STATUS[dat1$VACC_STATUS==2] =1



dat1$EDUC[dat1$EDUC==1] =1
dat1$EDUC[dat1$EDUC==2] =1

dat1$EDUC[dat1$EDUC==3] =2
dat1$EDUC[dat1$EDUC==4] =3
dat1$EDUC[dat1$EDUC==5] =3

dat1$EDUC[dat1$EDUC==6] =4
dat1$EDUC[dat1$EDUC==7] =4

###################################################################

dat1$INCOME[dat1$INCOME==1] =1

dat1$INCOME[dat1$INCOME==2] =2
dat1$INCOME[dat1$INCOME==3] =2

dat1$INCOME[dat1$INCOME==4] =3
dat1$INCOME[dat1$INCOME==5] =3
dat1$INCOME[dat1$INCOME==6] =3
dat1$INCOME[dat1$INCOME==7] =3
dat1$INCOME[dat1$INCOME==8] =3

















dat1$VACC_STATUS<-as.factor(dat1$VACC_STATUS) 

dat1$GENDER <- as.factor(dat1$GENDER)
dat1$RACE <- as.factor(dat1$RACE)
dat1$HISPANIC <-as.factor(dat1$HISPANIC)
dat1$EDUC <- as.factor(dat1$EDUC) 
dat1$INCOME <- as.factor(dat1$INCOME)
dat1$STATE <- as.factor(dat1$STATE)

dat1$REGION <- as.factor(dat1$REGION)


################ inference using glmer #########################################
Vec_M1_0 <- glmer(VACC_STATUS ~ GENDER + RACE + 
              + EDUC + INCOME + 
              (1+|STATE), family=binomial, data=dat1)



############## Bayesian inference using rstanarm #################################

Vacc_M1 <- stan_glmer(VACC_STATUS ~ GENDER + RACE + EDUC + INCOME + 
                             (1+REGION|STATE), 
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





