

#install.packages("igraph")
library(igraph) 

library(ggplot2)


D0<-read.csv('County_VaccinationsAA.csv')

names(D0)
head(D0) 



vec18<-D0$Series_Complete_18PlusPop_Pct
vec65<-D0$Series_Complete_65PlusPop_Pct

############### Clusters N 18 ###########################################

D18<-data.frame(D0$Series_Complete_18PlusPop_Pct,D0$Cluster_member_18)


cls1<-D18[D18$D0.Cluster_member_18==1,]$D0.Series_Complete_18PlusPop_Pct
cls2<-D18[D18$D0.Cluster_member_18==2,]$D0.Series_Complete_18PlusPop_Pct
cls3<-D18[D18$D0.Cluster_member_18==3,]$D0.Series_Complete_18PlusPop_Pct
cls4<-D18[D18$D0.Cluster_member_18==4,]$D0.Series_Complete_18PlusPop_Pct

#---- Cluster 1 ----------------
length(cls1)
mean(cls1)
sd(cls1)

#---- Cluster 2 ----------------
length(cls2)
mean(cls2)
sd(cls2)
#---- Cluster 3 ----------------
length(cls3)
mean(cls3)
sd(cls3)
#---- Cluster 4 ----------------
length(cls4)
mean(cls4)
sd(cls4)


############### Clusters N 65 ###########################################

D65<-data.frame(D0$Series_Complete_65PlusPop_Pct,D0$Cluster_member_65)


cls651<-D65[D65$D0.Cluster_member_65==1,]$D0.Series_Complete_65PlusPop_Pct
cls652<-D65[D65$D0.Cluster_member_65==2,]$D0.Series_Complete_65PlusPop_Pct
cls653<-D65[D65$D0.Cluster_member_65==3,]$D0.Series_Complete_65PlusPop_Pct
cls654<-D65[D65$D0.Cluster_member_65==4,]$D0.Series_Complete_65PlusPop_Pct

#---- Cluster 1 ----------------
length(cls651)
mean(cls651)
sd(cls651)

#---- Cluster 2 ----------------
length(cls652)
mean(cls652)
sd(cls652)
#---- Cluster 3 ----------------
length(cls3)
mean(cls3)
sd(cls3)
#---- Cluster 4 ----------------
length(cls654)
mean(cls654)
sd(cls654)











