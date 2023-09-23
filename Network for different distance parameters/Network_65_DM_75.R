#install.packages("igraph")
library(igraph) 

library(ggplot2)


D0<-read.csv('County_VaccinationsAA.csv')

names(D0)
head(D0) 

length(unique(D0$STATE))




boxplot1 = ggplot(D0, aes(x = state, y = Series_Complete_65PlusPop_Pct)) + geom_boxplot()+ 
  #coord_cartesian(ylim = c(0,30))+
  coord_flip()


boxplot1

#==========================
boxplot1 = ggplot(D0, aes(x = state, y = covid_cases_per_100k)) + geom_boxplot()+ 
  #coord_cartesian(ylim = c(0,30))+
  coord_flip()


boxplot1




###################### Distance matrix #############################

#================= Vec Distance matrix ===================

#Series_Complete_65PlusPop_Pct

vec65<-D0$Series_Complete_65PlusPop_Pct

summary(vec65)

#Min. 1st Qu.  Median    Mean  3rd Qu.    Max. 
#17.60   73.70   83.10   81.06   91.00   95.00 

#vec18[vec18>70]
#length(vec18)

DD1<-as.matrix(vec65)
n<-dim(DD1)[1]
n

##===============  Geo Distance matrix =======================

library(geosphere)

#fun = distGeo/distHaversine

M_Goe<-distm(D0[,c(7,8)],fun = distGeo)/ 1609.344 # Distance in miles
dim(M_Goe)

summary(as.numeric(M_Goe))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0   476.0   766.2   851.7  1128.3  2858.1 
#########################################################################

k2=62.46   # Mean or above of 18 netowrk for comparison

#k1=0.05 # Counties which vec rate similar dist less than k2 (less than median)


# Counties which are within "DM" miles are connected by an edge
DM=75 # 30, 75


A1 <- matrix(0, nrow = n, ncol =n);dim(A1)
A2 <- matrix(0, nrow = n, ncol =n);dim(A2)

############### Step 1 #####################################

for (i in 1:n){ #i=1
  for (j in 1:n){ #j=2
    
    
    if  (DD1[,1][i] > k2 & DD1[,1][j]> k2)  A2[i,j]=1
    
    if  (DD1[,1][i] > k2 & DD1[,1][j]> k2 & M_Goe[i,j] < DM)  A1[i,j]=1
    
    
  }
  
  print(i)
}



dim(A1)

dim(A2)


### AM- Adjacency Matrix  for nodetvec ################
#write.csv(A2,'AM_65_Dist_V2_DM75.csv')








###################### Network ###################



G0<- graph.adjacency(A1, mode="undirected",weighted=TRUE,diag=FALSE)

V0<-length(V(G0));V0 # 3192
E0<-length(E(G0));E0 #  37726

plot(G0, vertex.size=5, 
     vertex.label=NA,vertex.color="red",edge.arrow.size=0.9)





edlist_65_Dist<-cbind(get.edgelist(G0))


#edlist_65<-cbind(get.edgelist(G1))


#write.csv(edlist_65_Dist,'edlist_65_Dist_V2_aboveMean_DM75.csv')



############################################################################################
################### Statistical Analysis #############################




##########################################################################################

#set_vertex_attr(G1, 'vac18', index = V(G1), value=D0$Series_Complete_18PlusPop_Pct)



G11<-graph.adjacency(A1, mode="undirected",weighted=TRUE,diag=FALSE)%>% 
  set_vertex_attr("vac18", value = D0$Series_Complete_18PlusPop_Pct)%>% 
  set_vertex_attr("state", value = as.character(D0$state))%>% 
  set_vertex_attr("region", value = as.character(D0$Region_N))%>% 
  set_vertex_attr("Edu_ass_coll", value = D0$Percent.of.adults.completing.some.college.or.associate.s.degree..2016.20)%>% 
  set_vertex_attr("Income", value = D0$Median_Household_Income_2020)%>% 
  set_vertex_attr("Unemployment", value = D0$Unemployment_rate_2021)%>% 
  set_vertex_attr("Covid_cases", value = D0$covid_cases_per_100k)


Atr<-get.vertex.attribute(G11)
Atr



############ Average nearest neighbor degree #################

ND0<-knn(G11)
ND0

ND<-ND0$knn
ND[is.na(ND)] <- 0 # replace NA by 0
ND


vac18<-Atr$vac18
state1<-Atr$state
region1<-Atr$region
d0a<-data.frame(ND,state1,region1)




boxplot1 = ggplot(d0a, aes(x = state1, y = ND)) + geom_boxplot()+ 
  #coord_cartesian(ylim = c(0,30))+
  coord_flip()+ylab("Average nearest neighbor degree") +
  xlab("US State")

boxplot1



boxplot2 = ggplot(d0a, aes(x = region1, y = ND, fill=region1)) + geom_boxplot()+ 
  #coord_cartesian(ylim = c(0,30))+
  coord_flip()+ylab("Average nearest neighbor degree") +
  xlab("US Census Region")+
  theme(legend.position="none")

boxplot2

################# Degree distribution #############



degree1<-degree(G0)
summary(degree1)
table(degree1)


meanD1<-mean(degree1);meanD1 # 23.63784
sdD1<-sd(degree1);sdD1 # 15.77914


hist(degree1,xlab="Degree",col="blue",main="",first.panel=grid())
box()

summary(ND)

plot(degree1, ND, log="xy",
     col="goldenrod", xlab="Log node Degree",
     ylab="Log Average Neighbor Degree",first.panel=grid())








degree_D65<-degree.distribution(G0)

plot(degree_D65, type="o", lwd=2,ylab = "p(k)",xlab = "Number of edge (Degree), k",
     main = "",xlim = c(),col="blue",first.panel=grid())










########## APL and Diameter ###########################

AVPL<-average.path.length(G0);AVPL # 14.10547
Diameter<-diameter(G0) ;Diameter # 61





#################### Assortativity and homophily ###########################

assortativity_degree(G11)  #  0.8973962

assortativity(G11, V(G11)$vac18, directed=F) #  0.4003912


assortativity_nominal(G11, as.integer(as.factor(V(G11)$state), directed=F)) #  0.7019899


assortativity(G11, V(G11)$Covid_cases, directed=F) # 0.3121925


assortativity(G11, V(G11)$Edu_ass_coll, directed=F) # 0.4370742

assortativity(G11, V(G11)$Income, directed=F) # 0.4692469






