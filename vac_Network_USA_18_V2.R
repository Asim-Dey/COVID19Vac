#install.packages("igraph")
library(igraph) 

library(ggplot2)


D0<-read.csv('C:/Users/adey/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/County_VaccinationsAA.csv')

#D0<-read.csv('C:/Users/asimi/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/County_VaccinationsAA.csv')

names(D0)
head(D0) 

length(unique(D0$STATE))


 

boxplot1 = ggplot(D0, aes(x = state, y = Series_Complete_18PlusPop_Pct)) + geom_boxplot()+ 
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

#Series_Complete_18PlusPop_Pct
vec18<-D0$Series_Complete_18PlusPop_Pct



summary(vec18)
#Min. 1st Qu.  Median    Mean  3rd Qu.    Max. 
#13.90   53.30   61.80   62.46   70.80   95.00 

#vec18[vec18>70]
#length(vec18)

DD1<-as.matrix(vec18)
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

k2=62.46   # Mean or above
#k1=0.05 # Counties which vec rate similar dist less than k2 (less than median)
DM=50  # Counties which are within "DM" miles are connected by an edge


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

###################### Network ###################



G0<- graph.adjacency(A1, mode="undirected",weighted=TRUE,diag=FALSE)

V0<-length(V(G0));V0 # 3192
E0<-length(E(G0));E0 #  5794

plot(G0, vertex.size=5, 
     vertex.label=NA,vertex.color="red",edge.arrow.size=0.9)

 
# The proportion of present edges from all possible edges in the network.


edge_density(G0, loops=F) # 0.001137677



# G1<- graph.adjacency(A2, mode="undirected",weighted=TRUE,diag=FALSE)
# V<-length(V(G1));V # 3193
# E<-length(E(G1));E # 1168156


# plot(G1, vertex.size=5, 
#     vertex.label=NA,vertex.color="red",edge.arrow.size=0.9)




edlist_18_Dist<-cbind(get.edgelist(G0))

# edlist_18<-cbind(get.edgelist(G1))
# write.csv(edlist_18,'C:/Users/asimi/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/edlist_18_V2.csv')
# write.csv(edlist_18_Dist,'C:/Users/adey/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/edlist_18_Dist_V2_aboveMean.csv')


############################## Motif ###############################################
#---------------------------------------

node<-length(V(G0));node
edge<-length(E(G0));edge


#-------- Motif, size=3 ---------------

m01<-motifs(G0, 3)
m01[is.na(m01)] <- 0
m01

T2<-m01[4]


n21<-count_motifs(G0,3);n21

# Concentration
C_T1<-m01[3]/n21;C_T1    # 
C_T2<-m01[4]/n21;C_T2     # 

#drawing 3-nodes motifs method ---
par(mfrow=c(1,4))

for(i in 0:3){ 
  motifgraph <- graph.isocreate(size=3, number=i, directed=F)
  plot(motifgraph)
}

#-------- Motif, size=4 ---------------
par(mfrow=c(2,6))

for(i in 0:10){
  motifgraph <- graph.isocreate(size=4, number=i, directed=F)
  plot(motifgraph)
}


m02<-motifs(G0, 4)
m02[is.na(m02)] <- 0
m02

n22<-count_motifs(G0, 4);n22

V1_0<-m02[5];C_V1_0<-m02[5]/n22; C_V1_0[is.na(C_V1_0)] <-0
V2_0<-m02[7];C_V2_0<-m02[7]/n22; C_V2_0[is.na(C_V2_0)] <-0
V3_0<-m02[8];C_V3_0<-m02[8]/n22;C_V3_0[is.na(C_V3_0)] <-0
V4_0<-m02[9];C_V4_0<-m02[9]/n22;C_V4_0[is.na(C_V4_0)] <-0
V5_0<-m02[10];C_V5_0<-m02[10]/n22;C_V5_0[is.na(C_V5_0)] <-0
V6_0<-m02[11];C_V6_0<-m02[11]/n22;C_V6_0[is.na(C_V6_0)] <-0

Motif_counts<-c(T2,V1_0,V2_0,V3_0,V5_0,V6_0)
Motif_counts

Motif_concentration<-c(C_T2,C_V1_0,C_V2_0,C_V3_0,C_V4_0,C_V5_0,C_V6_0)
Motif_concentration


######### Only dense motifs ##########
Motif_counts2<-c(T2,V1_0,V2_0,V4_0,V5_0,V6_0)
Motif_counts2

Motif_concentration2<-c(C_T2,C_V1_0,C_V3_0,C_V4_0,C_V5_0,C_V6_0)
Motif_concentration2

#Motif_concentration
# 0.370866913 0.038725055 0.328102559 0.003358649 0.140127193 0.113031670



Motif_counts2/edge
# 2.549534  2.391957 23.265102  0.207456  8.655333  6.981705











#----------------------------------------------- 

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


ND0<-knn(G0)
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


meanD1<-mean(degree1);meanD1 #3.630326
sdD1<-sd(degree1);sdD1 # 5.824924


hist(degree1,xlab="Degree",col="blue",main="",first.panel=grid())
box()

summary(ND)

plot(degree1, ND, log="xy",
   col="goldenrod", xlab="Log node Degree",
     ylab="Log Average Neighbor Degree",first.panel=grid())



  

degree_D18<-degree.distribution(G0)


plot(degree_D18, type="o", lwd=2,xaxt="no",ylab = "p(k)",xlab = "Number of edge (Degree), k",
     main = "",xlim = c(),col="blue",first.panel=grid())


axis(1, at=c(1,  5,  10,  15,  20, 25, 30,35,40,45), 
     labels= c("0", "5", "10", "15", "20", "25","30","35","40","45"),
     cex=1.3)

######################################################################################

plot(degree_D18, type="o", lwd=2,xaxt="no",ylab = "p(k)",xlab = "Number of edge (Degree), k",
     main = "",xlim = c(0,45),col="blue",first.panel=grid())

lines(degree_D65,col="red",type="o", lwd=2)

axis(1, at=c(1,  5,  10,  15,  20, 25, 30,35,40,45), 
     labels= c("0", "5", "10", "15", "20", "25","30","35","40","45"),
     cex=1.3)




########## APL and Diameter ###########################

AVPL<-average.path.length(G0);AVPL # 20.8021
Diameter<-diameter(G0) ;Diameter # 65






#################### Assortativity and homophily ###########################

assortativity_degree(G11) # 0.8938759




assortativity(G11, V(G11)$vac18, directed=F) #  0.4579363


assortativity_nominal(G11, as.integer(as.factor(V(G11)$state), directed=F)) #  0.7693208


assortativity(G11, V(G11)$Covid_cases, directed=F) # 0.4831742



assortativity(G11, V(G11)$Edu_ass_coll, directed=F) # 0.5696295

assortativity(G11, V(G11)$Income, directed=F) # 0.5751481







