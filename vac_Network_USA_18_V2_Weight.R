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


###################### Distance matrix #############################

#================= Vec Distance matrix ===================

#Series_Complete_18PlusPop_Pct
vec18<-D0$Series_Complete_18PlusPop_Pct
summary(vec18)

#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#   13.90   53.30   61.80   62.47   70.80   95.00 


DD1<-as.matrix(vec18)
n<-dim(DD1)[1]
n


################################################################################
m<-dim(DD1)[1];m
M11<- matrix(data=NA,ncol=m,nrow=m)

#DD1[,1]
#DD1[,2] 



for(i in 1:m){     # i=1
  for(j in 1:m){   # j=2
    
    
    M11[i,j] = abs(DD1[i,]- DD1[j,])
    dim(M11)
    
    
    # print(i)
    # print(j)
    
  }
}



S1<-summary(as.numeric(M11))
S1

#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00    5.90   12.50   14.92   21.60   81.10 



#Mx<-max(as.numeric(M11))
#M_Vec<- (M11/Mx) 
#dim(M_Vec)
#summary(as.numeric(M_Vec))


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
A3 <- matrix(0, nrow = n, ncol =n);dim(A3)


############### Step 1 #####################################

for (i in 1:n){ #i=1
  for (j in 1:n){ #j=2
    
    
    if  (DD1[,1][i] > k2 & DD1[,1][j]> k2)  A2[i,j]=(1/M11[i,j])
    
    if  (DD1[,1][i] > k2 & DD1[,1][j]> k2 & M_Goe[i,j] < DM)  A1[i,j]= (1/M11[i,j])
    if  (DD1[,1][i] > k2 & DD1[,1][j]> k2 & M_Goe[i,j] < DM)  A3[i,j]= 1
    

  }
  
  print(i)
}
   

 
dim(A1)
dim(A2)
dim(A3)

#write.csv(A3,'C:/Users/adey/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/edlist_18_Dist_V2_aboveMean_WeightA3.csv')




 


###################### Network ###################



G0<- graph.adjacency(A1, mode="undirected",weighted=TRUE,diag=FALSE)

V0<-length(V(G0));V0 # 3192
E0<-length(E(G0));E0 #  5794

plot(G0, vertex.size=5, 
     vertex.label=NA,vertex.color="red",edge.arrow.size=0.9)


# The proportion of present edges from all possible edges in the network.



edlist_18_Dist<-cbind(get.edgelist(G0))


W1<-round(E(G0)$weight,5)
summary(W1)

W1[is.infinite(W1)]<-25 ## put a large value
summary(W1)


edlist_18_with_Weight<-cbind(get.edgelist(G0), W1)


# write.csv(edlist_18_with_Weight,'C:/Users/adey/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/edlist_18_Dist_V2_aboveMean_Weight.csv')


