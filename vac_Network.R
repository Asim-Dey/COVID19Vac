

D0<-read.csv('C:/Users/adey/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/County_VaccinationsAA.csv')

D0<-read.csv('C:/Users/asimi/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/County_VaccinationsAA.csv')

names(D0)
length(unique(D0$STATE))


###################### Distance matrix #############################

#================= Vec Distance matrix ===================

#Series_Complete_18PlusPop_Pct
vec18<-D0$Series_Complete_18PlusPop_Pct

vec65<-D0$Series_Complete_65PlusPop_Pct


DD1<-as.matrix(vec65)

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


#summary(as.numeric(M11))
Mx<-max(as.numeric(M11))

M_Vec<- (M11/Mx) 
dim(M_Vec)

summary(as.numeric(M_Vec))

# Min.    1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.00000 0.07508 0.16817 0.21136  0.31832   1.00000 



##===============  Geo Distance matrix =======================

library(geosphere)

M_Goe<-distm(D0[,c(7,8)],fun = distHaversine)/ 1609.344 # Distance in miles
dim(M_Goe)


#########################################################################


k1=0.01 # Counties which vec rate similar dist less than k2 (less than median)
DM=100  # Counties which are within "DM" miles are connected by an edge


n=dim(M_Vec)[1]
A <- matrix(0, nrow = n, ncol =n);dim(A)

A2 <- matrix(0, nrow = n, ncol =n);dim(A2)



for (i in 1:n){ #i=1
  for (j in 1:n){ #j=2
    
    
    #if  (M_Vec[i,j] < k1 & M_Goe[i,j] < DM)  A2[i,j]=1
    
    if  (M_Vec[i,j] < k1 & M_Goe[i,j] < DM)  A2[i,j]=M_Vec[i,j]
    
    #if  (d11[,1][i] > k2 & d11[,2][j]> k2 & M_Goe[i,j] < DM)  A2[i,j]=1
    # if  (M[i,j] < DM)  A[i,j]=1 
    
  }
}


###################### Network ###################

#install.packages("igraph")
library(igraph) 

G1<- graph.adjacency(A2, mode="undirected",weighted=TRUE,diag=FALSE)


V<-length(V(G1));V
E<-length(E(G1));E


edlist_65<-cbind(get.edgelist(G1) , (1/round(E(G1)$weight, 3)))
#write.csv(edlist_65,'C:/Users/adey/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/edlist_65.csv')

#--------------------------------------------------


G1<- graph.adjacency(A2, mode="undirected",weighted=TRUE,diag=FALSE)%>% 
 set_vertex_attr("LON", value = D0$LON)%>%
           set_vertex_attr("LAT", value = D0$LAT)%>% 
                    set_vertex_attr("name", value = D0$FIPS)


as_tbl_graph(G1)


Atr<-get.vertex.attribute(G1)

LON1<-Atr$LON
LAT1<-Atr$LAT
FIPS<-Atr$name

DD1<-data.frame(LON1,LAT1,FIPS)
dim(DD1)













plot(G1, layout=layout_with_fr, vertex.size=1, vertex.label=NA,vertex.color="red")



library(maps)
library(geosphere)


map("world", regions=c("usa"), fill=T, col="grey90",  ylim=c(21.0,50.0), 
    xlim=c(-130.0,-65.0))
points(D0$LON,D0$LAT, pch=19, cex=0.1, col="blue")







g <- make_ring(10) %>%
  set_vertex_attr("color", value = "red") %>%
  set_vertex_attr("label", value = letters[1:10])
vertex_attr(g, "label")
vertex_attr(g)
plot(g)








