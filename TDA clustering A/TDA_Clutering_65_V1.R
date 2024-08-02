
library(data.table)
library(dplyr)
library(partykit)
library(ggplot2)
library(plyr)
library(cluster)
  library(NbClust)
library(factoextra)
library(purrr)

library(cluster)
library(TDA)
library(igraph)

##################################################################

setwd("C:/Users/asimi/OneDrive/COVID-19 Vac Network/Data and Code/R Codes/TDA Clustering/")

data0 <-read.table("node2vec_65.txt", quote="\"", comment.char="")
dim(data0) #  3192   50

# So, we now have a 50-dimensional vector associated with each node 
# in our network

data <- as.matrix(dist(data0))


#summary(as.numeric(data))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.000   4.326   4.554   4.533   4.817   8.464  


###########################################################################

cap=2# dimension cap
nodes=nrow(data) # number of points
nNghb=30 # number of neighbors around each data point
delta=0.280 # step size for filtration
filt_len=30 # filtration length.


# When distances are transformed, delta=0.01 and filt_len=L 
#implies that we account for L% of variation in distances
#data <- as.matrix(data)
# Running Perseus
#betti_0=betti_1=c()
wasser_dist_0 = wasser_dist_1 =wasser_dist_2 = c()
Persistance =c()
index=c()


for (i in 1:nodes){ #i=1
  
  ind=which(rank(data[i,],ties.method = 'first')<=nNghb) # indices of points around point i
  
  # writing data into file M.txt
  cat(length(ind),file='M.txt',append=F,sep = '\n')
  cat(paste(0,delta,filt_len,cap,sep = ' '),file = 'M.txt',append = T,sep = '\n') # threshold: 0, 0.1, 0.2,...,1
  cat(data[ind,ind],file='M.txt',append = T)
  
  system("perseusWin.exe distmat M.txt OutputFile") # for Windows
  print(i) 
  index=rbind(index,ind)
  betti_data=as.matrix(read.table('OutputFile_betti.txt'))
  # dim=0
  persist_data = as.matrix(read.table('OutputFile_0.txt'))
  persist_data[persist_data[,2] == -1, 2] = filt_len + 1
  persist_data = persist_data/(filt_len + 1)
  P = cbind(rep(0, nrow(persist_data)), persist_data)
  
  # dim=1
  if (file.info('OutputFile_1.txt')$size>0)
  { 
    persist_data = as.matrix(read.table('OutputFile_1.txt', blank.lines.skip = T))
    persist_data[persist_data[,2] == -1, 2] = filt_len + 1
    persist_data = persist_data/(filt_len + 1)
    P = rbind(P, cbind(rep(1, nrow(persist_data)), persist_data))
    
  }
  
  if (file.info('OutputFile_2.txt')$size>0)
  { 
    persist_data = as.matrix(read.table('OutputFile_2.txt', blank.lines.skip = T))
    persist_data[persist_data[,2] == -1, 2] = filt_len + 1
    persist_data = persist_data/(filt_len + 1)
    P = rbind(P, cbind(rep(2, nrow(persist_data)), persist_data))
    
  }
  
  Persistance[[i]] <- P
  
}








###########################################################################################
wasser_dist_2 = matrix(0,nodes,nodes)

for(i in 1:nodes){
  for(j in 1:nodes){
    wasser_dist_2[i,j] = wasserstein(Persistance[[i]], Persistance[[j]], dimension = c(0,1))
  }
  print(i)
}


# summary(as.vector(wasser_dist_2))
#   Min.   1st Qu.  Median    Mean     3rd Qu.    Max. 
# 0.0000  0.2581     0.5161  0.7708    1.0645  2.8226  




#write.csv(wasser_dist_2,'C:/Users/asimi/OneDrive/COVID-19 Vac Network/Data and Code/Data/COVID-19 Vaccinations in the US/wasser_dist_2_65_d50_V1.csv')

wasser_dist_2<-read.csv('C:/Users/asimi/OneDrive/COVID-19 Vac Network/Data and Code/Data/COVID-19 Vaccinations in the US/wasser_dist_2_65_d50_V1.csv',header=FALSE)
#wasser_dist_2<-read.csv('C:/Users/adey/OneDrive/COVID-19 Vac Network/Data and Code/Data/COVID-19 Vaccinations in the US/wasser_dist_2_18_d50_V1.csv',header=FALSE)
head(wasser_dist_2)
dim(wasser_dist_2)



############################################################################################

quantile(wasser_dist_2,c(0.1,0.2,0.25,0.3,0.35,0.4))
#cutoff_2 = quantile(wasser_dist_2,0.25)

## Elbow Plot#####

summary_CPD <- data.frame(cutoff = double(),
                          WithinSumSquares= double(),
                          BetweenSumSquares= double(),
                          NoOfClusters = integer())


############################################################################################
########################## Begin Loop ##################################################

cutoff = seq(0.1,2.5,length.out=20)
test =c()

for (c in cutoff){ #c=1.5
  cutoff_2 = c
  A=matrix(0,ncol = nodes,nrow = nodes) # Forming adjacency matrix
  
  for (i in 1:nodes)  {
    index_2 = which(wasser_dist_2[i,]<=cutoff_2)
    A[i,index_2]=1
  }
  
  
  
  # Clustering
  g=graph_from_adjacency_matrix(A,"directed") # form a graph from adj. matrix A
  clstrs=clusters(g)  # clusters are strongly connected components of adj. matrix A
  
  center =c()
  total_wss = 0
  
  for(i in 1:clstrs$no) {
    if(clstrs$csize[i] > 1) {
      indices <- which(clstrs$membership == i,arr.ind = T) ## all the indices in the cluster
      wasser_0 <- wasser_dist_2[indices,indices] #subset dataset for those points
      index <- which.min(rowSums(wasser_0)) # find centroid of the cluster
      point_0 <- indices[index] #put centroid in point_0
      wss <- sum((wasser_0[index,])^2)
      wss <- min(rowSums(wasser_dist_2[indices,indices]))
      
    }
    else
    {
      point_0 <- which(clstrs$membership == i,arr.ind = T)
      wss <- 0
      
    }
    
    center <- c(center,point_0)
    total_wss = total_wss + wss
  }
  
  wasser_cen <- wasser_dist_2[center,center]
  total_bss = sum(wasser_cen)/2
  
  test <- c(cutoff_2,total_wss,total_bss,clstrs$no)
  #test <- c(cutoff_2,avg_sil,clstrs$no)
  #print(test)
  summary_CPD <- rbind(test,summary_CPD)
}




############################# End Loop ####################################################
###########################################################################################

summaryCPD=summary_CPD
colnames(summaryCPD) <- c("cutOff","WithinSumSquares","BetweenSumSquares","NoOfClusters")
summaryCPD


summaryCPD$WBRatio = summaryCPD$WithinSumSquares/summaryCPD$BetweenSumSquares
summaryCPD

  # cutOff    WithinSumSquares   BetweenSumSquares NoOfClusters WBRatio

#12 1.1105263        1061.2258          2.241935            2 473.35252
#13 0.9842105        1061.2258          2.241935            2 473.35252
#14 0.8578947        1061.2258          2.241935            2 473.35252
#15 0.7315789        1061.2258          2.241935            2 473.35252

#16 0.6052632        1061.2258          2.241935            2 473.35252
#17 0.4789474         940.3871          5.016129            3 187.47267
#18 0.3526316         939.3065         10.435484            4  90.01082
#19 0.2263158         848.9839         57.822581           10  14.68257
#20 0.1000000         702.0484       1716.709677           60   0.40895


summaryCPD$WBRatio = summaryCPD$WithinSumSquares/summaryCPD$BetweenSumSquares


 
plot(summaryCPD$cutOff,summaryCPD$WBRatio)
plot(summaryCPD$NoOfClusters,summaryCPD$WBRatio,xlim=c(0,40))

plot(summaryCPD$NoOfClusters,summaryCPD$WBRatio, main="CPD",cex.main=0.8,
     type='o',pch = 1,cex = .7,col='green',
     xlab="Number of cluster", ylab="WB-ratio",xlim=c(),first.panel=grid())

#plot(summaryCPD$cutOff,summaryCPD$BetweenSumSquares)
#ggplot(data = summaryCPD,aes(x=cutOff,y=WithinSumSquares))+geom_line()
plot(summaryCPD$NoOfClusters,summaryCPD$WithinSumSquares)
#plot(summaryCPD$NoOfClusters,summaryCPD$cutOff)



library(readr)
#Performance_N18 <- read.table("TDA Clustering/Performance_N18.txt")
Performance_N65 <-read.table("~/Library/CloudStorage/OneDrive-Personal/Mixed Effect Model/Data and Code/R Codes/TDA Clustering/Performance_N65.txt",header=TRUE)


Performance_N65$WBRatio<-Performance_N65$WithinSumSquares/Performance_N65$BetweenSumSquares
Performance_N65

Performance_N65<-Performance_N65[1:8,]
Performance_N65



plot(Performance_N65$NoOfClusters,Performance_N65$WBRatio, main="",
     type='o',pch = 1,cex = .7,col='blue',
     xlab="Number of cluster", ylab="WB-ratio",xlim=c(2,10),first.panel=grid())

abline(v=4, col="green", lty=2)


plot(Performance_N65$NoOfClusters,Performance_N65$WithinSumSquares, main="",
     type='o',pch = 1,cex = .7,col='blue',
     xlab="Number of cluster", ylab="WSS",xlim=c(2,10),first.panel=grid())

abline(v=4, col="green", lty=2)



data.frame(Performance_N65$cutOff,Performance_N65$WithinSumSquares,Performance_N65$BetweenSumSquares,
           Performance_N65$WBRatio, Performance_N65$NoOfClusters)



#summary_CPD[summaryCPD$NoOfClusters ==4,]
#cutOff         WithinSumSquares   BetweenSumSquares  NoOfClusters    WB-Ratio      
# 0.10795918      91.25806452         44.096774           13       2.069495e+00

########################################################################################

cutoff_A = 0.353
A=matrix(0,ncol = nodes,nrow = nodes) # Forming adjacency matrix

for (i in 1:nodes)  {
  index_2 = which(wasser_dist_2[i,]<=cutoff_A)
  A[i,index_2]=1
}



# Clustering
g=graph_from_adjacency_matrix(A,"directed") # form a graph from adj. matrix A
clstrs_CPD=clusters(g)  # clusters are strongly connected components of adj. matrix A

cluster_member_CPD<-clstrs_CPD$membership


#write.csv(cluster_member_CPD,"cluster_member_TDA_65.csv")











##########################################################################################
#################################### Non TDA clusterings ###############################
#########################################################################################

############################# K-means #########################################


#data <- as.matrix(dist(data0))


#summary(as.numeric(data))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0000  0.1371  3.9669  3.2006  4.5080  8.5873 



fviz_nbclust(scale(data), FUN = kmeans, method = "wss",k.max=20)


################ SS ####################################################

kmean <- data.frame(WithinSumSquares= double(),
                    BetweenSumSquares= double(),
                    NoOfClusters = integer())

cutoff = seq(2,10) # number of clusters




test =c()


################################# Loop Start ##################################

for (c in cutoff) { # c=18
  
  # Clustering
  
  model1 <- kmeans(scale(data), centers = c, nstart = 8)
  
  #Kmedoids= pam(wasser_dist_2, k= c, diss = TRUE,medoids = NULL, stand = FALSE, 
  #             cluster.only = FALSE,do.swap = FALSE, pamonce = FALSE, trace.lev = 0)
  
  center =c()
  total_wss = 0
  
  for(i in 1:c){ # i=3
    
    
    indices <- which(model1$cluster == i,arr.ind = T) ## all the indices in the cluster
    wasser_0 <- wasser_dist_2[indices,indices] #subset dataset for those points
    index <- which.min(rowSums(wasser_0)) # find centroid of the cluster
    point_0 <- indices[index] #put centroid in point_0
    wss <- sum((wasser_0[index,])^2)
    wss <- min(rowSums(wasser_dist_2[indices,indices]))
    #print(wss)
    center <- c(center,point_0)
    #print(center)
    total_wss = total_wss + wss
    
  }
  
  wasser_cen <- wasser_dist_2[center,center]
  total_bss = sum(wasser_cen)/2
  
  test <- c(total_wss,total_bss,c)
  kmean <- rbind(test,kmean)
}


################################# Loop End ##################################




colnames(kmean) <- c("WithinSumSquares","BetweenSumSquares","NoOfClusters")


kmean$WBRatio = kmean$WithinSumSquares/kmean$BetweenSumSquares


plot(kmean$NoOfClusters,kmean$WBRatio)

plot(kmean$NoOfClusters,kmean$WBRatio, main="k-mean",cex.main=0.8,
     type='o',pch = 1,cex = .7,col='green',
     xlab="Number of cluster", ylab="WB-ratio",first.panel=grid())


plot(kmean$NoOfClusters,kmean$WithinSumSquares)

 
#WithinSumSquares BetweenSumSquares NoOfClusters    WBRatio
#         745.0161        39.3064516           10   18.95404
#         754.6774        33.9354839            9   22.23859
#         757.2581        28.8064516            8   26.28779
#         740.3387        23.6935484            7   31.24643

#         743.3548        18.8064516            6   39.52659
#         768.1452        14.1129032            5   54.42857
#         858.5968         9.3225806            4   92.09862
#        1001.2581         5.4193548            3  184.75595
#        1537.5806         0.5322581            2 2888.78788

###############################################################################
################ K Medoids --- Euclidian distance metrics########################



ElbowPAM <- data.frame(WithinSumSquares= double(),
                       BetweenSumSquares= double(),
                       NoOfClusters = integer())

cutoff = seq(2,10) # number of clusters

test =c()


################################# Loop Start ##################################

for (c in cutoff) { # c=2
  
  # Clustering
  Kmedoids= pam(data, k= c, diss = TRUE,medoids = NULL, stand = FALSE, 
                cluster.only = FALSE,do.swap = FALSE, pamonce = FALSE, trace.lev = 0)
  
  center =c()
  total_wss = 0
  
  for(i in 1:c){ # i=1
    
    indices <- which(Kmedoids$clustering == i,arr.ind = T) ## all the indices in the cluster
    wasser_0 <- wasser_dist_2[indices,indices] #subset dataset for those points
    index <- which.min(rowSums(wasser_0)) # find centroid of the cluster
    point_0 <- indices[index] #put centroid in point_0
    wss <- sum((wasser_0[index,])^2)
    wss <- min(rowSums(wasser_dist_2[indices,indices]))
    center <- c(center,point_0)
    total_wss = total_wss + wss
  }
  
  wasser_cen <- wasser_dist_2[center,center]
  total_bss = sum(wasser_cen)/2
  
  test <- c(total_wss,total_bss,c)
  ElbowPAM <- rbind(test,ElbowPAM)
}


################################# Loop End ##################################

colnames(ElbowPAM) <- c("WithinSumSquares","BetweenSumSquares","NoOfClusters")
ElbowPAM$WBRatio = ElbowPAM$WithinSumSquares/ElbowPAM$BetweenSumSquares

plot(ElbowPAM$NoOfClusters,ElbowPAM$WBRatio)
plot(ElbowPAM$NoOfClusters,ElbowPAM$WBRatio, main="K-medoids (Euclidean)",cex.main=0.8,
     type='o',pch = 1,cex = .7,col='green',
     xlab="Number of cluster", ylab="WB-ratio",first.panel=grid())



#   WithinSumSquares BetweenSumSquares NoOfClusters   WBRatio
#         876.5323         28.758065           10  30.47953
#         889.5806         25.032258            9  35.53737
#         904.6129         20.822581            8  43.44384
#         924.0484         17.516129            7  52.75414
#         933.5645         14.032258            6  66.52989
#         944.0645         10.580645            5  89.22561

#         970.5484          7.693548            4 126.15094
#         997.5968          4.774194            3 208.95608
#        1061.2258          2.241935            2 473.35252









###############################################################################################
############################# Spectral Clustering #########################################

nm=4

data<-wasser_dist_2
#data <- as.matrix(dist(data0))


S <- as.matrix(wasser_dist_2)
dim(S)

## Create Degree matrix
D <- matrix(0, nrow=nrow(data), ncol = nrow(data)) # empty nxn matrix

for (i in 1:nrow(data)) {
  
  # Find top 10 nearest neighbors using Euclidean distance
  index <- order(S[i,])[2:11]
  
  # Assign value to neighbors
  D[i,][index] <- 1 
  
  print(i)
}

# find mutual neighbors
D = D + t(D) 
D[ D == 2 ] = 1

# find degrees of vertices
degrees = colSums(D) 
n = nrow(D)

## Compute Laplacian matrix
# Since k > 2 clusters (3), we normalize the Laplacian matrix:
laplacian = ( diag(n) - diag(degrees^(-1/2)) %*% D %*% diag(degrees^(-1/2)) )
dim(laplacian)



## Compute eigenvectors
eigenvectors = eigen(laplacian, symmetric = TRUE)
n = nrow(laplacian)
eigenvectors = eigenvectors$vectors[,(n - 2):(n - 1)]


set.seed(1748)
## Run Kmeans on eigenvectors
sc = kmeans(eigenvectors, nm)
sc

## Pull clustering results
sc_results = cbind(data, cluster = as.factor(sc$cluster))
head(sc_results)




################ SS ####################################################

SC2 <- data.frame(WithinSumSquares= double(),
                  BetweenSumSquares= double(),
                  NoOfClusters = integer())


test =c()


################################# Loop Start ##################################


 
c<-4
center =c()
total_wss = 0

for(i in 1:c){ # i=1
  
  
  indices <- which(sc$cluster == i,arr.ind = T) ## all the indices in the cluster
  wasser_0 <- wasser_dist_2[indices,indices] #subset dataset for those points
  index <- which.min(rowSums(wasser_0)) # find centroid of the cluster
  point_0 <- indices[index] #put centroid in point_0
  wss <- sum((wasser_0[index,])^2)
  wss <- min(rowSums(wasser_dist_2[indices,indices]))
  #print(wss)
  center <- c(center,point_0)
  #print(center)
  total_wss = total_wss + wss
  
}



wasser_cen <- wasser_dist_2[center,center]
total_bss = sum(wasser_cen)/2

test <- c(total_wss,total_bss,c)
SC2 <- rbind(test,SC2)
 
 
   
################################# Loop End ##################################




colnames(SC2) <- c("WithinSumSquares","BetweenSumSquares","NoOfClusters")
SC2$WBRatio = SC2$WithinSumSquares/SC2$BetweenSumSquares


SC2

# WithinSumSquares BetweenSumSquares NoOfClusters  WBRatio
#        1092.516          6.870968            4 159.0047

