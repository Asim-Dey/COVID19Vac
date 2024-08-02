
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

setwd("C:/Users/adey/OneDrive/COVID-19 Vac Network/Data and Code/R Codes/TDA Clustering/")
setwd("~/Library/CloudStorage/OneDrive-Personal/Mixed Effect Model/Data and Code/R Codes/TDA Clustering/")


data0 <-read.table("node2vec_18.txt", quote="\"", comment.char="")
dim(data0) #  3192   50

#write.csv(data0,"dataAB.csv")


# So, we now have a 50-dimensional vector associated with each node 
# in our network

data <- as.matrix(dist(data0))  # matrix of pairwise distance 


#summary(as.numeric(data))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0000  0.1371  3.9669  3.2006  4.5080  8.5873 


###########################################################################

cap=2# dimension cap
nodes=nrow(data) # number of points
nNghb=30 # number of neighbors around each data point
delta=0.286 # step size for filtration
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




dim(wasser_dist_2) # 3192 3192

# summary(as.vector(wasser_dist_2))
# Min.   1st Qu.  Median    Mean     3rd Qu.    Max. 
# 0.000   0.000   1.355    1.197     2.194   2.645 



#write.csv(wasser_dist_2,'C:/Users/adey/OneDrive/COVID-19 Vac Network/Data and Code/Data/COVID-19 Vaccinations in the US/wasser_dist_2_18_d50_V1.csv')


wasser_dist_2<-read.csv('~/Library/CloudStorage/OneDrive-Personal/COVID-19 Vac Network/Data and Code/Data/COVID-19 Vaccinations in the US/wasser_dist_2_18_d50_V1.csv',header=FALSE)
wasser_dist_2<-read.csv('C:/Users/adey/OneDrive/COVID-19 Vac Network/Data and Code/Data/COVID-19 Vaccinations in the US/wasser_dist_2_18_d50_V1.csv',header=FALSE)
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

cutoff = seq(0.1,2,length.out=20)
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

#    cutOff WithinSumSquares BetweenSumSquares NoOfClusters

#10    1.1        3120.7258          0.000000            1
#11    1.0        3120.7258          0.000000            1
#12    0.9        3120.7258          0.000000            1
#13    0.8        3120.7258          0.000000            1
#14    0.7        3120.7258          0.000000            1
#15    0.6        3120.7258          0.000000            1
#16    0.5         689.8548          5.000000            3
#17    0.4         678.1452          8.241935            4
#18    0.3         676.5161         18.806452            6
#19    0.2         495.9677        118.967742           15
#20    0.1         318.7903       3035.225806           83

 
# WB-rartio for k=4
# 678.1452/8.241935=82.27985



 
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
Performance_N18 <- read.table("~/Library/CloudStorage/OneDrive-Personal/Mixed Effect Model/Data and Code/R Codes/TDA Clustering/Performance_N18.txt",header=TRUE)
Performance_N18

Performance_N18$WBRatio<-Performance_N18$WithinSumSquares/Performance_N18$BetweenSumSquares
Performance_N18



Performance_N18<-Performance_N18[1:10,]

plot(Performance_N18$NoOfClusters,Performance_N18$WithinSumSquares, main="",
     type='o',pch = 1,cex = .7,col='blue',
     xlab="Number of cluster", ylab="WSS",xlim=c(1,15),first.panel=grid())

abline(v=4, col="green", lty=2)


plot(Performance_N18$NoOfClusters,Performance_N18$WBRatio, main="",cex.main=0.8,
     type='o',pch = 1,cex = .7,col='blue',
     xlab="Number of cluster", ylab="WB Ratio",xlim=c(1,15),first.panel=grid())

abline(v=4, col="green", lty=2) 


plot(Performance_N18$WithinSumSquares,type="o")
#--------------------------------------------------------------

plot(Performance_N18$NoOfClusters,Performance_N18$BetweenSumSquares, main="CPD",cex.main=0.8,
     type='o',pch = 1,cex = .7,col='green',
     xlab="Number of cluster", ylab="WB-ratio",xlim=c(0,15),first.panel=grid())



data.frame(Performance_N18$cutOff,Performance_N18$WithinSumSquares,Performance_N18$BetweenSumSquares,
           Performance_N18$WBRatio, Performance_N18$NoOfClusters)


#summary_CPD[summaryCPD$NoOfClusters ==4,]
#cutOff         WithinSumSquares   BetweenSumSquares  NoOfClusters    WB-Ratio      
# 0.10795918      91.25806452         44.096774           13       2.069495e+00

########################################################################################

cutoff_A = 0.4
A=matrix(0,ncol = nodes,nrow = nodes) # Forming adjacency matrix

for (i in 1:nodes)  {
  index_2 = which(wasser_dist_2[i,]<=cutoff_A)
  A[i,index_2]=1
}



# Clustering
g=graph_from_adjacency_matrix(A,"directed") # form a graph from adj. matrix A
clstrs_CPD=clusters(g)  # clusters are strongly connected components of adj. matrix A



cluster_member_CPD<-clstrs$clstrs_CPD


#write.csv(cluster_member,"cluster_member_TDA_18.csv")




 






##########################################################################################
#################################### Non TDA clusterings ###############################
#########################################################################################

############################# K-means #########################################


#data <- as.matrix(dist(data0))


#summary(as.numeric(data))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0000  0.1371  3.9669  3.2006  4.5080  8.5873 

data<-wasser_dist_2

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
#         574.5968         58.790323           10   9.773663
#         510.5806         45.951613            9  11.111267
#         572.3710         36.693548            8  15.598681
#         587.4839         27.419355            7  21.425882
#         572.3710         19.661290            6  29.111567
#         587.4839         13.000000            5  45.191067
#         743.0323          8.612903            4  86.269663
#         713.8871          5.000000            3 142.777419
#         858.1935          2.193548            2 391.235294


###############################################################################
################ K Medoids --- Euclidian distance metrics ########################



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



#WithinSumSquares BetweenSumSquares NoOfClusters   WBRatio
#         767.8226         36.854839           10  20.83370
#         773.3065         30.774194            9  25.12841
#         778.1774         25.758065            8  30.21102
#         791.9516         19.032258            7  41.61102
#         810.9839         14.403226            6  56.30571
#         825.8871         10.629032            5  77.70106
#         833.9355          7.322581            4 113.88546
#         840.4677          4.564516            3 184.13074
#         858.1935          2.193548            2 391.23529










###############################################################################################
############################# Spectral Clustering #########################################

nm=4

data<-wasser_dist_2
#data <- as.matrix(dist(data0))


S <- wasser_dist_2
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

# WithinSumSquares BetweenSumSquares   NoOfClusters  WBRatio
#     858.1935          6.580645            4      130.4118





