
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
#setwd("~/Library/CloudStorage/OneDrive-Personal/Mixed Effect Model/Data and Code/R Codes/TDA Clustering/")

setwd("C:/Users/asimi/OneDrive/Mixed Effect Model/Data and Code/R Codes/TDA Clustering/")

data0 <-read.table("node2vec_18_DM30.txt", quote="\"", comment.char="")
dim(data0) #  3192   50

#write.csv(data0,"dataAB.csv")


# So, we now have a 50-dimensional vector associated with each node 
# in our network

data <- as.matrix(dist(data0))


#summary(as.numeric(data))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.0000  0.1371  0.6522  0.5771  0.8549  0.9540 


###########################################################################

cap=2# dimension cap
nodes=nrow(data) # number of points
nNghb=30 # number of neighbors around each data point
delta=0.031 # step size for filtration
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
# 0.0000    0.3226  1.1290    3.7832    7.2258  7.9355 


############################################################################################

quantile(wasser_dist_2,c(0.1,0.2,0.25,0.3,0.35,0.4 ))
#cutoff_2 = quantile(wasser_dist_2,0.25)

## Elbow Plot#####

summary_CPD <- data.frame(cutoff = double(),
                          WithinSumSquares= double(),
                          BetweenSumSquares= double(),
                          NoOfClusters = integer())


############################################################################################
########################## Begin Loop ##################################################

cutoff = seq(0.1,1.5,length.out=20)
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

cutOff WithinSumSquares BetweenSumSquares NoOfClusters
1  1.5000000         821.1129          7.193548            2
2  1.4263158         821.1129          7.193548            2
3  1.3526316         821.1129          7.193548            2
4  1.2789474         821.1129          7.193548            2
5  1.2052632         821.1129          7.193548            2
6  1.1315789         821.1129          7.193548            2
7  1.0578947         821.1129          7.193548            2
8  0.9842105         821.1129          7.193548            2
9  0.9105263         821.1129          7.193548            2
10 0.8368421         821.1129          7.193548            2
11 0.7631579         821.1129          7.193548            2
12 0.6894737         821.1129          7.193548            2
13 0.6157895         821.1129          7.193548            2
14 0.5421053         821.1129          7.193548            2
15 0.4684211         558.2258         30.161290            4
16 0.3947368         558.2258         30.161290            4
17 0.3210526         558.2258         30.161290            4
18 0.2473684         558.2258         30.161290            4
19 0.1736842         557.6774         45.838710            5
20 0.1000000         541.8710        455.967742           23


summaryCPD$WBRatio = summaryCPD$WithinSumSquares/summaryCPD$BetweenSumSquares



plot(summaryCPD$cutOff,summaryCPD$WBRatio)
plot(summaryCPD$NoOfClusters,summaryCPD$WBRatio,xlim=c())

plot(summaryCPD$NoOfClusters,summaryCPD$WBRatio, main="CPD",cex.main=0.8,
     type='o',pch = 1,cex = .7,col='green',
     xlab="Number of cluster", ylab="WB-ratio",xlim=c(),first.panel=grid())

#plot(summaryCPD$cutOff,summaryCPD$BetweenSumSquares)
#ggplot(data = summaryCPD,aes(x=cutOff,y=WithinSumSquares))+geom_line()
plot(summaryCPD$NoOfClusters,summaryCPD$WithinSumSquares)
#plot(summaryCPD$NoOfClusters,summaryCPD$cutOff)








#summary_CPD[summaryCPD$NoOfClusters ==4,]
#cutOff         WithinSumSquares   BetweenSumSquares  NoOfClusters    WB-Ratio      
# 0.10795918      91.25806452         44.096774           13       2.069495e+00

########################################################################################

cutoff_A = 0.47
A=matrix(0,ncol = nodes,nrow = nodes) # Forming adjacency matrix

for (i in 1:nodes)  {
  index_2 = which(wasser_dist_2[i,]<=cutoff_A)
  A[i,index_2]=1
}



# Clustering
g=graph_from_adjacency_matrix(A,"directed") # form a graph from adj. matrix A
clstrs_CPD=clusters(g)  # clusters are strongly connected components of adj. matrix A



cluster_member_CPD<-clstrs_CPD$membership


write.csv(cluster_member_CPD,"cluster_member_TDA_18_30DM.csv")











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










###############################################################################
################### Hierarchical Clustering ###################################

fviz_nbclust(data, FUN = hcut, method = "wss",k.max=20)

model2 <- hclust(scale(data))
plot(model2)

# function to compute average silhouette for k clusters
avg_sil <- function(k) {
  heir.res <- hclust(data)
  clust_heir <- cutree(heir.res,k=k)
  ss <- cluster::silhouette(clust_heir, dist(data))
  mean(ss[, 3])
}


k.values <- 2:20

# extract avg silhouette for 2-15 clusters
avg_sil_values1 <- map_dbl(k.values, avg_sil)


plot(k.values, avg_sil_values1,
     type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters K",
     ylab = "Average Silhouettes")



k_hc=10 # 10

clust_heir <- cutree(model2,k=k_hc)  # based on Eulidian Distance
rect.hclust(model2, k = k_hc, border = 2:k_hc)

cluster_size_heir <- table(clust_heir)
center_heir =c()
total_wss_heir = 0
test = c()

for(i in 1:length(cluster_size_heir))
{
  if(cluster_size_heir[i] > 1)
  {
    indices <- which(clust_heir == i,arr.ind = T) ## all the indices in the cluster
    wasser_0 <- wasser_dist_2[indices,indices] #subset dataset for those points
    index <- which.min(rowSums(wasser_0)) # find centroid of the cluster
    point_0 <- indices[index] #put centroid in point_0
    wss <- sum((wasser_0[index,])^2)
    wss <- min(rowSums(wasser_dist_2[indices,indices]))
    
  }
  else
  {
    point_0 <- which(clust_heir == i,arr.ind = T)
    wss <- 0
    
  }
  center_heir <- c(center_heir,point_0)
  total_wss_heir = total_wss_heir + wss
  print(i)
}

wasser_cen_heir <- wasser_dist_2[center_heir,center_heir]
total_bss_heir = sum(wasser_cen_heir)/2

test <- c(total_wss_heir,total_bss_heir)
#test <- c(cutoff_2,avg_sil,clstrs$no)
print(test)
#summary <- rbind(test,summary)
wb_ratio_heir <- test[1]/test[2]
#wb_ratio_heir



SS_kheir<-c(test[1],test[2],wb_ratio_heir)
names(SS_kheir)<-c("WSS", "BSS", "WB-Ratio")
SS_kheir

# k=10

#WSS        BSS   WB-Ratio 
#119.177419  14.096774   8.454233 

