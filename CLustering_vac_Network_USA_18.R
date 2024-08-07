

D0<-read.csv('C:/Users/asimi/OneDrive/COVID-19 Vac Network/Data and Code/Data/COVID-19 Vaccinations in the US/County_VaccinationsAA.csv')
head(D0)

state<-unique(D0$state)
m<-length(state);m

names(D0)
state_nm<-state


D1<-D0[,c(5,12)]
head(D1)

k=20
M11<- matrix(data=NA,ncol=k,nrow=m)

for (j in 1:m){ #j=8
 M11[j,]<-D1[,2][D1[,1]==state[j]][1:k]

}





## Replace NA with the last value of the row---in the csv file manually 


D10<-cbind(state_nm,M11)
#write.csv(D10,'C:/Users/adey/OneDrive/Mixed Effect Model/Data and Code/Data/COVID-19 Vaccinations in the US/county_vex_rate_statwise18.csv')

  

 
############################################################################################################################################################
############################################################################################################################################################


library(factoextra)
library(cluster)


D11<-read.csv('C:/Users/asimi/OneDrive/COVID-19 Vac Network/Data and Code/Data/COVID-19 Vaccinations in the US/county_vex_rate_statwise18.csv',header=TRUE)
head(D11)


D12<-data.frame(D11, row.names = 1)
head(D12)

#scale the dataset
D12 <- scale(D12)
head(D12)




# Correlation-based distance method
res.dist_Vec <- get_dist(D12, method = "pearson")


head(round(as.matrix(res.dist_Vec), 2))#[, 1:6]


# Visualize the dissimilarity matrix
fviz_dist(res.dist_Vec, lab_size = 5)


###########################################################################
###################### Clustering ######################################

# Hierarchical clustering


gap_stat_HC <- clusGap(D12, FUN = hcut,K.max = 10)
fviz_gap_stat(gap_stat_HC)




res.hc <- eclust(D12, FUNcluster="hclust", k=3) 
fviz_dend(res.hc, rect = TRUE) # dendrogam
res.hc <- hclust(res.dist_Vec, method = "ward.D2")



# Compute hierarchical clustering
res.hc <- D12 %>%
  scale() %>%                    # Scale the data
  dist(method = "euclidean") %>% # Compute dissimilarity matrix
  hclust(method = "ward.D2")     # Compute hierachical clustering

# Visualize using factoextra
# Cut in k groups and color by groups
fviz_dend(res.hc, k = 3, # Cut in four groups
          cex = 0.7, # label size
          k_colors = c("red", "blue","black"),
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE # Add rectangle around groups
)



fviz_cluster(res.hc)

hc.cut <- hcut(D12, k = 3, hc_method = "complete")
# Visualize dendrogram
fviz_dend(hc.cut, show_labels = TRUE, rect = TRUE)



# Visualize cluster
fviz_cluster(hc.cut, ellipse.type = "convex")


 

###########################################################################################
# Compute and visualize k-means clustering



fviz_nbclust(D12, kmeans,
             method = "gap_stat")

set.seed(123)
km.res <- kmeans(D12, 3, nstart = 25)

# Visualize
library("factoextra")
fviz_cluster(km.res, data = D12,
             ellipse.type = "convex",
             palette = "jco",
             ggtheme = theme_minimal())




 

# Compute and visualize PAM clustering

fviz_nbclust(D12, pam, method = "gap_stat")
fviz_nbclust(D12, pam, method = "silhouette")



pam.res <- pam(D12, 3)
# Visualize
fviz_cluster(pam.res)


