
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt



from node2vec import Node2Vec
from sklearn.cluster import KMeans 
from sklearn import cluster
from sklearn import metrics




import plotly

plotly.tools.set_credentials_file(username='moniyuv', api_key='w2NBdksdtxziaW5GHpbG')


import plotly.plotly as py




############################# Implementing the Node2vec algorithm #################################

data = pd.read_csv('AM_65_Dist_V2_DM30.csv') # Load data

data.head()
data.shape #  (3192, 3192) # Adjacency Matrix


A=np.matrix(data)
A.shape  # (3192, 3192)


# Create a graph
G=nx.from_numpy_matrix(A)

G.number_of_nodes() # 3192

####################################################################################################################################

# Precompute probabilities and generate walks - **ON WINDOWS ONLY WORKS WITH workers=1**
#node2vec = Node2Vec(G, dimensions=50, walk_length=50, num_walks=200, weight_key="V3")

node2vec = Node2Vec(G, dimensions=50, walk_length=50, num_walks=200)

#node2vec = Node2Vec(graph, dimensions=64, walk_length=30, num_walks=200, workers=4)  # Use temp_folder for big graphs



# Embed nodes
model = node2vec.fit(window=10, min_count=0, batch_words=4)  # Any keywords acceptable by gensim.Word2Vec can be passed, `dimensions` and `workers` are automatically passed (from the Node2Vec constructor)


########################################################################################################################################


emb_df = (
    pd.DataFrame(
        [model.wv.get_vector(str(n)) for n in G.nodes()],
        index = G.nodes
    )
)



emb_df.shape # 


################### Save e values #########################################################################################################
 
import os
 #os.getcwd()
 #os.chdir("/COVID-19 Vaccinations in the US")
 
 
 os.chdir("/TDA Clustering/")
 os.getcwd()


 ##### matrix in a text file ##########################
 
 # matrix in a text file
  mat = np.matrix(emb_df)
  
  with open('node2vec_65_DM30.txt','a') as f:
   for line in mat:
    np.savetxt(f, line, fmt='%.2f')
        



   







