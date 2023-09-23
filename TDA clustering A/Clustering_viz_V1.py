

import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns



from node2vec import Node2Vec
from sklearn.cluster import KMeans
from sklearn import cluster
from sklearn import metrics




import plotly
plotly.tools.set_credentials_file(username='moniyuv', api_key='w2NBdksdtxziaW5GHpbG')
import plotly.plotly as py




############################# Implementing the Node2vec algorithm #################################

dataAB = pd.read_csv('dataAB.csv') # Load data



#dataAB = pd.read_csv('C:/Users/adey/OneDrive/Mixed Effect Model/Data and Code/R Codes/TDA Clustering/dataAB.csv') # Load data

dataAB.head()
dataAB.shape #  (3191, 50) # So, we now have a 50-dimensional vector associated with each node in our network.




from umap import UMAP

umapper = UMAP(metric="cosine")
umap_node2vec_embedding = umapper.fit_transform(dataAB)


fig, ax = plt.subplots(1, 1, figsize=(10, 10))

sns.scatterplot(
    x=umap_node2vec_embedding[:, 0],
    y=umap_node2vec_embedding[:, 1],
    ax=ax,
    s=3,
    alpha=.5,
    linewidth=0.1,
)
_ = ax.axis("off")

