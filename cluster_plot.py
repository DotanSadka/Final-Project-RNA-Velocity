import pandas as pd
import numpy as np
import matplotlib as plt
import plotly.express as px
import plotly.graph_objs as go


import scanpy as sc
from sklearn.neighbors import NearestNeighbors
import louvain
import plotly.express as px
import matplotlib.pyplot as plt

cd8a_data = {}

for t in [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0]:
    cd8a_data[t] = pd.read_csv(f"isCloseT/df/T-cell_CD8A_{t}.csv", index_col=0)
    cd8a_percent = 100 * (cd8a_data[t]["isClose"] == 1).sum() / len(cd8a_data[t])


    #s_all_PCs_subset = pd.read_csv('isCloseT/df/T-cell_CD8A_{t}.csv', index_col=0)
    s_all_PCs_subset = cd8a_data[t][['pc1', 'pc2', 'pc3']]

    # Create an AnnData object from the PCA results DataFrame
    adata = sc.AnnData(s_all_PCs_subset)

    # Find nearest neighbors
    sc.pp.neighbors(adata, n_neighbors=180) #change number of n_neighbors according to the number of clusters you would like to obtain

    # Perform Louvain clustering
    sc.tl.louvain(adata)

    # Access the cluster labels
    cluster_labels = adata.obs['louvain']

    # Print the first 10 cluster labels
    #print(cluster_labels[:10])
    cluster_labels = cluster_labels.reset_index(drop=True)
    s_all_PCs_subset = s_all_PCs_subset.reset_index()#drop=True)
    s_all_PCs_subset['Cluster'] = cluster_labels.astype(str)

    # type(s_all_PCs_subset['Cluster'])
    second_column = list(s_all_PCs_subset['Cluster'])
    second_column = [int(x) for x in second_column]

    # Create the 3D scatter plot PC1 PC2
    names = cd8a_data[t]['Unnamed: 0']
    # points = [(x,y) for x,y in zip(s_all_PCs_subset['pc1'],s_all_PCs_subset['pc2'])]
    plt.scatter(s_all_PCs_subset['pc1'],s_all_PCs_subset['pc2'], s=200, edgecolors='black', c=second_column)
    plt.title(f'T-cell CD8A - PC1 PC2_{t}')
    plt.xlabel('PC1')
    plt.ylabel('PC2')

    plt.savefig(f'after cluster/PC1 PC2/PC1_PC2_{t}.png')
    plt.clf()
    # Show the plot
    #plt.show()

    #Create the 3D scatter plot PC1 PC3
    plt.scatter(s_all_PCs_subset['pc1'],s_all_PCs_subset['pc3'], s=200, edgecolors='black', c=second_column)
    plt.title(f'T-cell CD8A - PC1 PC3_{t}')
    plt.xlabel('PC1')
    plt.ylabel('PC3')

    plt.savefig(f'after cluster/PC1 PC3/PC1_PC3_{t}.png')
    plt.clf()


    ##Create the 3D scatter plot PC2 PC3

    plt.scatter(s_all_PCs_subset['pc2'],s_all_PCs_subset['pc3'], s=200, edgecolors='black', c=second_column)
    plt.title(f'T-cell CD8A - PC2 PC3_{t}')
    plt.xlabel('PC2')
    plt.ylabel('PC3')

    plt.savefig(f'after cluster/PC2 PC3/PC2_PC3_{t}.png')
    plt.clf()
