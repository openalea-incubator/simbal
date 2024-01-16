import numpy as np
from numpy import corrcoef
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.widgets import Slider, Button
from sklearn import model_selection
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc, accuracy_score, confusion_matrix
from sklearn.preprocessing import LabelEncoder
from sklearn.decomposition import PCA, NMF
from sklearn.manifold import TSNE
import seaborn as sns ; sns.set()
import numpy as np
import plotly.express as px
import umap

#### Datas analysis ####

###Data cleaning
impulse = pd.read_csv('C:/Users/FARRA/AppData/Local/GitHubDesktop/app-3.2.3/src/caracol/Stats/IMPULSE_training.csv', sep=';', encoding='latin-1', decimal='.')
impulse_df = pd.DataFrame(data=impulse)
print(impulse_df)
impulse_df.info()

var = ['GDD_rape', 'GDD_faba', 'DAS_rape', 'DAS_faba','Plant_position', 'TTT', 'Scale', 'Position', 'QN_roots', 'QC_roots', 'n_greenleaves_mainstem','n_dead_leaves_mainstem', 'n_emerged_leaves_mainstem', 'n_Branches','DM_dead_leaves_plant', 'DM_dead_leaf','DM_greenleaves_mainstem','DM_dead_leaf_or_litter','f_DM_aerial_to_whole_plant', 'QN_internode', 'QN_stem_branch', 'QN_greenleaf', 'QN_greenleaves', 'QN_mainstem', 'QN_dead_leaf', 'QN_dead_leaves', 'QN_dead_leaves_and_litter', 'QN_whole_plant_without_DL_or_litter', 'QC_internode', 'QC_greenleaf', 'QC_greenleaves', 'QC_mainstem', 'QC_dead_leaf', 'QC_dead_leaves']
features = pd.DataFrame(impulse_df[var])
print(features)
features.info()

features_drop = features.drop(features[(features['Position'] == 100)].index)
features_drop = features_drop.drop(features_drop[(features_drop['Position'] == 112)].index)
features_drop = features_drop.drop(features_drop[(features_drop['Position'] == 113)].index)
features_drop = features_drop.drop(features_drop[(features_drop['Position'] == 114)].index)
features_drop = features_drop.drop(features_drop[(features_drop['Position'] == 115)].index)
features_drop = features_drop.drop(features_drop[(features_drop['Position'] == 116)].index)
features_drop = features_drop.drop(features_drop[(features_drop['Position'] == 117)].index)
features_drop = features_drop.drop(features_drop[(features_drop['Scale'] == 3)].index)
features_drop = features_drop.drop(features_drop[(features_drop['Scale'] == 4)].index)
features_drop = features_drop.drop(features_drop[(features_drop['Scale'] == 5)].index)
              
features_clean = features_drop.apply(pd.to_numeric)
features_clean.info()

###Correlation matrix
corr_matrix = features_clean.corr(method='pearson', min_periods=1, numeric_only=False)
print(corr_matrix)

f, ax = plt.subplots(figsize=(12, 10))
sns.heatmap(corr_matrix, ax=ax, cmap="viridis", linewidths=0.1)
plt.show()

sns.pairplot(corr_matrix,hue=['Position'], vars=['GDD_rape','Plant_position', 'TTT', 'Scale', 'Position','QN_greenleaf', 'QN_greenleaves','QN_dead_leaf', 'QN_dead_leaves','QC_greenleaf', 'QC_greenleaves','QC_dead_leaf', 'QC_dead_leaves'])
plt.show()

####Clustering
###Non negative Matrix Factorization NMF
#data matrix X
X = features_clean.drop(features_clean[(features_clean['Position'] == 0)].index)#drop plant scale
X.fillna(0, inplace=True)#NaN imputation to 0

# Initializing model
nmf_model = NMF(n_components=9)
nmf_model.fit(X.values)

print (nmf_model.reconstruction_err_)

#looking clusters
n_top_words = 10

def print_top_words(model, feature_names, n_top_words):
    for topic_idx, topic in enumerate(model.components_):
        print("\nCluster #%d:" % topic_idx)
        print(" ".join([feature_names[i]
                        for i in topic.argsort()[:-n_top_words - 1:-1]]))
    print()

print("Clusters found via NMF:")
tfidf_feature_names = X.columns.tolist()
print_top_words(nmf_model, tfidf_feature_names, n_top_words)

#looking in which leaf position the cluster are "activated" the most.
nmf_embedding = nmf_model.transform(X)
nmf_embedding = (nmf_embedding - nmf_embedding.mean(axis=0))/nmf_embedding.std(axis=0)
top_idx = np.argsort(nmf_embedding,axis=0)[-3:]

count = 0
for idxs in top_idx.T: 
    print("\nCluster {}:".format(count))
    for idx in idxs:
        print(X.iloc[idx]['Position'])
    count += 1

#t-SNE
tsne = TSNE(n_components = 3, random_state=42, n_iter=1000)
tsne_embedding = tsne.fit_transform(nmf_embedding)
tsne_embedding = pd.DataFrame(tsne_embedding,columns=['x','y','z'])
#tsne_embedding = tsne.fit_transform(X)# or nmf_embedding ?
tsne_embedding['hue'] = nmf_embedding.argmax(axis=1)
tsne.kl_divergence_

tsne_embedding['GDD_rape'] = X['GDD_rape']
tsne_embedding['Position'] = X['Position']


print(tsne_embedding)

y = X['Position']


fig = px.scatter_3d(x=tsne_embedding['x'], y=tsne_embedding['y'], z=tsne_embedding['z'], color=y)
fig.update_layout(
    title="t-SNE visualization",
    xaxis_title="First t-SNE",
    yaxis_title="Second t-SNE"
)
fig.show()



topics = ['Cluster 1',
          'Cluster 2',
          'Cluster 3',
          'Cluster 4',
          'Cluster 5',
          'Cluster 6',
          'Cluster 7',
          'Cluster 8',
          'Cluster 9']

#code used to create the first plot for getting the colors 
plt.style.use('ggplot')

data = tsne_embedding
legend_list = []

fig, axs = plt.subplots(1,1, figsize=(5, 5), facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace = .1, wspace=.001)
scatter = plt.scatter(data=data,x='x',y='y',s=6,c=data['hue'],cmap="Set1")
plt.axis('off')
plt.show()

colors = []
for i in range(len(topics)):
    idx = np.where(data['hue']==i)[0][0]
    color = scatter.get_facecolors()[idx]
    colors.append(color)
    legend_list.append(mpatches.Ellipse((0, 0), 1, 1, fc=color))

print(colors)

colors = np.array([[0.89411765, 0.10196078, 0.10980392, 1.        ],
                   [0.21568627, 0.49411765, 0.72156863, 1.        ],
                   [0.30196078, 0.68627451, 0.29019608, 1.        ],
                   [0.59607843, 0.30588235, 0.63921569, 1.        ],
                   [1.        , 0.49803922, 0.        , 1.        ],
                   [1.        , 1.        , 0.2       , 1.        ],
                   [0.65098039, 0.3372549 , 0.15686275, 1.        ],
                   [0.96862745, 0.50588235, 0.74901961, 1.        ],
                   [0.6       , 0.6       , 0.6, 1.               ]])

legend_list = []
for i in range(len(topics)):   
    color = colors[i]
    legend_list.append(mpatches.Ellipse((0, 0), 1, 1, fc=color))

matplotlib.rc('font',family='monospace')
plt.style.use('ggplot')

fig, axs = plt.subplots(2,2, figsize=(10, 15), facecolor='w', edgecolor='k')
fig.subplots_adjust(hspace = .1, wspace=0)

axs = axs.ravel()

count = 0
legend = []
for year, idx in zip([633.4, 653.9, 921.7, 943.6 ], range(4)):
    data = tsne_embedding[tsne_embedding['GDD_rape']<=year]
    scatter = axs[idx].scatter(data=data,x='x',y='y',s=6,c=data['hue'],cmap="Set1")
    axs[idx].set_title('Growing Degree-day {}'.format(year),**{'fontsize':'10'})
    axs[idx].axis('off')

plt.suptitle("Leaf position clusters ",**{'fontsize':'14','weight':'bold'})
plt.figtext(.51,0.95,'unsupervised cluster modeling with NMF + 2D-embedding with t-SNE:', **{'fontsize':'10','weight':'light'}, ha='center')

fig.legend(legend_list,topics,loc=(0.1,0.89),ncol=3)
plt.subplots_adjust(top=0.85)

plt.show()

#UMAP 3D
umap_3d = umap.UMAP(n_components=3, init='random', random_state=42)
proj_3d = umap_3d.fit_transform(X)

fig_3d = px.scatter_3d(
    proj_3d, x=0, y=1, z=2,
    color=y, labels={'color': 'position'}
)
fig_3d.update_traces(marker_size=5)
fig_3d.update_layout(title="UMAP visualization", 
                     xaxis={'visible': False, 'showticklabels': False},
                     yaxis={'visible': False, 'showticklabels': False})
    
fig_3d.show()







#Application de modèle (Régression si non colinéaire ou Arbre de décision)