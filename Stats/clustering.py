import matplotlib.pyplot as plt
from sklearn.datasets import make_moons
from sklearn.metrics import pairwise_distances_argmin
from sklearn.metrics import pairwise_distances
import numpy as np


X, y = make_moons(n_samples=100, random_state=123)

sc = X.reshape(-1)
print(sc)

def find_clusters(sc, n_clusters, rseed=3,metric='euclidean'):
    # 1. Randomly choose clusters
    rng = np.random.RandomState(rseed)
    i = rng.permutation(X.shape[0])[:n_clusters]
    centers = sc[i]
    n = 0
    labels_iter = []
    centers_iter = [centers]
    while True:
        # 2a. Assign labels based on closest center
        labels = pairwise_distances_argmin(sc, centers, metric=metric)
        labels_iter.append(labels)

        # 2b. Find new centers from means of points
        new_centers = np.array([sc[labels == i].mean(0)
                                for i in range(n_clusters)])

        # 2c. Check for convergence
        if np.all(centers == new_centers):
            break
        centers = new_centers
        centers_iter.append(centers)

    return centers_iter, labels_iter

centers_iter, labels_iter = find_clusters(sc, 4,rseed=3,metric='euclidean')
plt.scatter(sc[:, 0], X[:, 1], c=labels_iter[-1], s=50, cmap='viridis')
plt.scatter(centers_iter[-1][:, 0], centers_iter[-1][:, 1], c='purple', s=200)

