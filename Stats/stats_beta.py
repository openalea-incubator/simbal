import numpy as np
from numpy import corrcoef
import pandas as pd
import matplotlib.pyplot as plt
from sklearn import model_selection
from sklearn.datasets import make_classification
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, auc, accuracy_score, confusion_matrix
from sklearn.preprocessing import LabelEncoder
import seaborn as sns

#Analyse exploratoire

impulse = pd.read_csv('C:/Users/FARRA/AppData/Local/GitHubDesktop/app-3.2.3/src/caracol/Stats/IMPULSE_training.csv', sep=';', encoding='latin-1')
impulse_df = pd.DataFrame(data=impulse)
print(impulse_df)
impulse_df.info()

matrix = impulse.corr(method='pearson', min_periods=1, numeric_only=False)
print(matrix)


f, ax = plt.subplots(figsize=(12, 10))
sns.heatmap(matrix, ax=ax, cmap="YlGnBu", linewidths=0.1)
plt.show()


#Sélection des features sensibles

#Application de modèle (Régression si non colinéaire ou Arbre de décision)