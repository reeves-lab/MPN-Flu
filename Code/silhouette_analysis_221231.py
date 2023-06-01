import pandas as pd
import numpy as np
import re
from sklearn.metrics import silhouette_score

# dat = pd.read_csv('./Data/omiqData/data_export_I_15x10.csv')
dat = pd.read_csv('./Data/omiqData/data_export_I_10x10.csv')
dat['OmiqFilter'] = [int(re.sub('cl-', '', i)) for i in dat['OmiqFilter']]
dat = dat.sample(n=400000)

X = dat.drop(columns={'In115Di :: CD3', 'Cd111Di :: CD19', 'OmiqFilter', 'OmiqFileIndex'})
X = np.array(X)
cluster_labels = np.array(dat['OmiqFilter'])
silhouette_avg = silhouette_score(X, cluster_labels)


