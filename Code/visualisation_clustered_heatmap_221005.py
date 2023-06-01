from sklearn.preprocessing import MinMaxScaler
import matplotlib.pyplot as plt
from datetime import datetime
import seaborn as sns
import pandas as pd
import re

today = datetime.now().strftime("%y%m%d")
phenotypes = pd.read_csv("Data/omiqData_formatted/230324_cluster_phenotypes.csv")
phenotypes = phenotypes.drop(columns={'sample_name'}).groupby('cluster').median().reset_index()
annotation = pd.read_excel('Data/Cluster_phenotypes.xlsx')

for cell in ['B', 'I', 'T']:
    # cell = 'T'
    per = phenotypes[phenotypes['cluster'].str.contains(cell)].set_index('cluster')
    per.columns = [re.sub(' median', '', i) for i in per.columns]
    per.columns = [re.sub('^.*?_(.*?$)', '\\1', i) for i in per.columns]
    scaler = MinMaxScaler()
    per_scl = scaler.fit_transform(per)
    per_scl = pd.DataFrame(per_scl, columns=per.columns, index=per.index)
    per_scl = pd.merge(annotation[['Cluster', 'Cell type']], per_scl, left_on='Cluster', right_index=True)
    per_scl = per_scl.sort_values(by='Cell type')
    per_scl['Cell type with cluster'] = per_scl['Cluster'] + ' (' + per_scl['Cell type'] + ')'
    per_scl = per_scl.drop(columns={'Cluster', 'Cell type'}).set_index('Cell type with cluster')

    sns.clustermap(per_scl, cmap='viridis', yticklabels=True, xticklabels=True,
                   col_cluster=False, row_cluster=False, figsize=(10, 7))
    plt.xlabel('')
    plt.ylabel('')
    plt.tight_layout()
    plt.savefig(f'Graphs/clustered_heatmap_{cell}.jpeg', dpi=300)
    plt.close()
