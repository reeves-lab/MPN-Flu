import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import glob
import re


key_word = 'HAI'  # or HAI
files = glob.glob(f'./Data/results/Boruta/*{key_word}*')
files.extend(glob.glob(f'./Data/results/Elastic Net/*validation*{key_word}*'))
# dac = pd.read_csv('./Data/results/DAC_results.csv')
annotation = pd.read_excel('./Data/Cluster_phenotypes.xlsx')

li = ['HAI_Healthy', 'HAI_MPN']
# li = ['PV_ET', 'MF']
for group in li:
    # group = 'HAI_Healthy'
    group_files = [i for i in files if group in i]
    boruta = pd.read_excel(f'{[i for i in group_files if "Boruta" in i][0]}')
    en = pd.read_excel(f'{[i for i in group_files if "Elastic" in i][0]}')

    # subset features are common in Boruta and EN
    robust_df = en[en['Cluster'].isin(boruta['Cluster'])].copy()
    robust_df['Cluster'] = [re.sub(' Cluster ', '', i) for i in robust_df['Cluster']]

    # merge with phenotypes dataframe
    robust_df = pd.merge(left=robust_df, right=annotation[['Cluster', 'Count', 'Cell type']], on='Cluster')
    robust_df['Cell type with cluster'] = robust_df['Cell type'] + ' (' + robust_df['Cluster'] + ')'
    robust_df = robust_df[robust_df['Count'] > 1000]
    robust_df['cell'] = [re.sub('^(.).*?$', '\\1', i) for i in robust_df['Cluster']]
    robust_df = robust_df[robust_df['Cell type'] != '-']  # - marks clusters with fewer than 1000 cells
    robust_df = robust_df[~robust_df['Cell type'].str.contains('Unassigned')]

    # plot
    col_dict = {'B': '#0944a8', 'I': '#008000', 'T': '#a11313'}
    plt.figure(figsize=(5, 5))
    sns.barplot(x="Coef", y="Cell type with cluster", data=robust_df, hue='cell', dodge=False, palette=col_dict)
    plt.legend(handles=[], labels=[], frameon=False)
    plt.xlabel('change in LogOdds per unit SD')
    plt.ylabel('')
    plt.suptitle(f'{group}')
    plt.tight_layout()
    plt.savefig(f"./Graphs/Consistent_features_{group}.jpg", dpi=300)
    plt.close()

