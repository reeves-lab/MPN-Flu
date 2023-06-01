import Code.en_class_model_230310 as en
import matplotlib.pyplot as plt
from datetime import datetime
import seaborn as sns
import pandas as pd
import numpy as np
import os
import re
os.chdir('./Code')


def format_data():
    clinical = pd.read_excel("../Data/221027_metadata.xlsx", sheet_name="Combined", dtype={'Sample ID': str})
    clinical['Diagnosis'] = clinical['Diagnosis'].fillna('Healthy')
    clinical['Diagnosis'] = [re.sub('MPN', 'PV_ET', j) for j in clinical['Diagnosis']]
    df = pd.read_csv("../Data/omiqData_formatted/230324_cluster_abundances.csv")
    df[['Sample ID', 'Stimulation']] = df['sample_name'].str.split(pat='_', n=1, expand=True)
    df = df[df["Stimulation"] == 'Unstim'].drop(columns={'Stimulation', 'sample_name'})
    df = pd.merge(df, clinical[['Diagnosis', 'Sample ID']], on="Sample ID")
    comparisons_ = [["Healthy", "PV_ET"], ["Healthy", "MF"]]
    return df, comparisons_


def plot_waterfall(df):
    if len(df) < 10:
        plt.figure(figsize=(4, 3))
    else:
        plt.figure(figsize=(6, 4))
    col_dict = {'B': '#0944a8', 'I': '#008000', 'T': '#a11313'}
    df = df.sort_values(by='Coef', ascending=False)
    df['cell'] = [re.sub('^(.).*?$', '\\1', i) for i in df['Cluster']]
    df['Cell_type_w_cluster'] = df['Cell type'] + ' (' + df['Cluster'] + ')'
    sns.barplot(x="Coef", y="Cell_type_w_cluster", data=df, hue='cell', dodge=False, palette=col_dict)
    plt.legend(handles=[], labels=[], frameon=False)
    plt.xlabel('change in LogOdds per unit SD')
    plt.ylabel('')
    plt.suptitle(f'{name}')
    plt.tight_layout()
    plt.savefig(f"../Graphs/Elastic Net/{name}.jpg", dpi=300)
    plt.close()


today = datetime.now().date().strftime("%y%m%d")
startTime = datetime.now()
data, comparisons = format_data()
annotation = pd.read_excel('../Data/Cluster_phenotypes.xlsx')

for comparison in comparisons:
    # comparison = comparisons[0]
    name = "_".join(comparison)
    to_model = data.loc[data['Diagnosis'].isin(comparison)]
    to_model = to_model.drop(columns={'Sample ID'})

    en_test = en.FeatureSelection(run="test", classification='Diagnosis')
    en_test.run_en(to_model=to_model)
    out = en_test.output
    features_from_test = np.append(np.unique(out.Cluster.to_list()), 'Diagnosis')

    en_validation = en.FeatureSelection(run="validation", classification='Diagnosis')
    en_validation.run_en(to_model=to_model[features_from_test], test_score=en_test.best_model_score)

    coef_df = en_validation.output
    coef_df = pd.merge(coef_df, annotation[['Cluster', 'Cell type']], on='Cluster')
    coef_df = coef_df[coef_df['Cell type'] != '-']  # - marks clusters with fewer than 1000 cells
    plot_waterfall(df=coef_df)

with open(f'../Data/en_runtime.txt', 'a') as runtime_file:
    runtime_file.write(f'{datetime.now()}: {datetime.now() - startTime}\n')
