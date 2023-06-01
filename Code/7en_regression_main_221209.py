import Code.en_regression_model_221209 as en
import matplotlib.pyplot as plt
from datetime import datetime
import seaborn as sns
import pandas as pd
import numpy as np
import re
import os
os.chdir('./Code')


def abundances_with_HAI():
    """
    Add UNSTIM cluster abundances and metadata to formatted HAI titer data
    :return: Formatted abundances with HAI and Diagnosis column
    """
    serology = pd.read_csv('../Data/230326_formatted_HAI_means.csv', dtype={'Sample ID': str})
    serology = serology.drop(columns={'Diagnosis', 'vaccine_year'})
    serology = serology.pivot(columns='group', index='Sample ID')
    serology.columns = [re.sub('^', 'HAI_', i[1]) for i in serology.columns]

    df = pd.read_csv("../Data/omiqData_formatted/230324_cluster_abundances.csv")
    # df = pd.read_csv("./Data/230321_manual gating analysis.csv")
    df['sample_name'] = [re.sub('.*?-(.*?) (.*?)_Hobbs.*?$', '\\1 \\2', i) for i in df['sample_name']]
    df[['Sample ID', 'Stimulation']] = df['sample_name'].str.split(pat='_', n=1, expand=True)
    df = df[df["Stimulation"] == 'Unstim'].drop(columns={'Stimulation'})
    df = df.drop(columns={'sample_name'})
    df = pd.merge(df, serology, left_on="Sample ID", right_index=True)
    return df


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
    plt.suptitle(f'{group}')
    plt.tight_layout()
    plt.savefig(f"../Graphs/Elastic Net/{today}_{group}.jpg", dpi=300)
    plt.close()


today = datetime.now().date().strftime("%y%m%d")
dir_name = f'{today}_regression'
startTime = datetime.now()
data = abundances_with_HAI()
annotation = pd.read_excel('../Data/Cluster_phenotypes.xlsx')

np.random.seed(92487)

for group in ['HAI_MPN', 'HAI_Healthy']:
    # group = 'HAI_MPN'
    to_model = data.dropna(subset=[group])
    en_test = en.FeatureSelection(reg_name=group, run="test", alpha_th=(.001, 1), l1_th=(.4, 0.6))
    en_test.run_en(to_model=to_model)
    out = en_test.output

    fts_from_test = np.append(np.unique(out.Cluster.to_list()), group)
    en_valid = en.FeatureSelection(reg_name=group, run="validation", alpha_th=(.001, 1), l1_th=(.4, 0.6))
    en_valid.run_en(to_model=to_model[fts_from_test], test_score=en_test.best_model_score)

    coef_df = en_valid.output
    coef_df = pd.merge(coef_df, annotation[['Cluster', 'Cell type']], on='Cluster')
    coef_df = coef_df[coef_df['Cell type'] != '-']  # - marks clusters with fewer than 1000 cells
    plot_waterfall(df=coef_df)

with open(f'runtime.txt', 'a') as runtime_file:
    runtime_file.write(f'{datetime.now()}: {datetime.now() - startTime}\n   regression')
