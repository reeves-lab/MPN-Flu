from Code.Manual_gating_dict import manual_gates_dict
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

    df = pd.read_csv("../Data/230327_manual gating analysis.csv")
    df = df.rename(columns={'file': 'sample_name'})
    df['sample_name'] = [re.sub('.*?-(.*?) (.*?)_Hobbs.*?$', '\\1_\\2', i) for i in df['sample_name']]
    df[['Sample ID', 'Stimulation']] = df['sample_name'].str.split(pat='_', n=1, expand=True)
    df = df[df["Stimulation"] == 'Unstim'].drop(columns={'Stimulation'})
    df = df.drop(columns={'sample_name'})
    df = pd.merge(df, serology, left_on="Sample ID", right_index=True)
    return df


def plot_waterfall(df):
    col_dict = {'B': '#0944a8', 'I': '#008000', 'T': '#a11313'}
    if len(df) < 10:
        plt.figure(figsize=(4, 3))
    else:
        plt.figure(figsize=(6, 6))
    sns.barplot(x="Coef", y="pretty_name", data=df, hue='cell_type', dodge=False, palette=col_dict)
    plt.legend(handles=[], labels=[], frameon=False)
    plt.xlabel('change in LogOdds per unit SD')
    plt.ylabel('')
    plt.suptitle(f'{group}')
    plt.tight_layout()
    plt.savefig(f"../Graphs/Elastic Net/{today}_{group}_MG.jpg", dpi=300)
    plt.close()


today = datetime.now().date().strftime("%y%m%d")
startTime = datetime.now()
data = abundances_with_HAI()
data['HAI_All'] = data['HAI_MPN'].fillna(data['HAI_Healthy'])
total = data.filter(regex='total')
group = 'HAI_All'

np.random.seed(92487)

for group in ['HAI_MPN', 'HAI_Healthy']:
    # group = 'HAI_MPN'
    data = pd.merge(data[[group]], total, left_index=True, right_index=True)
    to_model = data.dropna(subset=[group])
    en_test = en.FeatureSelection(reg_name=group, run="test", alpha_th=(.01, 0.5), l1_th=(0.6, 0.8), id='MG')
    en_test.run_en(to_model=to_model)
    out = en_test.output

    fts_from_test = np.append(np.unique(out.Cluster.to_list()), group)
    en_valid = en.FeatureSelection(reg_name=group, run="validation", alpha_th=(.01, 0.5), l1_th=(0.6, 0.8), id='MG')
    en_valid.run_en(to_model=to_model[fts_from_test], test_score=en_test.best_model_score)

    coef_df = en_valid.output
    coef_df['Cluster'] = [re.sub('[|].*? :: (.*?)[|].*?$', '__\\1', i) for i in coef_df['Cluster']]
    coef_df['Cluster'] = [re.sub('[|]', '__', i) for i in coef_df['Cluster']]
    coef_df['Cluster'] = [re.sub('[|]', '__', i) for i in coef_df['Cluster']]
    coef_df['Cluster'] = [re.sub(' ', '_', i) for i in coef_df['Cluster']]
    coef_df['Cluster'] = [re.sub('[(](.*?)[)]', '\\1', i) for i in coef_df['Cluster']]

    coef_df[['Gate', 'Marker']] = coef_df['Cluster'].str.split(pat='__', n=1, expand=True)
    coef_df['Marker'] = [re.sub('^.*?_', '', i) for i in coef_df['Marker']]
    coef_df[['cell_type', 'pretty_name']] = coef_df['Gate'].map(manual_gates_dict).apply(pd.Series)
    coef_df['pretty_name'] = coef_df['pretty_name'] + ' (' + coef_df['Marker'] + ')'
    coef_df = coef_df.sort_values(by='Coef', ascending=False)

    plot_waterfall(df=coef_df)

with open(f'runtime.txt', 'a') as runtime_file:
    runtime_file.write(f'{datetime.now()}: {datetime.now() - startTime}\n   regression')
