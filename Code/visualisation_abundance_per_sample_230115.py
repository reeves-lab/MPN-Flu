import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import re


def format_hai():
    data = pd.read_excel('./Data/221110_HAI_titers.xlsx', sheet_name='formatted', skiprows=1, dtype={'Sample ID': str})
    data = data.melt(id_vars=['idx', 'Diagnosis', 'Sample ID', 'Date', 'vaccine_year'])
    data[['strain', 'replicate']] = data['variable'].str.split(pat='-', n=1, expand=True)
    data['Diagnosis'] = [re.sub('MPN', 'PV/ET', i) for i in data['Diagnosis']]
    data.drop(columns=['idx', 'variable', 'Date'], inplace=True)
    strain_dict = {'2016-2017': {'on-year': ['California', 'HK'], 'off-year': ['Michigan', 'Singapore']},
                   '2017-2018': {'on-year': ['Michigan', 'HK'], 'off-year': ['California', 'Singapore']},
                   '2018-2019': {'on-year': ['Michigan', 'Singapore'], 'off-year': ['California', 'HK']}}

    for year in data['vaccine_year'].unique():
        mask_year = data['vaccine_year'] == year
        for on_off in ['on-year', 'off-year']:
            mask_strains = data['strain'].isin(strain_dict[year][on_off])
            data.loc[mask_year & mask_strains, 'immunity'] = on_off

    data = data[data['immunity'] == 'on-year']
    return data


meta = pd.read_excel("./Data/221027_metadata.xlsx", sheet_name="Combined", dtype={'Sample ID': str})
meta['Diagnosis'] = meta['Diagnosis'].fillna('Healthy')
meta['Diagnosis'] = [re.sub('MPN', 'PV/ET', i) for i in meta['Diagnosis']]

hai = format_hai()[['Sample ID', 'value']]
hai = pd.merge(hai, meta[['Sample ID', 'Diagnosis']], on='Sample ID', how='right')
hai['Diagnosis'] = pd.Categorical(hai['Diagnosis'], categories=['Healthy', 'PV/ET', 'MF'], ordered=True)
hai = hai.sort_values(by=['Diagnosis', 'value']).reset_index(drop=True)
hai_idx = list(hai['Sample ID'].unique())

df = pd.read_csv('./Data/omiqData_formatted/230324_cluster_abundances.csv')
df[['Sample ID', 'Stim']] = df['sample_name'].str.split(pat='_', n=1, expand=True)
df = df[df['Stim'] == 'Unstim']
df = pd.merge(df, meta[['Diagnosis', 'Sample ID']], on='Sample ID')
df = df.drop(columns={'sample_name', 'Stim'}).melt(id_vars=['Diagnosis', 'Sample ID'])
df['Diagnosis'] = pd.Categorical(df['Diagnosis'], categories=['Healthy', 'PV/ET', 'MF'], ordered=True)
df['Sample ID'] = pd.Categorical(df['Sample ID'], categories=hai_idx, ordered=True)
df = df.sort_values(by=['Diagnosis', 'Sample ID', 'variable'])


pal = ['#ed4135', '#252160', '#05cfc2']
for cluster in df['variable'].unique():
    # cluster='B01'
    per_cluster = df[(df['variable'] == cluster)].copy()

    plt.figure(figsize=(12, 10))
    plt.subplot(2, 1, 1)
    ax = sns.boxplot(data=hai, x='Sample ID', y="value", linewidth=0.6, hue='Diagnosis', palette=pal, dodge=False)
    for patch in ax.patches:
        r, g, b, a = patch.get_facecolor()
        patch.set_facecolor((r, g, b, .6))
    sns.stripplot(data=hai, x='Sample ID', y="value", s=5,  hue='Diagnosis', palette=pal, jitter=0.2)
    plt.xlabel('')
    plt.ylabel('HAI titer')
    plt.legend(frameon=False, loc='center left', bbox_to_anchor=(1, 0.5))

    plt.subplot(2, 1, 2)
    # plt.figure(figsize=(12, 5))
    sns.barplot(data=per_cluster, x='Sample ID', y='value', hue='Diagnosis', palette=pal, dodge=False)
    plt.xlabel('')
    plt.ylabel('')
    plt.legend(frameon=False, loc='center left', bbox_to_anchor=(1, 0.5))
    plt.tight_layout()
    plt.suptitle(f'{cluster}', y=1.005)
    plt.savefig(f'./Graphs/Abundance per sample/{cluster}.jpg', dpi=300, bbox_inches='tight')
    plt.close()
