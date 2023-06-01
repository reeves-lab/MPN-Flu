import statsmodels.stats.multicomp as mc
from statsmodels.formula.api import ols
import matplotlib.pyplot as plt
from datetime import datetime
from tabulate import tabulate
import statsmodels.api as sm
from functools import reduce
from scipy import cluster
from scipy import stats
import seaborn as sns
import pandas as pd
import numpy as np
import re

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


def format_hai(path):
    # path = './Data/221110_HAI_titers.xlsx'
    data = pd.read_excel(path, sheet_name='formatted', skiprows=1, dtype={'Sample ID': str})
    data = data.melt(id_vars=['idx', 'Diagnosis', 'Sample ID', 'Date', 'vaccine_year'])
    data[['strain', 'replicate']] = data['variable'].str.split(pat='-', n=1, expand=True)
    data['Diagnosis'] = [re.sub('MPN', 'PV_ET', i) for i in data['Diagnosis']]
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
    data['group'] = [re.sub('PV_ET|MF', 'MPN', i) for i in data['Diagnosis']]
    return data


def format_anova_tbl(aov):
    """
    Calculates eta and omega squared and add them as columns to anova table
    :param aov: anova table
    :return: formatted anova table
    """
    aov['mean_sq'] = aov[:]['sum_sq'] / aov[:]['df']
    aov['eta_sq'] = aov[:-1]['sum_sq'] / sum(aov['sum_sq'])
    aov['omega_sq'] = (aov[:-1]['sum_sq'] - (aov[:-1]['df'] * aov['mean_sq'][-1])) / (
            sum(aov['sum_sq']) + aov['mean_sq'][-1])

    columns = ['sum_sq', 'df', 'mean_sq', 'F', 'PR(>F)', 'eta_sq', 'omega_sq']
    aov = aov[columns]
    return aov


def format_spearmanr_output(df, suffix):
    df = pd.DataFrame(df, columns=names, index=names)
    df = df.iloc[0:cl_len, cl_len:len(df)]
    df.columns = [re.sub('$', suffix, col) for col in df.columns]
    return df


def plot_sub(year, n):
    # n = 3
    # year = '2018-2019'
    strain_dict = {
        'California': 'California/7/2009 (H1N1)',
        'Michigan': 'Michigan/45/2015 (H1N1)',
        'HK': 'Hong Kong/4801/2014 (H3N2)',
        'Singapore': 'Singapore/INFIMH-16-0019/2016 (H3N2)'
    }
    plt.subplot(1, 3, n)
    per_year = HAI_means[HAI_means['vaccine_year'] == year]
    per_year = per_year.replace({'strain': strain_dict})
    if year == '2017-2018':
        per_year['strain'] = pd.Categorical(per_year['strain'], ordered=True,
                                            categories=['Michigan/45/2015 (H1N1)', 'Hong Kong/4801/2014 (H3N2)'])
        per_year = per_year.sort_values(by='strain')
    ax = sns.lineplot(data=per_year, x='strain', y="value", hue='Diagnosis', units="Sample ID",
                      estimator=None, palette=['#ed4135', '#252160', '#05cfc2'], alpha=0.5)
    sns.stripplot(data=per_year, x='strain', y="value", hue='Diagnosis', palette=['#ed4135', '#252160', '#05cfc2'],
                  jitter=True, legend=False)
    sns.despine()
    plt.xlabel('')
    plt.ylabel('')
    plt.ylim(0, 70)
    plt.title(f'{year}')
    plt.xticks(rotation=30)

    if n == 1:
        plt.ylabel('On-year HAI titer')
    if n == 3:
        plt.legend(frameon=False)
    else:
        ax.get_legend().remove()


# TODO 1: per-sample HAI visualisation and ANOVA
today = datetime.now().strftime('%y%m%d')
hai = format_hai(path='./Data/221110_HAI_titers.xlsx')

# get means per sample (from the two on-year strains)
cols = ['Sample ID', 'group', 'Diagnosis', 'vaccine_year']
HAI_means = hai[cols + ['value']].groupby(cols).mean().reset_index()
HAI_means['group'] = pd.Categorical(HAI_means['group'], ordered=True, categories=['Healthy', 'MPN'])
HAI_means['Diagnosis'] = pd.Categorical(HAI_means['Diagnosis'], ordered=True, categories=['Healthy', 'PV_ET', 'MF'])
HAI_means.to_csv('Data/230326_formatted_HAI_means.csv', index=False)

# plot boxplot of HAI in groups with stripplot of HAI in diagnosis
plt.figure(figsize=(3.6, 5))
sns.set(font_scale=1.3, style='ticks')
ax = sns.boxplot(data=HAI_means, x='group', y="value", linewidth=1.2, width=0.5,
                 palette=['#ed4135', '#196285'])
for patch in ax.patches:
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r, g, b, .6))
# get diagnosis-specific dataframes for stripplot
Healthy = HAI_means[HAI_means['Diagnosis'] == 'Healthy']
PVET = HAI_means[HAI_means['Diagnosis'] == 'PV_ET']
MF = HAI_means[HAI_means['Diagnosis'] == 'MF']
# plot stripplots
sns.stripplot(data=Healthy, x='group', y="value", s=8, color='#ed4135', jitter=0.2, marker='D')
sns.stripplot(data=PVET, x='group', y="value", s=8, color='#252160', jitter=0.2, marker='o')
sns.stripplot(data=MF, x='group', y="value", s=10, color='#05cfc2', jitter=0.2, marker='^')
sns.despine()
plt.xlabel('')
plt.ylabel('On-year HAI titer')
lgd = plt.legend(frameon=False, bbox_to_anchor=(1, 0.5))
plt.savefig('./Graphs/HAI_boxplot.jpeg', dpi=300, bbox_extra_artists=(lgd,), bbox_inches='tight')

# run ANOVA and post-hoc test
model = ols('value ~ C(Diagnosis)', data=HAI_means).fit()
# model = ols('value ~ C(group)', data=HAI_means).fit()
aov_table = sm.stats.anova_lm(model, typ=2)
aov_table = format_anova_tbl(aov_table)
comp = mc.MultiComparison(HAI_means['value'], HAI_means['Diagnosis'])
post_hoc_res = comp.tukeyhsd()
tukey_tbl = post_hoc_res.summary()
print(tukey_tbl)
with open('./Data/results/HAI_ANOVA.txt', 'a') as file:
    aov_table = pd.DataFrame(np.vstack([aov_table.columns, aov_table]))
    file.write(f'ANOVA TABLE:\n{tabulate(aov_table)}\n\nPOSTHOC:\n{tabulate(tukey_tbl.data)}\n')


# TODO 2: per-sample per strain HAI visualisation
# get means per sample per strain
cols = ['Sample ID', 'group', 'Diagnosis', 'vaccine_year', 'strain']
HAI_means = hai[cols + ['value']].groupby(cols).mean().reset_index()
HAI_means['group'] = pd.Categorical(HAI_means['group'], ordered=True, categories=['Healthy', 'Disease'])
HAI_means['Diagnosis'] = pd.Categorical(HAI_means['Diagnosis'], ordered=True, categories=['Healthy', 'MPN', 'MF'])
# plot
plt.figure(figsize=(12, 7))
sns.set(font_scale=1, style='ticks')
plot_sub('2016-2017', n=1)
plot_sub('2017-2018', n=2)
plot_sub('2018-2019', n=3)
plt.tight_layout()
plt.savefig('./Graphs/HAI/HAI_boxplot_strains.jpeg', dpi=300)


# TODO 3: HAI correlation to cluster abundance
HAI_means = HAI_means.pivot(index=['Sample ID', 'Diagnosis'], columns='group', values='value').reset_index()
HAI_means = HAI_means.rename(columns={'Disease': 'HAI_Disease', 'Healthy': 'HAI_Healthy'})
HAI_means['HAI_All'] = HAI_means['HAI_Disease'].fillna(HAI_means['HAI_Healthy'])
# read in cluster abundance dataframe to correlate
abund = pd.read_csv('./Data/omiqData_formatted/230101_cluster_abundances_subs_corr.csv')
abund['Sample ID'] = [re.sub('^([0-9]*?)_.*?$', '\\1', i) for i in abund['sample_name']]
abund = abund[abund['sample_name'].str.contains('Unstim')].drop(columns='sample_name')
HAI_w_abund = pd.merge(HAI_means, abund, on='Sample ID', how='left').set_index('Sample ID')

clusters = HAI_w_abund.filter(regex='.*Cluster.*')
cl_len = len(clusters.columns)
HAI_cols = ['HAI_Healthy', 'HAI_Disease', 'HAI_All']
names = clusters.columns.tolist() + HAI_cols
# spearmanr correlation
rho, pval = stats.spearmanr(clusters, HAI_w_abund[HAI_cols], nan_policy='omit')
rho = format_spearmanr_output(rho, '_coef')
pval = format_spearmanr_output(pval, '_pval')
abund_means = HAI_w_abund.drop(columns=HAI_cols).groupby('Diagnosis').agg(func='mean').transpose()

# merge mean abundance in diagnosis with correlation coefs and pvalues
merged = reduce(
    lambda left, right: pd.merge(left, right, left_index=True, right_index=True),
    [abund_means, rho, pval])
merged['cell'] = [re.sub('^(.) Cl.*?$', '\\1', i) for i in merged.index]

# get clusters with significant correlation coefs
melt = merged.drop(columns=['Healthy', 'MPN', 'MF', 'cell']).reset_index().melt(id_vars='index')
melt[['hai', 'corr_to', 'val_type']] = melt['variable'].str.split(n=3, pat='_', expand=True)
melt = melt.drop(columns=['hai', 'variable'])
melt = melt.pivot(columns='val_type', values='value', index=['index', 'corr_to']).reset_index()
melt = melt[melt['pval'] < 0.05]
melt.to_csv('./Data/for_circos_significant_clusters.csv')


# TODO 4: Export cell-specific dataframe for the circos script
cell = 'T'
per_cell = merged.filter(regex=f'{cell}.*?', axis=0)
ax = sns.clustermap(per_cell.filter(like='coef'), cmap='coolwarm',
                    method="ward", yticklabels=1, xticklabels=1, figsize=(6, 6))
plt.savefig(f'./Graphs/HAI/heatmap_{cell}.jpeg', dpi=300)

Z = cluster.hierarchy.ward(np.array(per_cell.filter(like='coef')))
cutree = pd.DataFrame(cluster.hierarchy.cut_tree(Z, n_clusters=2), columns=['id'])
per_cell = pd.merge(per_cell.reset_index(), cutree, left_index=True, right_index=True).set_index('index')
per_cell.to_csv(f'./Data/for_circos_{cell}.csv')
