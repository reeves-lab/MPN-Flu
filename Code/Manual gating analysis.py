from statannotations.Annotator import Annotator
import statsmodels.stats.multicomp as mc
from statsmodels.formula.api import ols
from Code.Manual_gating_dict import *
import matplotlib.pyplot as plt
from datetime import datetime
import statsmodels.api as sm
import seaborn as sns
import pandas as pd
import re

today = datetime.today().strftime('%y%m%d')


# TODO 1: format manually gated features (marker medians and gate abundances)
# read in clinical data
clinical = pd.read_excel("Data/221027_metadata.xlsx", sheet_name="Combined", dtype={'Sample ID': str})
clinical['Diagnosis'] = clinical['Diagnosis'].fillna('Healthy')
clinical['Diagnosis'] = [re.sub('MPN', 'PV/ET', i) for i in clinical['Diagnosis']]
# tidy up omiq data
df = pd.read_csv("Data/230327_manual gating analysis.csv")
df.columns = [re.sub('[|].*? :: (.*?)[|].*?$', '__\\1', i) for i in df.columns]
df.columns = [re.sub('[|]', '__', i) for i in df.columns]
df.columns = [re.sub(' ', '_', i) for i in df.columns]
df.columns = [re.sub('-', '_', i) for i in df.columns]
df.columns = [re.sub('[(](.*?)[)]', '\\1', i) for i in df.columns]
df['file'] = [re.sub('.*?-(.*?) (.*?)_Hobbs.*?$', '\\1 \\2', i) for i in df['file']]
df[['Sample ID', 'Stimulation']] = df['file'].str.split(pat=' ', n=1, expand=True)
df = df[df["Stimulation"] == 'Unstim']
# merge with clinical
df = pd.merge(df, clinical[['Diagnosis', 'Sample ID']], on="Sample ID")
df = df.drop(columns={'file', 'Sample ID', 'Stimulation'})


# TODO 2: statistical analysis: ANOVA and Tukey HSD on each feature
lst = list(df.columns)
lst.remove('Diagnosis')
aov_res = {'feature': [], 'pvalue': [], 'MC res': []}
for ft in lst:
    # ft = 'Effector_helper_T__CD197_CCR7'
    mod = ols('{} ~ Diagnosis'.format(ft), data=df).fit()
    aov_table = sm.stats.anova_lm(mod, typ=2)
    p = aov_table['PR(>F)'].values[0]
    if p < 0.05:
        comp = mc.MultiComparison(df[ft], df['Diagnosis'])
        post_hoc_res = comp.tukeyhsd()
        tukey_tbl = pd.DataFrame(post_hoc_res.summary().data)
        tukey_tbl.columns = [str(i) for i in tukey_tbl.iloc[0]]
        tukey_tbl.drop(tukey_tbl.index[0], inplace=True)
        sig = tukey_tbl[tukey_tbl['reject']].copy().reset_index(drop=True)
        sig['groups'] = sig['group1'] + '_vs_' + sig['group2']
        MC_res = {}
        for j in range(len(sig)):
            MC_res[f'{sig["groups"][j]}'] = f'{sig["p-adj"][j]}'

        aov_res['feature'].append(ft)
        aov_res['pvalue'].append(p)
        aov_res['MC res'].append(MC_res)

# TODO 3: format ANOVA results
res = pd.DataFrame.from_dict(aov_res, orient='columns')
res[['gate', 'marker']] = res['feature'].str.split('__', n=1, expand=True)
# remove rows where Tukey didn't yield any H0 rejections
non_zero_di_index = [i for i in range(len(res)) if len(res['MC res'][i]) > 0]
res = res.loc[non_zero_di_index]
# explode dictionary column with Tukey's multiple comp results
di_df = res['MC res'].apply(pd.Series)
res = pd.merge(res, di_df, left_index=True, right_index=True)
# removing unnecessary columns
res = res.drop(columns={'MC res', 'pvalue'})
# melt by dictionary comparisons and remove nan values
res = res.melt(id_vars=['feature', 'gate', 'marker']).dropna()
res = res.rename(columns={'variable': 'comparison'})

# TODO 4: Benjamini-Hochberg correction
# 1. order by p-value column
# 2. calculate BH critical value for each p-value: rank (l)/ n of comparisons (m) * FDR (Q)
# 3. p-values lower than their critical threshold are significant
res = res.sort_values(by='value').reset_index(drop=True)
res['q-val'] = [(i + 1) / 850 * 0.2 for i in range(len(res))]
res_s = res[res['value'].astype(float) < res['q-val']].copy()
res_s[['cell', 'pretty_name']] = res_s['gate'].map(manual_gates_dict).apply(pd.Series)
better_order = ['cell', 'gate', 'pretty_name', 'marker', 'feature', 'comparison', 'value', 'q-val']
res_s = res_s[better_order].sort_values(by=['marker', 'cell', 'comparison']).reset_index(drop=True)
res_s['marker'] = [re.sub('^.*?_', '', i) for i in res_s['marker']]

res_s_B = res_s[res_s['cell'] == 'B'].sort_values(by=['marker', 'gate'])
res_s_T = res_s[res_s['cell'] == 'T'].sort_values(by=['marker', 'gate', 'comparison'])
res_s_I = res_s[res_s['cell'] == 'I'].sort_values(by=['marker', 'gate', 'comparison'])

di_chosen_markers = {
    'B': ['CXCR5', 'CXCR3', 'CCR7', 'total'],
    'I': ['CD11b', 'CCR2', 'CD57', 'DR', 'CD45RA', 'CD123', 'total'],
    'T': ['CXCR3', 'CCR2', 'CD27', 'CD44', 'total']
}

manual_gating_res = pd.DataFrame()
for cell in ['B', 'I', 'T']:
    per = res_s[(res_s['cell'] == cell) & (res_s['marker'].isin(di_chosen_markers[cell]))].copy()
    per = per.drop(columns={'feature', 'gate'})
    manual_gating_res = pd.concat([manual_gating_res, per])

manual_gating_res = manual_gating_res.sort_values(by=['cell', 'pretty_name', 'marker', 'comparison'])
manual_gating_res = manual_gating_res.rename(columns={'pretty_name': 'gated population'})

manual_gating_res.to_csv('manual_gating_stats_results.csv', index=False)


# TODO 5: Visualisation
melt = df.melt(id_vars='Diagnosis')
melt[['gate', 'marker']] = melt['variable'].str.split('__', n=1, expand=True)
melt = melt[['Diagnosis', 'gate', 'marker', 'value']]
melt['marker'] = [re.sub('^.*?_', '', i) for i in melt['marker']]
melt[['cell', 'pretty_name']] = melt['gate'].map(manual_gates_dict).apply(pd.Series)

cell = 'I'
ft = 'CD57'
# ylab = '% of total cells'
ylab = 'marker intensity'
per = melt[(melt['marker'] == ft) & (melt['cell'] == cell)].copy()
cell_ord = gate_order[cell]['order'].copy()
# gates_to_remove = ['Effector helper T', 'Effector cytotoxic T', 'DN T', 'NK T']
# cell_ord = [i for i in cell_ord if i not in gates_to_remove]

# gates_to_keep = ['CL mono', 'INT mono', 'NC mono', 'DC']
gates_to_keep = ['NK (56+16+)', 'NK (56+16-)', 'NK (56++16-)']
cell_ord = [i for i in cell_ord if i in gates_to_keep]
# cell_ord.insert(0, 'Innate cells')
per = per[per['pretty_name'].isin(cell_ord)]
# cell_ord.insert(0, gate_order[cell]['remove'])
per['pretty_name'] = pd.Categorical(per['pretty_name'], ordered=True, categories=cell_ord)

sns.set(font_scale=1.2, style='ticks', rc={'figure.figsize': (4.5, 5)})
plotting_parameters = {
    'data': per,
    'x': 'pretty_name',
    'y': 'value',
    'palette': ['#ed4135', '#252160', '#05cfc2'],
    'hue': 'Diagnosis'
}
ax = sns.boxplot(**plotting_parameters, linewidth=1.2, fliersize=0)
for patch in ax.patches:
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r, g, b, .6))
sns.stripplot(**plotting_parameters, dodge=True, s=5, linewidth=0.2, edgecolor='w')
plt.title('')
plt.ylabel(ylab)
plt.xlabel('')
plt.xticks(rotation=40)
# plt.ylim(-1, 21)
lgd = plt.legend(title='', frameon=False, bbox_to_anchor=(1, 0.5))
plt.tight_layout()
sns.despine()
plt.savefig(f"Graphs/Manual gating analysis/{cell}_{ft}.png", dpi=250, bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.close()


# PLOTTING BY GATE
gate = 'Lin_neg'
# per = melt[(melt['gate'] == gate) & (melt['marker'].isin(res_s_I[res_s_I['gate'] == gate]['marker']))].copy()
per = melt[(melt['gate'] == gate) & (melt['marker'].isin(['CD45RA', 'DR', 'CD123']))].copy()
sns.set(font_scale=1.2, style='ticks', rc={'figure.figsize': (4.5, 5)})
plotting_parameters = {
    'data': per,
    'x': 'marker',
    'y': 'value',
    'palette': ['#ed4135', '#252160', '#05cfc2'],
    'hue': 'Diagnosis'
}
ax = sns.boxplot(**plotting_parameters, linewidth=1.2, fliersize=0)
for patch in ax.patches:
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r, g, b, .6))
sns.stripplot(**plotting_parameters, dodge=True, s=5, linewidth=0.2, edgecolor='w')
plt.title('')
plt.ylabel(ylab)
plt.xlabel('')
plt.xticks(rotation=40)
# plt.ylim(-2, 50)
lgd = plt.legend(title='', frameon=False, bbox_to_anchor=(1, 0.5))
plt.tight_layout()
sns.despine()
plt.savefig(f"Graphs/Manual gating analysis/{gate}.png", dpi=250, bbox_extra_artists=(lgd,), bbox_inches='tight')
plt.close()





# marker = 'CCR7'
# cell_type = 'T'
#
# for j in melt[graph_of].unique():
#     with sns.plotting_context("notebook"):
#         # metric = 'CD185_CXCR5'
#         per = melt[melt[graph_of] == j]
#         relevant_vars = per[variable].unique()
#         size = (4 + (len(relevant_vars) - 1) * 1.5, 6)
#         fig, ax = plt.subplots(1, 1, figsize=size)
#         plotting_parameters = {
#             'data': per,
#             'x': variable,
#             'y': 'value',
#             'palette': ['#ed4135', '#252160', '#05cfc2'],
#             'hue': 'Diagnosis'
#         }
#
#         pairs = []
#         for x in relevant_vars:
#             li = [
#                 [(x, 'Healthy'), (x, 'PV/ET')],
#                 [(x, 'Healthy'), (x, 'MF')],
#                 [(x, 'PV/ET'), (x, 'MF')]
#             ]
#             pairs.extend(li)
#
#         sns.boxplot(**plotting_parameters)
#         plt.title(re.sub('_', ' ', j))
#         plt.legend(title='', frameon=False)
#         plt.xlabel('')
#         plt.tight_layout()
#
#         annotator = Annotator(ax, pairs, **plotting_parameters)
#         annotator.configure(test="Mann-Whitney")
#         _, corrected_results = annotator.apply_and_annotate()
#
#         plt.savefig(f"Graphs/Manual gating analysis/{j}.png", dpi=300)
#         plt.close()


