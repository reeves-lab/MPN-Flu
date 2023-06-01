import plotly.graph_objects as go
import matplotlib.pyplot as plt
import plotly.offline as pyo
from functools import reduce
import seaborn as sns
import pandas as pd
import numpy as np
import math
import re


def format_data(counts_path, perc_path, as_percent_total):
    """
    Formats subsampling corrected counts dataframe, by calculating percentage total of each cluster
    using the perc_path dataframe. This step can be skipped. Next it merges the counts with metadata
    - diagnosis and stimulation information.
    """
    # counts_path = './Data/omiqData_formatted/230101_cluster_abundances_subs_corr.csv'
    dat = pd.read_csv(counts_path)
    dat[['Sample ID', 'stimulation']] = dat['sample_name'].str.split(n=2, pat='_', expand=True)
    dat = dat[dat['stimulation'] == 'Unstim'].drop(columns=['stimulation', 'sample_name'])
    dat.set_index('Sample ID', inplace=True)

    # perc_path = './Data/230108_gated_populations_percentages.csv'
    perc_total = pd.read_csv(perc_path)
    cols = ['Sample ID', 'variable', 'Diagnosis', 'normalised to sample']
    perc_total = perc_total[cols].pivot(columns='variable', index=['Sample ID', 'Diagnosis'])
    perc_total.columns = perc_total.columns.droplevel(0)

    cluster_n = {'B': 22, 'Innate': 17, 'T': 32}
    perc_frames = []
    for cell_type in ['B', 'Innate', 'T']:
        per_cell = perc_total[[f'{cell_type} cells']]
        # per_cell.set_index("sample_name", inplace=True)
        per_cell = per_cell[[col for col in per_cell.columns for i in range(cluster_n[cell_type])]]
        perc_frames.append(per_cell)
    perc_total = reduce(lambda left, right: pd.merge(left, right, left_index=True, right_index=True), perc_frames)
    perc_total = perc_total.reset_index()
    perc_total.set_index('Sample ID', inplace=True)
    perc_total.index = perc_total.index.astype('str')
    diag = perc_total[['Diagnosis']]
    perc_total.drop(columns='Diagnosis', inplace=True)
    if as_percent_total:
        perc_total = perc_total.reindex(dat.index)
        # perc_total.index.equals(dat.index)
        dat = dat.astype(float) * perc_total.to_numpy().astype(float) * 100
    dat = pd.merge(diag, dat, right_index=True, left_index=True)
    return dat


df_perc_gated = format_data(counts_path='./Data/omiqData_formatted/230101_cluster_abundances_subs_corr.csv',
                            perc_path='./Data/230108_gated_populations_percentages.csv',
                            as_percent_total=False)
df = df_perc_gated.groupby(['Diagnosis']).agg('median')
df = df.transpose()
df.to_csv('./Data/230120_for_circos_no_correlation.csv')

# TODO: PLOT RADAR PLOTS OF GATED POPULATIONS
for cell in ['B', 'I', 'T']:
    df = df.filter(like=cell)
    df.columns = [re.sub('Cluster ', '', i) for i in df.columns]
    df = df.transpose().reset_index()

    cells = df['index']
    cells = [*cells, cells[0]]
    healthy = df['Healthy'] * 100
    mpn = df['MPN'] * 100
    mf = df['MF'] * 100
    healthy = [*healthy, healthy[0]]
    mpn = [*mpn, mpn[0]]
    mf = [*mf, mf[0]]

    r_max = (max(df.select_dtypes(include=[np.number]).max()) + 0.01) * 100
    fig = go.Figure(
        data=[
            go.Scatterpolar(r=healthy, theta=cells, name='Healthy', fill='toself', line=dict(color='#ed4135', width=1)),
            go.Scatterpolar(r=mpn, theta=cells, name='MPN', fill='toself', line=dict(color='#252160', width=1)),
            go.Scatterpolar(r=mf, theta=cells, name='MF', fill='toself', line=dict(color='#05cfc2', width=1))
        ],
        layout=go.Layout(
            title=go.layout.Title(text=f'{cell} abundances'),
            polar={'radialaxis': {'visible': True, 'range': [-3, r_max], 'nticks': 3, 'tickfont': {'size': 10}},
                   'angularaxis': {'tickfont': {'size': 10}}},
            showlegend=True
        )
    )

    # fig.update_polars(radialaxis_tickfont=dict(size=10))
    # pyo.plot(fig)
    fig.write_image(f"./Graphs/Cell repertoire/radar/{cell}.png", scale=5)


# TODO: PLOT GLOBAL RADAR BAR CHARTS
df_perc_total = format_data(counts_path='./Data/omiqData_formatted/230101_cluster_abundances_subs_corr.csv',
                            perc_path='./Data/230108_gated_populations_percentages.csv',
                            as_percent_total=True)

group = 'MF'
df = df_perc_total[df_perc_total['Diagnosis'] == group].copy()
df.drop(columns=['Diagnosis'], inplace=True)
df = df.agg('mean').reset_index()
df = df.rename(columns={0: 'mean abundance'})
df['cell'] = [re.sub(' Cluster.*?$', '', i) for i in df['index']]
df['index'] = [re.sub(' Cluster ', '', i) for i in df['index']]
cell_cols = {'B': '#0944a8', 'I': '#008000', 'T': '#a11313'}
df['cols'] = df['cell'].map(cell_cols)
df = df.sort_values(by='mean abundance')
df = df.iloc[len(df)-20:len(df), :].reset_index(drop=True)


plt.gcf().set_size_inches(5, 5)
sns.set_style('darkgrid')
max_val = 18
ax = plt.subplot(projection='polar')

for i in range(len(df)):
    ax.barh(i, df['mean abundance'][i]*2*np.pi/max_val,
            label=df['index'][i], color=cell_cols[df['cell'][i]])

ax.set_theta_zero_location('N')
ax.set_theta_direction(1)
ax.set_rlabel_position(0)
lines_th, labels_th = ax.set_thetagrids(range(0, 360, 40), labels=[str(i) + ' %' for i in range(0, 18, 2)], fontsize=8)
lines, labels = ax.set_rgrids([i-0.5 for i in range(len(df))],
                              labels=df['index'], fontfamily='Courier', fontsize=8)
for i in lines:
    i.set_alpha(0.2)
plt.savefig(f"./Graphs/Abundances in groups/global_{group}.png", dpi=300)

