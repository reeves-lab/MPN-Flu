from distinctipy import distinctipy
import matplotlib.pyplot as plt
from datetime import datetime
import seaborn as sns
import pandas as pd
import umap
import re


def format_data(omiq_path, meta_path):
    omiq_df = pd.read_csv(omiq_path)
    omiq_df["sample_name"] = [re.sub('[0-9]{2}-([0-9]*?) (.*?)_.*?$', '\\1_\\2', i) for i in omiq_df['OmiqFileIndex']]
    omiq_df[['Sample ID', 'Stimulation']] = omiq_df['sample_name'].str.split(pat='_', n=1, expand=True)
    omiq_df.drop(columns='OmiqFileIndex', inplace=True)
    omiq_df['OmiqFilter'] = [re.sub('cl-', '', i).zfill(2) for i in omiq_df['OmiqFilter'].astype(str)]
    omiq_df.rename(columns={'OmiqFilter': 'cluster'}, inplace=True)
    omiq_df.columns = [re.sub('^.*?:: ', '', i) for i in omiq_df.columns]
    # get metadata
    metadata = pd.read_excel(meta_path, sheet_name='Combined', dtype={'Sample ID': str})
    metadata['Diagnosis'] = metadata['Diagnosis'].fillna('Healthy')
    # merge omiq frame and metadata
    merged = pd.merge(omiq_df, metadata[['Sample ID', 'Diagnosis']], on='Sample ID', how='left')
    return merged


def balance_and_umap(df, n_cells, features, n_neigh, min_dist):
    balanced = pd.DataFrame()
    for sample in df['sample_name'].unique():
        single_sample = df.loc[df['sample_name'] == sample]
        if single_sample['sample_name'].value_counts().values[0] < n_cells:
            balanced = pd.concat([balanced, single_sample])
        else:
            balanced = pd.concat([balanced, single_sample.sample(n_cells)])

    # umap dataframe and parameters
    balanced = balanced.reset_index(drop=True)
    if len(features) != 0:
        for_umap = balanced[features].values
    else:
        cols_to_drop = ['cluster', 'gate', 'Sample ID', 'Stimulation', 'Diagnosis', 'sample_name']
        for_umap = balanced.drop(columns=cols_to_drop).values

    reducer = umap.UMAP(n_neighbors=n_neigh, n_epochs=100, min_dist=min_dist)
    start_time = datetime.now()
    embedding = reducer.fit_transform(for_umap)
    print(f'Runtime: {datetime.now() - start_time}\nShape of embedding: {embedding.shape}')

    embedding = pd.merge(pd.DataFrame(embedding), balanced, right_index=True, left_index=True)
    embedding.rename(columns={0: 'u1', 1: 'u2'}, inplace=True)
    return embedding


def plot_global(u_df, plot_markers, hue_, **kw):
    sns.set_style('white')
    plt.figure(figsize=(4, 4.5))

    if hue_ == 'minor':
        clust_n = len(u_df[hue_].unique())
        colors = distinctipy.get_colors(clust_n)
        ncol_ = 2

    else:
        colors = kw.get('colors')
        ncol_ = 1

    u_df[hue_] = pd.Categorical(u_df[hue_], ordered=True, categories=kw.get('sorted_hue'))
    ax = sns.scatterplot(data=u_df, x="u1", y="u2", s=2, hue=u_df[hue_], linewidth=0, palette=colors)
    lgd = ax.legend(loc='center left', ncol=ncol_, frameon=False, bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.savefig(f"Graphs/umaps/{today}_{cell}_clusters_{hue_}.jpg", bbox_extra_artists=(lgd,), bbox_inches='tight',
                dpi=500)
    plt.close()

    if plot_markers:
        for marker in u.columns[2:43]:
            plt.figure(figsize=(3, 3))
            ax = sns.scatterplot(data=u, x='u1', y='u2', s=2, hue=marker, linewidth=0, palette='viridis')
            norm = plt.Normalize(-1, 1)
            sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
            sm.set_array([])
            ax.get_legend().remove()
            cbar = ax.figure.colorbar(sm, ax=ax)
            cbar.set_ticks([-1, -0.5, 0, 0.5, 1])
            plt.savefig(f"Graphs/umaps/{cell} markers/{marker}.jpg", dpi=300)
            plt.close()


def plot_umap(u_df, plot_markers, hue_, **kw):
    sns.set_style('white')
    plt.figure(figsize=(4, 4.5))
    clust_n = len(u_df[hue_].unique())

    if hue_ == 'cluster':
        u_df.sort_values(hue_, inplace=True)
        ncol_ = 2
        colors = distinctipy.get_colors(clust_n)
    else:
        u_df = u_df[u_df[hue_].isin(kw.get('sorted_hue'))]
        u_df[hue_] = pd.Categorical(u_df[hue_], ordered=True, categories=kw.get('sorted_hue'))
        ncol_ = 1
        colors = sns.color_palette("bright", len(kw.get('sorted_hue')))

    ax = sns.scatterplot(data=u_df, x="u1", y="u2", s=2, hue=u_df[hue_], linewidth=0, palette=colors)
    lgd = ax.legend(loc='center left', ncol=ncol_, frameon=False, bbox_to_anchor=(1, 0.5), fontsize=13)
    plt.savefig(f"Graphs/umaps/{today}_{cell}_{hue_}.jpg", bbox_extra_artists=(lgd,), bbox_inches='tight',
                dpi=300)
    plt.close()

    if plot_markers:
        for marker in u.columns[2:43]:
            plt.figure(figsize=(4, 3))
            ax = sns.scatterplot(data=u, x='u1', y='u2', s=2, hue=marker, linewidth=0, palette='viridis')
            norm = plt.Normalize(-1, 1)
            sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
            sm.set_array([])
            ax.get_legend().remove()
            cbar = ax.figure.colorbar(sm, ax=ax)
            cbar.set_ticks([-1, -0.5, 0, 0.5, 1])
            plt.savefig(f"Graphs/umaps/{cell} markers/{marker}.jpg", dpi=300)
            plt.close()


def plot_contours(u_df):
    u_df.Diagnosis = pd.Categorical(u_df.Diagnosis, categories=['Healthy', 'MPN', 'MF'])
    u_df = u_df.sort_values('Diagnosis')
    cmap = ['#ed4135', '#252160', '#05cfc2']
    sns.set_style('white')
    sns.displot(data=u_df, x="u1", y="u2", hue=u_df['Diagnosis'], kind="kde", palette=cmap, bw_adjust=.7, alpha=0.8)
    plt.savefig(f"Graphs/umaps/{today}_{cell}_contours.jpg", dpi=300)
    plt.close()


# TODO 1: Major gate-specific UMAPs (B, I, T)
cell = 'T'
today = datetime.now().strftime('%y%m%d')
data = format_data(omiq_path=f"Data/omiqData/for UMAP {cell} clusters.csv", meta_path="Data/221027_metadata.xlsx")
data_cells = format_data(omiq_path=f"Data/omiqData/for UMAP {cell} gates.csv", meta_path="Data/221027_metadata.xlsx")
li = list(data.columns)
li.remove('cluster')
data = pd.merge(data, data_cells, on=li)
data = data[data['Stimulation'] == 'Unstim']
data = data.rename(columns={'cluster_x': 'cluster', 'cluster_y': 'gate'})
u = balance_and_umap(df=data, features=[], n_cells=5000, n_neigh=50, min_dist=0.3)

order = {
    'B': ['Naive B', 'Plasmablast (38hi27hi)', 'Memory B', 'other B'],
    'I': ['Lin neg', 'NK cell (56hi16neg)', 'NK cells (56hi16hi)', 'NK cells (56hi16lo)',
          'Classical monocytes', 'Nonclassical monocytes', 'Intermediate monocytes', 'Dendritic cells'],
    'T': ['NK T cells', 'DN T cells',
          'Naive helper T', 'Effector helper T', 'CM helper T', 'EM helper T',
          'Naive cytotoxic T', 'Effector cytotoxic T', 'CM cytotoxic T', 'EM cytotoxic T']
}

plot_umap(u_df=u, plot_markers=False, hue_='gate', sorted_hue=order[cell])
plot_umap(u_df=u, plot_markers=True, hue_='cluster')
plot_contours(u_df=u)


# TODO 2: Global UMAP
cell = 'global'
today = datetime.now().strftime('%y%m%d')
order = {
    'B': ['Naive B', 'Plasmablast (38hi27hi)', 'Memory B', 'other B'],
    'I': ['Lin neg', 'NK cell (56hi16neg)', 'NK cells (56hi16hi)', 'NK cells (56hi16lo)',
          'Classical monocytes', 'Nonclassical monocytes', 'Intermediate monocytes', 'Dendritic cells'],
    'T': ['NK T cells', 'DN T cells',
          'Naive helper T', 'Effector helper T', 'CM helper T', 'EM helper T',
          'Naive cytotoxic T', 'Effector cytotoxic T', 'CM cytotoxic T', 'EM cytotoxic T']
}
data = format_data(omiq_path="Data/omiqData/for UMAP global major gates.csv", meta_path="Data/221027_metadata.xlsx")
data2 = format_data(omiq_path="Data/omiqData/for UMAP global minor gates.csv", meta_path="Data/221027_metadata.xlsx")
li = list(data.columns)
li.remove('cluster')
data = pd.merge(data, data2, on=li)
data = data[data['Stimulation'] == 'Unstim']
data = data.rename(columns={'cluster_x': 'major', 'cluster_y': 'minor'})
u = balance_and_umap(df=data, features=['CD16', 'CD11c', 'CD14', 'CD19', 'CD56', 'CD3', 'CD11b', 'CD4', 'CD8'],
                     n_cells=5000, n_neigh=50, min_dist=0.3)

u['major'] = [re.sub(' [(]CD14neg[)]', '', i) for i in u['major']]
sorted_ = order['B'] + order['I'] + order['T']
plot_global(u_df=u, plot_markers=False, hue_='major',
            colors=['#0944a8', '#008000', '#a11313'], sorted_hue=['B cells', 'Innate cells', 'T cells'])
plot_global(u_df=u, plot_markers=False, hue_='minor', sorted_hue=sorted_)
plot_contours(u_df=u)

# u['group_clusters'] = u['cluster']
# u['hai_clusters'] = u['cluster']
# u.loc[~u['cluster'].isin(group_clusters[cell]), 'group_clusters'] = 'ns'
# u.loc[~u['cluster'].isin(hai_clusters[cell]), 'hai_clusters'] = 'ns'

# plot_umap(u_df=u,
#           plot_markers=True,
#           col_by_cluster=True,
#           hue_col='group_clusters',
#           only_relevant_clusters=True,
#           suffix='group_clusters')

