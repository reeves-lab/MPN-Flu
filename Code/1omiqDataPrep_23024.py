from datetime import datetime
import pandas as pd
import re


def format_omiq_data(df):
    df["sample_name"] = [re.sub('[0-9]{2}-([0-9]*?) (.*?)_.*?$', '\\1_\\2', i) for i in df['file']]
    df.drop(columns=["file"], inplace=True)
    df = pd.melt(df, id_vars=['sample_name'])
    df["cluster"] = df["variable"].replace(to_replace="cl-([0-9]{2}).*?$", value=f"{cell}\\1", regex=True)
    df["marker"] = df["variable"].replace(to_replace="^.*?:: (.*?)\|.*?$", value="\\1 median", regex=True)
    df = df[df['sample_name'].str.contains('Unstim')]

    abundances = df[df["variable"].str.contains("percent_parent")][["sample_name", "cluster", "value"]]
    abundances.rename(columns={"value": "abundance"}, inplace=True)
    medians = df[df["variable"].str.contains("median")][["sample_name", "cluster", "marker", "value"]]

    final = pd.merge(left=medians, right=abundances, on=["sample_name", "cluster"])
    final = final.pivot(index=['sample_name', 'cluster', 'abundance'], columns='marker', values='value').reset_index()
    df_li.append(final)


today = datetime.now().strftime("%y%m%d")

df_li = []
for cell in ['B', 'I', 'T']:
    # cell = 'B'
    data = pd.read_csv(f'Data/omiqData/230324_markers_perc_{cell}.csv', index_col=None, header=0, dtype=object)
    format_omiq_data(df=data)


dat = pd.concat(df_li)  # concatenate the frames together (B, I, T)
# making separate files for cluster counts/abundances (pivot table with shape: files x clusters)
# and for cluster phenotype, which is still a melted, long format df, but without cell count per cluster
cluster_counts = dat.iloc[:, 0:3].pivot(index="sample_name", columns='cluster', values='abundance')
cluster_phenotypes = dat.drop(["abundance"], axis=1)
cluster_phenotypes.to_csv(f"Data/omiqData_formatted/{today}_cluster_phenotypes.csv", index=False)
cluster_counts.to_csv(f"Data/omiqData_formatted/{today}_cluster_abundances.csv")
