import pandas as pd
import numpy as np
from datetime import datetime
import re
import glob
import json

today = datetime.now().strftime("%y%m%d")
files = glob.glob("./Data/omiqData/*markers*.csv")
file_list = []

for file in files:
    # file = files[0]
    cell = re.sub("^.*?markers_(.).csv", "\\1", file)
    file = pd.read_csv(f"{file}")
    file = file.rename(columns={"Unnamed: 0": "cluster"})
    file["cluster"] = [re.sub("concat.*?cl-", f"{cell} Cluster ", i) for i in file["cluster"]]
    file.columns = [re.sub("^.*? :: ", "", i) for i in file.columns]
    file.columns = [re.sub("^.*?_", "", i) for i in file.columns]
    file = file.reindex(sorted(file.columns), axis=1)
    file = file.sort_values(by="cluster")
    file_list.append(file)

data = pd.concat(file_list)
data = data.reset_index(drop=True)
data = data.replace({'0': np.nan, 0: np.nan})

for col in data.columns[0:41]:
    data[col] = data[col].rank(pct=True) * 100

data = data.fillna(0)

annot = pd.concat([data.copy()] * 2, axis=1)
new_list = data.columns.tolist()
new_list.extend([f"{col}_annot" for col in data.columns])
annot.columns = new_list

for col in annot.columns:
    # cd_ = "CCL4"
    if "annot" not in col and "cluster" not in col:
        annot.loc[annot[col] <= 25, f"{col}_annot"] = 0
        annot.loc[(annot[col] > 25) & (annot[col] <= 50), f"{col}_annot"] = f"{col} lo"
        annot.loc[(annot[col] > 50) & (annot[col] <= 75), f"{col}_annot"] = f"{col} mid"
        annot.loc[(annot[col] > 75), f"{col}_annot"] = f"{col} hi"

annot_data = annot[[col for col in annot.columns if "annot" in col]]
annot_data.columns = [re.sub("_annot", "", col) for col in annot_data.columns]
correct_order = ["cluster"]
correct_order.extend(annot_data.columns[0:41].tolist())

annot_data = annot_data[correct_order]

li = []
for index, row in annot_data.iterrows():
    row_ = row.tolist()
    row_ = [x for x in row_ if x != 0 and 'lo' not in x]
    li.append(row_)

items = {}
for line in li:
    key, value = line[0], ', '.join(line[1:])
    items[key] = value

with open(f'./Data/{today}_cluster_annotation.json', 'a') as f:
    json.dump(items, f)


for i in li:
    with open("test.txt", "a") as f:
        row = ", ".join(e for e in i)
        f.writelines(f"{row}\n")

# annot_data.to_csv("cluster_annotation.csv", index=False)
