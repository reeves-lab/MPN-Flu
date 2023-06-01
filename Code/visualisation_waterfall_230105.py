import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
import glob


def plot_waterfall(df, dir):
    if len(df) < 10:
        plt.figure(figsize=(4, 3))
    else:
        plt.figure(figsize=(6, 4))
    sns.barplot(x="coef", y="cluster", data=df, hue='cell', dodge=False, palette=col_dict)
    plt.legend(handles=[], labels=[], frameon=False)
    plt.xlabel('change in LogOdds per unit SD')
    plt.ylabel('')
    plt.suptitle(f'{name}')
    plt.tight_layout()
    plt.savefig(f"{dir}/{name}.jpg", dpi=300)
    plt.close()


col_dict = {'B': '#0944a8', 'I': '#008000', 'T': '#a11313'}
# TODO 1: Boruta pretty waterfall plots
files = glob.glob('./Data/results/Boruta/*')
for path in files:
    # path = files[0]
    file = pd.read_excel(path)
    name = re.sub('^.*?230104_(.*?)_coef.*?.xlsx', '\\1', path)
    file['cell'] = [re.sub('^(.) Clust.*?$', '\\1', i) for i in file['cluster']]
    file = file.sort_values(by='coef', ascending=False)

    plot_waterfall(file, './Graphs/Boruta')

# TODO 2: Elastic Net waterfall plots (effect of stimulation)
file = pd.read_csv('./Data/results/Elastic Net/effectofStim_validation_coefficients.csv')
file.columns = [i.lower() for i in file.columns]
file['comparison'] = [re.sub('^.*?Unstim_', 'Unstim vs Stim ', i) for i in file['classes']]
file['coef'] *= -1  # only for effect of stim
# in en the comparison was done with stim vs unstim (rather than the other way round)
file['cell'] = [re.sub('^(.) Clust.*?$', '\\1', i) for i in file['cluster']]

for name in file['comparison'].unique():
    per = file[file['comparison'] == name].copy()
    per = per.sort_values(by='coef', ascending=False)
    plot_waterfall(per, './Graphs/Elastic Net')

# TODO 3: Elastic Net waterfall plots (stim and unstim comparisons)
stimulation = 'Stim'
file = pd.read_csv(f'./Data/results/Elastic Net/{stimulation}_validation_coefficients.csv')
file.columns = [i.lower() for i in file.columns]
file['cell'] = [re.sub('^(.) Clust.*?$', '\\1', i) for i in file['cluster']]

for name in file['classes'].unique():
    per = file[file['classes'] == name].copy()
    per = per.sort_values(by='coef', ascending=False)
    name = name + '_' + stimulation
    plot_waterfall(per, './Graphs/Elastic Net')

