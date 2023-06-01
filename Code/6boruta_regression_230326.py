from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import *
import matplotlib.pyplot as plt
from collections import Counter
from itertools import compress
from datetime import datetime
import statsmodels.api as sm
from boruta import BorutaPy
import seaborn as sns
import pandas as pd
import numpy as np
import re
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)


def abundances_with_HAI():
    """
    Add UNSTIM cluster abundances and metadata to formatted HAI titer data
    :return: Formatted abundances with HAI and Diagnosis column
    """
    serology = pd.read_csv('./Data/230326_formatted_HAI_means.csv', dtype={'Sample ID': str})
    serology = serology.drop(columns={'Diagnosis', 'vaccine_year'})
    serology = serology.pivot(columns='group', index='Sample ID')
    serology.columns = [re.sub('^', 'HAI_', i[1]) for i in serology.columns]

    df = pd.read_csv("./Data/omiqData_formatted/230324_cluster_abundances.csv")
    df['sample_name'] = [re.sub('.*?-(.*?) (.*?)_Hobbs.*?$', '\\1_\\2', i) for i in df['sample_name']]
    df[['Sample ID', 'Stimulation']] = df['sample_name'].str.split(pat='_', n=1, expand=True)
    df = df[df["Stimulation"] == 'Unstim'].drop(columns={'Stimulation'})
    df = df.drop(columns={'sample_name'})
    df = pd.merge(df, serology, left_on="Sample ID", right_index=True)
    return df


def pred_vs_og_plot():
    # plot predicted vals vs original values
    plt.scatter(y, y_pred)
    plt.axline([0, 0], [1, 1], c="#aaf0d1")
    plt.xlabel('original values')
    plt.ylabel('predicted values')
    plt.title(group)
    plt.text(x=y.max()*0.65, y=y_pred.max()*0.15, s=f'best score = {best_score}\nrmse = {rmse}\nr2 = {r2}')
    plt.savefig(f'./Graphs/boruta/{today}_{group}_pred_vs_original.jpeg', dpi=200)
    plt.close()


def plot_waterfall(df):
    col_dict = {'B': '#0944a8', 'I': '#008000', 'T': '#a11313'}
    df['Cell'] = [re.sub('^(.).*?$', '\\1', i) for i in df['Cluster']]
    df = df.sort_values(by='Coef', ascending=False)
    if len(df) < 10:
        plt.figure(figsize=(4, 3))
    else:
        plt.figure(figsize=(6, 4))
    sns.barplot(x="Coef", y="Cell type with cluster", data=df, dodge=False, palette=col_dict, hue='Cell')
    plt.legend(handles=[], labels=[], frameon=False)
    plt.xlabel('change in LogOdds per unit SD')
    plt.ylabel('')
    plt.suptitle(f'{group}')
    plt.tight_layout()
    plt.savefig(f"./Graphs/boruta/{today}_{group}_{suffix}.jpg", dpi=300)
    plt.close()


today = datetime.now().strftime('%y%m%d')
annotation = pd.read_excel('./Data/Cluster_phenotypes.xlsx')
dat = abundances_with_HAI()

group = 'HAI_MPN'
to_model = dat.dropna(subset=[group])
X = to_model.drop(columns={'Sample ID', 'HAI_MPN', 'HAI_Healthy'})
y = to_model[group]

forest = RandomForestRegressor(random_state=42)
perc = 90
n_est = 500
max_iter = 200
boruta = BorutaPy(verbose=2, estimator=forest, n_estimators=n_est, perc=perc, max_iter=max_iter)
boruta.fit(np.array(X), np.array(y))  # train Boruta
X_filtered = boruta.transform(np.array(X))  # filter X to contain selected features

param_grid = {
    'n_estimators': [100, 200, 300, 400, 500],
    'max_features': ['sqrt', 'log2', 0.33, 10],
    'max_depth': [3, 4, 5, 6, 7, 8],
}

CV_rf = GridSearchCV(estimator=forest, param_grid=param_grid, cv=2)
CV_rf.fit(X_filtered, y)
param = CV_rf.best_params_
# fit an RF with filtered X
y_pred = CV_rf.predict(X_filtered)  # predict y based on the RF model
# metrics
mse = ((y_pred - y) ** 2).mean()
rmse = round(np.sqrt(mse), 3)
r2 = round(r2_score(y, y_pred), 3)
best_score = round(CV_rf.best_score_, 3)
print(f'Scores:\nbest score = {best_score}\nrmse = {rmse}\nr2 = {r2}')
pred_vs_og_plot()

scaler = StandardScaler()
x_scaled = scaler.fit_transform(X_filtered)
y_scaled = scaler.fit_transform(np.array(y).reshape(-1, 1))
model = sm.OLS(y_scaled, x_scaled)
results = model.fit()
print(results.summary())
selected_features = list(compress(X.columns, boruta.support_))  # get selected features
coef_df = pd.DataFrame({'Cluster': selected_features, 'Coef': results.params})

coef_df = pd.merge(left=coef_df, right=annotation[['Cluster', 'Count', 'Cell type']], on='Cluster')
coef_df['Cell type with cluster'] = coef_df['Cell type'] + ' (' + coef_df['Cluster'] + ')'
coef_df = coef_df[coef_df['Count'] > 1000]
coef_df['Cell'] = [re.sub('^(.).*?$', '\\1', i) for i in coef_df['Cluster']]

plot_waterfall(coef_df)

params = pd.DataFrame.from_dict({'name': group,
                                 'n': list(Counter(y).values()),
                                 'RF_max_depth': param['max_depth'],
                                 'RF_max_ft': param['max_features'],
                                 'RF_B_n_estimators': param['n_estimators'],
                                 'B_percentile': perc,
                                 'B_max_iter': max_iter,
                                 'Best model score': CV_rf.best_score_,
                                 'n_features': len(selected_features)}, orient='index')

writer = pd.ExcelWriter(f'./Data/results/Boruta/{today}_{group}.xlsx')
coef_df.to_excel(writer, sheet_name='Logistic Regression coefs', index=False)
params.to_excel(writer, sheet_name='Boruta params')
writer.close()
