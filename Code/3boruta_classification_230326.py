from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from collections import Counter
import matplotlib.pyplot as plt
from itertools import compress
from datetime import datetime
from sklearn.metrics import *
from boruta import BorutaPy
import seaborn as sns
import pandas as pd
import numpy as np
import re
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)


def format_data():
    clinical = pd.read_excel("./Data/221027_metadata.xlsx", sheet_name="Combined", dtype={'Sample ID': str})
    clinical['Diagnosis'] = clinical['Diagnosis'].fillna('Healthy')
    clinical['Diagnosis'] = [re.sub('MPN', 'PV_ET', j) for j in clinical['Diagnosis']]
    df = pd.read_csv("./Data/omiqData_formatted/230324_cluster_abundances.csv")
    df[['Sample ID', 'Stimulation']] = df['sample_name'].str.split(pat='_', n=1, expand=True)
    df = df[df["Stimulation"] == 'Unstim'].drop(columns={'Stimulation', 'sample_name'})
    df = pd.merge(df, clinical[['Diagnosis', 'Sample ID']], on="Sample ID")
    comparisons_ = [["Healthy", "PV_ET"], ["Healthy", "MF"]]
    return df, comparisons_


def plot_waterfall(df):
    if len(df) < 10:
        plt.figure(figsize=(4, 3))
    else:
        plt.figure(figsize=(6, 4))
    col_dict = {'B': '#0944a8', 'I': '#008000', 'T': '#a11313'}
    df = df.sort_values(by='coef', ascending=False)
    df['cell'] = [re.sub('^(.).*?$', '\\1', i) for i in df['cluster']]
    df['Cell_type_w_cluster'] = df['Cell type'] + ' (' + df['cluster'] + ')'
    sns.barplot(x="coef", y="Cell_type_w_cluster", data=df, hue='cell', dodge=False, palette=col_dict)
    plt.legend(handles=[], labels=[], frameon=False)
    plt.xlabel('change in LogOdds per unit SD')
    plt.ylabel('')
    plt.suptitle(f'{name}')
    plt.tight_layout()
    plt.savefig(f"Graphs/Boruta/{name}.jpg", dpi=300)
    plt.close()


today = datetime.now().date().strftime("%y%m%d")
startTime = datetime.now()
data, comparisons = format_data()
annotation = pd.read_excel('./Data/Cluster_phenotypes.xlsx')

for comparison in comparisons:
    # comparison = comparisons[0]
    name = "_".join(comparison)
    to_model = data.loc[data['Diagnosis'].isin(comparison)]
    X = to_model.drop(columns={'Sample ID', 'Diagnosis'})
    y = to_model['Diagnosis'].values.ravel()
    x_train, x_test, y_train, y_test = train_test_split(np.array(X), y, test_size=.3, random_state=42)

    forest = RandomForestClassifier(random_state=42)
    perc = 90
    max_iter = 100
    boruta = BorutaPy(verbose=2, estimator=forest, perc=perc, max_iter=max_iter)
    boruta.fit(x_train, y_train)
    selected_features = list(compress(X.columns, boruta.support_))  # get selected features
    x_train2 = boruta.transform(x_train)
    param_grid = {
        'n_estimators': [100, 200, 300, 400, 500],
        'max_features': ['sqrt', 'log2', 0.33, 10],
        'max_depth': [4, 5, 6, 7, 8],
        'criterion': ['gini', 'entropy']
    }
    CV_rf = GridSearchCV(estimator=forest, param_grid=param_grid, cv=3)
    CV_rf.fit(x_train2, y_train)
    param = CV_rf.best_params_
    # get accuracy for training data
    y_pred = CV_rf.predict(x_train2)
    accuracy_FS_train = accuracy_score(y_train, y_pred)
    # get accuracy for test data
    x_test2 = boruta.transform(x_test)
    y_pred = CV_rf.predict(x_test2)
    accuracy_FS_test = accuracy_score(y_test, y_pred)
    print(f'Accuracy w/ feature selection, for train: {accuracy_FS_train}; for test: {accuracy_FS_test}')
    print(f'Best model score: {CV_rf.best_score_}')

    # run logistic regression on test data to get coefficients of selected features
    scaler = StandardScaler()
    x_train2_scaled = scaler.fit_transform(x_train2)
    x_test2_scaled = scaler.fit_transform(x_test2)

    clf = LogisticRegression(random_state=0).fit(x_train2_scaled, y_train)
    accuracy_LR = clf.score(x_test2_scaled, y_test)
    print(accuracy_LR)

    class_count_train = list(Counter(y_train).values())
    class_count_test = list(Counter(y_test).values())

    params = pd.DataFrame.from_dict({'name': name, 'log_reg_score': accuracy_LR,
                                     'n_in_train_classes': class_count_train,
                                     'n_in_test_classes': class_count_test,
                                     'RF_max_depth': param['max_depth'],
                                     'RF_max_ft': param['max_features'],
                                     'RF_B_n_estimators': param['n_estimators'],
                                     'B_percentile': perc,
                                     'B_max_iter': max_iter,
                                     'accuracy_test': accuracy_FS_test,
                                     'Best model score': CV_rf.best_score_,
                                     'n_features': len(selected_features)}, orient='index')
    coef_df = pd.DataFrame({'cluster': selected_features, 'coef': clf.coef_[0]})
    coef_df = pd.merge(coef_df, annotation[['Cluster', 'Cell type']], left_on='cluster', right_on='Cluster')
    coef_df = coef_df[coef_df['Cell type'] != '-']  # - marks clusters with fewer than 1000 cells
    plot_waterfall(df=coef_df)

    writer = pd.ExcelWriter(f'./Data/results/Boruta/{today}_{name}.xlsx')
    coef_df.to_excel(writer, sheet_name='Logistic Regression coefs', index=False)
    params.to_excel(writer, sheet_name='Boruta params')
    writer.close()
    print('DONE!')
