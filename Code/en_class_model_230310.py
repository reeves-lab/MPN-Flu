from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import StandardScaler
import matplotlib.pyplot as plt
from datetime import datetime
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib
import re
matplotlib.use('Agg')


class FeatureSelection:
    def __init__(self, **kw):
        self.run = kw.get('run')  # test or validation
        self.classification = kw.get('classification')  # which classification is run
        self.classify_on_name = ''
        self.best_estimator = 0
        self.best_model_score = 0
        self.output = pd.DataFrame()
        self.model_stats = pd.DataFrame()
        self.summary = pd.DataFrame()
        self.score_comparison = None
        self.ft_n_list = []
        self.today = datetime.today().strftime('%y%m%d')

    def run_en(self, to_model, **kw):

        features = to_model.columns[to_model.columns != self.classification]
        y = to_model[self.classification]
        x = to_model[to_model.columns[to_model.columns != self.classification]]
        x_scaler = StandardScaler()
        x_scaler.fit(x)
        x = x_scaler.transform(x)

        log_en = LogisticRegression(penalty="elasticnet", max_iter=10000, warm_start=True, solver="saga")
        param_grid = [{"C": np.linspace(0.01, 1, num=100), 'l1_ratio': np.linspace(0.2, 1, num=21)}]
        clf = GridSearchCV(
            log_en,
            param_grid,
            cv=4,
            scoring="balanced_accuracy",
            n_jobs=-1,
            verbose=False
        )
        clf.fit(x, y)
        self.best_model_score = clf.best_score_
        self.model_stats = pd.DataFrame(clf.cv_results_)
        self.model_stats.sort_values(by="rank_test_score", inplace=True, ignore_index=True)
        self.best_estimator = clf.best_estimator_
        print(self.best_estimator)
        print(self.best_estimator.classes_)

        self.classify_on_name = " vs ".join(self.best_estimator.classes_)
        regression_coef = self.best_estimator.coef_[0]
        feature_num = len(np.nonzero(regression_coef)[0])
        nonzero_coef_arr = np.nonzero(regression_coef)[0]
        if feature_num:  # execute the next steps if >0 fts are selected in validation
            self.ft_n_list.append(feature_num)
            self.output = pd.DataFrame(
                [features[nonzero_coef_arr], regression_coef[nonzero_coef_arr]]).transpose()
            self.output.columns = ["Cluster", "Coef"]
            self.output = self.output.sort_values(by="Coef", ascending=False).reset_index(
                drop=True)
            self.output["Classification"] = self.classification
            self.output["Classes"] = self.classify_on_name

            # make a summary dataframe
            self.summary = pd.DataFrame.from_dict({
                'Name': self.classification,
                'Feature number': feature_num,
                'Score': self.best_model_score,
                'C': self.best_estimator.C,
                'l1_ratio': self.best_estimator.l1_ratio
            }, orient='index')

            if self.run != "test":
                if self.best_model_score > kw.get('test_score'):
                    self.summary.loc["score_comparison"] = "better"
                else:
                    self.summary.loc["score_comparison"] = "same_or_worse"

            writer = pd.ExcelWriter(f"../Data/results/Elastic Net/{self.today}_{self.run}_{self.classify_on_name}.xlsx")
            self.output.to_excel(writer, sheet_name='Cluster output', index=False)
            self.summary.to_excel(writer, sheet_name='Model summary')
            self.model_stats.to_excel(writer, sheet_name='Model stats', index=False)
            writer.close()
