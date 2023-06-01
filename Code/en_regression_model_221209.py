from sklearn.model_selection import GridSearchCV
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import ElasticNet
import matplotlib.pyplot as plt
from datetime import datetime
import pandas as pd
import numpy as np


class FeatureSelection:
    def __init__(self, **kw):
        self.id = kw.get('id')
        self.a_th = kw.get('alpha_th')
        self.l1_th = kw.get('l1_th')
        self.run = kw.get('run')
        self.reg_name = kw.get('reg_name')
        self.today = datetime.today().strftime("%y%m%d")
        self.best_estimator = 0
        self.best_model_score = 0
        self.output = pd.DataFrame()
        self.model_stats = pd.DataFrame()
        self.summary = pd.DataFrame()
        self.score_comparison = None
        self.ft_n_list = []

    def run_en(self, to_model, **kw):
        y = to_model[self.reg_name]  # subset response column and normalise
        # y_scaler = StandardScaler()
        # y_scaler.fit(y.values.reshape(-1, 1))
        # y = y_scaler.transform(y.values.reshape(-1, 1))

        if self.run == "test":  # subset features and normalise
            # feature_df = to_model.drop(columns={'Sample ID', 'HAI_Healthy', 'HAI_MPN'})
            feature_df = to_model.drop(columns={self.reg_name})
            feature_df = feature_df.dropna(axis=1)  # remove columns with nans
        else:
            feature_df = to_model.drop(self.reg_name, axis=1)
        features = feature_df.columns
        x_scaler = StandardScaler()
        x = x_scaler.fit_transform(feature_df)


        log_en = ElasticNet(max_iter=10000, warm_start=True)
        param_grid = [
            {'alpha': np.linspace(self.a_th[0], self.a_th[1], num=100),
             'l1_ratio': np.linspace(self.l1_th[0], self.l1_th[1], num=21)}
        ]
        clf = GridSearchCV(log_en, param_grid, cv=3, scoring="neg_mean_absolute_error", n_jobs=-1, verbose=False)
        clf.fit(x, y)

        # pred_y = y_scaler.inverse_transform(clf.predict(x).reshape(-1, 1))  # generate predicted values
        pred_y = clf.predict(x)  # generate predicted values
        plt.scatter(to_model[self.reg_name], pred_y)  # plot predicted vals vs original values
        plt.axline([0, 0], [1, 1], c="r")
        plt.xlabel("Original values")
        plt.ylabel("Predicted values")
        plt.savefig(f"../Graphs/Elastic Net/{self.today}_{self.reg_name}_pred_vs_og_{self.id}.png", dpi=200)
        plt.close()

        # model statistics
        self.model_stats = pd.DataFrame(clf.cv_results_)
        self.model_stats.sort_values(by="rank_test_score", inplace=True, ignore_index=True)

        # performance statistics
        self.best_model_score = clf.best_score_
        self.best_estimator = clf.best_estimator_
        coef = self.best_estimator.coef_
        feature_num = len(np.nonzero(coef)[0])

        # make feature output
        self.output = pd.DataFrame([features[np.nonzero(coef)], coef[np.nonzero(coef)]]).transpose()
        self.output.columns = ["Cluster", "Coef"]
        self.output.Coef = self.output.Coef.astype(float)
        self.output = self.output.sort_values(by='Coef', ascending=False).reset_index(drop=True)
        self.output["Regression"] = self.reg_name

        # make a summary dataframe
        self.summary = pd.DataFrame.from_dict({
            'Name': self.reg_name,
            'Feature number': feature_num,
            'Score': self.best_model_score,
            'Alpha': self.best_estimator.alpha,
            'l1_ratio': self.best_estimator.l1_ratio
        }, orient='index')

        if self.run != "test":
            if self.best_model_score > kw.get('test_score'):
                self.summary.loc["score_comparison"] = "better"
            else:
                self.summary.loc["score_comparison"] = "same_or_worse"

        writer = pd.ExcelWriter(f"../Data/results/Elastic Net/{self.today}_{self.run}_{self.reg_name}_{self.id}.xlsx")
        self.output.to_excel(writer, sheet_name='Cluster output', index=False)
        self.summary.to_excel(writer, sheet_name='Model summary')
        self.model_stats.to_excel(writer, sheet_name='Model stats', index=False)
        writer.close()

