import pandas as pd
from statsmodels.formula.api import ols
import statsmodels.stats.multicomp as mc
import statsmodels.api as sm


def format_anova_tbl(aov):
    """
    Calculates eta and omega squared and add them as columns to anova table
    :param aov: anova table
    :return: formatted anova table
    """
    aov['mean_sq'] = aov[:]['sum_sq'] / aov[:]['df']
    aov['eta_sq'] = aov[:-1]['sum_sq'] / sum(aov['sum_sq'])
    aov['omega_sq'] = (aov[:-1]['sum_sq'] - (aov[:-1]['df'] * aov['mean_sq'][-1])) / (
            sum(aov['sum_sq']) + aov['mean_sq'][-1])

    columns = ['sum_sq', 'df', 'mean_sq', 'F', 'PR(>F)', 'eta_sq', 'omega_sq']
    aov = aov[columns]
    return aov


meta = pd.read_excel("./Data/221027_metadata.xlsx", sheet_name="Combined", dtype={'Sample ID': str})

model = ols('Age ~ C(Group)', data=meta).fit()
aov_table = sm.stats.anova_lm(model, typ=2)
aov_table = format_anova_tbl(aov_table)


gated_perc = pd.read_csv('./Data/230108_gated_populations_percentages.csv')
gated_perc.rename(columns={'normalised to sample': 'norm_value'}, inplace=True)
per_cell = gated_perc[gated_perc['variable'] == 'Innate cells']
model = ols('norm_value ~ C(Diagnosis)', data=per_cell).fit()
aov_table = sm.stats.anova_lm(model, typ=2)
aov_table = format_anova_tbl(aov_table)
print(aov_table)

comp = mc.MultiComparison(per_cell['norm_value'], per_cell['Diagnosis'])
post_hoc_res = comp.tukeyhsd()
tukey_tbl = post_hoc_res.summary()
print(tukey_tbl)
