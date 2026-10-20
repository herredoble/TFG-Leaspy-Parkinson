# https://disease-progression-modelling.github.io/pages/notebooks/disease_course_mapping/TP1_LMM.html
# https://gitlab.com/icm-institute/aramislab/disease-course-mapping-tutorials/-/tree/master/solutions/


import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import copy
from scipy import stats
#get_ipython().run_line_magic('matplotlib', 'inline')

sns.set_style('whitegrid')
pd.options.display.float_format = '{:.4g}'.format

from leaspy import Leaspy, IndividualParameters, AlgorithmSettings, Data, __watermark__
#from leaspy.io.logs.visualization.plotter import Plotter
from leaspy.io.logs.visualization.plotting import Plotting


# Watermark trace with all important packages versions
#__watermark__

# from leaspy.datasets import Loader
# df = Loader.load_dataset('parkinson-putamen-train_and_test')
 
# ### I.1. Prepare your data

data_path = os.path.join(os.getcwd(),'inputs')
df_train = pd.read_csv(os.path.join(data_path, 'LSDSMA_train.csv'))
newColum = pd.DataFrame({'SPLIT': ['train'] * df_train.shape[0]})
df_train = pd.concat([df_train, newColum], axis=1)

df_train['ID'] = df_train['ID'].astype(str)
df_train = df_train.set_index(['ID', 'TIME'], verify_integrity=True).sort_index()
df_train.head(15)

df_test = pd.read_csv(os.path.join(data_path, 'LSDSMA_test.csv'))
newColum = pd.DataFrame({'SPLIT': ['test'] * df_test.shape[0]})
df_test = pd.concat([df_test, newColum], axis=1)

df_test['ID'] = df_test['ID'].astype(str)
df_test = df_test.set_index(['ID', 'TIME'], verify_integrity=True).sort_index()
df_test.head(15)

df=pd.concat([df_train, df_test], axis=0)
df = df.set_index('SPLIT', append=True)

cols_norm=['LNS','SDMT','HVLTRDLY','SFT','MEANCAU','ABETA']
cols_orig=['LNS_orig','SDMT_orig','HVLTRDLY_orig','SFT_orig','MEANCAU_orig','ABETA_orig']

quantiles = [.05, .25, .5, .75, .95]
df.reset_index().describe(percentiles=quantiles).round(2).T


# plt.figure(figsize=(14,6))
# sns.scatterplot(data=df.xs('test', level='SPLIT').reset_index(),
#                 x='TIME', y='SDMT', alpha=.5, s=60)
# plt.title('SDMT - Raw data')
# plt.show()

######################################################################
# Linear regresion
# help(stats.linregress)
x = df_train.index.get_level_values('TIME')
y = df_train['SDMT']
result = stats.linregress(x, y)
slope, intercept = result.slope, result.intercept

df['Linear Regression'] = intercept + slope * df.index.get_level_values('TIME')

def compute_rmse(df: pd.DataFrame, model_name: str) -> float:
    """Compute RMSE between PUTAMEN column and the <model_name> column of df"""
    y = df['SDMT']
    y_hat = df[model_name]
    err = y_hat - y
    return (err ** 2).mean() ** .5

def store_rmse_train_test(df: pd.DataFrame, overall_results: pd.DataFrame, model_name: str) -> None:
    """Compute RMSE (per data split) and store it in <overall_results> dataframe"""
    for split_name, df_split in df.groupby('SPLIT'):
        overall_results.loc[model_name, split_name] = compute_rmse(df_split, model_name)

overall_results = pd.DataFrame({'train': [], 'test': []})
store_rmse_train_test(df, overall_results, 'Linear Regression')

overall_results


def get_title(overall_results: pd.DataFrame, model_name: str) -> str:
    """Precise model's name and its RMSE train & test"""
    rmse_train = overall_results.loc[model_name, 'train']
    rmse_test = overall_results.loc[model_name, 'test']
    title = f'SDMT - Raw data vs {model_name}\n'
    title += fr'$RMSE_{{train}}$ = {rmse_train:.1%}     $RMSE_{{test}}$ = {rmse_test:.1%}'
    return title

def plot_individuals(df: pd.DataFrame, overall_results: pd.DataFrame, model_name: str, **kwargs):
    # ---- Input manager
    kind = kwargs.get('kind', 'lines')
    sublist = kwargs.get('sublist', None)
    highlight_test = kwargs.get('highlight_test', True)
    ax = kwargs.get('ax', None)

    if ax is None:
        fig, ax = plt.subplots(figsize=(14, 6))
    else:
        plt.figure(figsize=(14, 6))

    # ---- Select subjects
    if sublist is None:
        sublist = df.index.unique('ID')

    # ---- If too many subject, do not display them in the legend
    display_id_legend = len(sublist) <= 10

    # ---- Plot
    if kind == 'scatter':
        # -- Used for Linear Regression
        sns.scatterplot(data=df.reset_index(), x='TIME', y='SDMT', alpha=.5,
                        s=60, label='Observations', ax=ax)
        ax.plot(df.index.get_level_values('TIME').values,
                df[model_name].values,
                label=model_name, c='C1')
        ax.legend(title='LABEL')
        if highlight_test:
            test = df.xs('test', level='SPLIT').loc[sublist].reset_index()
            sns.scatterplot(data=test, x='TIME', y='SDMT', legend=None, ax=ax)

    elif kind == 'lines':
        # -- Used for the other models
        # - Stack observations & reconstructions by the model
        df_stacked = df[['SDMT', model_name]].copy()
        df_stacked.rename(columns={'SDMT': 'Observations'}, inplace=True)
        df_stacked = df_stacked.stack().reset_index().set_index(['ID', 'SPLIT'])
        df_stacked.columns = ['TIME', 'LABEL', 'SDMT']

        # - Plot
        sns.lineplot(data=df_stacked.loc[sublist].reset_index(),
                     x='TIME', y='SDMT', hue='ID', style='LABEL',
                     legend=display_id_legend, ax=ax)
        if highlight_test:
            test = df.xs('test', level='SPLIT').loc[sublist].reset_index()
            sns.scatterplot(data=test, x='TIME', y='SDMT', hue='ID', legend=None, ax=ax)
        if display_id_legend:
            ax.legend(title='LABEL', bbox_to_anchor=(1.05, 1), loc='upper left')
    else:
        raise ValueError('<kind> input accept only "scatter" and "lines".'
                         f' You gave {kind}')

    ax.set_title(get_title(overall_results, model_name))
    return ax

# plot_individuals(df, overall_results, 'Linear Regression',
#                  kind='scatter', highlight_test=False)
# plt.show()



# plot_individuals(df, overall_results, 'Linear Regression',
#                  kind='lines', highlight_test=True)
# plt.title('Longitudinal aspect')
# plt.show()






###############################################################################
# Individual linear regression

def compute_individual_parameters(df_patient: pd.DataFrame) -> pd.Series:
    x = df_patient.index.get_level_values('TIME')
    y = df_patient['SDMT']
    # -- Individual Linear regression
    result = stats.linregress(x, y)
    return pd.Series({'INTERCEPT': result.intercept, 'SLOPE_TIME': result.slope})


subject_idx = '3003' # for example
compute_individual_parameters(df_train.loc[[subject_idx]])



individual_parameters = df_train.groupby('ID').apply(compute_individual_parameters)
individual_parameters.head()

def compute_individual_reconstruction(df: pd.DataFrame, individual_parameters: pd.DataFrame) -> pd.Series:

    df_with_individual_parameters = df.join(individual_parameters, on='ID')
    df_with_individual_parameters['TIME'] = df.index.get_level_values('TIME')

    # what is the column name for the slope? and what is the corresponding feature
    slope_column = individual_parameters.drop(columns='INTERCEPT').columns.item()
    slope_feature = slope_column.replace('SLOPE_', '')

    slope = df_with_individual_parameters[slope_column]
    intercept = df_with_individual_parameters['INTERCEPT']

    return intercept + slope * df_with_individual_parameters[slope_feature]

df['Individual Linear Regression'] = compute_individual_reconstruction(df, individual_parameters)
df.sample(10, random_state=1111).sort_index()

store_rmse_train_test(df, overall_results, 'Individual Linear Regression')
overall_results



visits_per_subjects = df_train.groupby('ID').count().sort_values('SDMT')
sublist = visits_per_subjects.tail(5).index.tolist()
sublist += visits_per_subjects.head(5).index.tolist()

plot_individuals(df, overall_results, 'Individual Linear Regression', sublist=sublist)
plt.show()



############################################################################
# lme
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.regression.mixed_linear_model import MixedLMParams

lmm = smf.mixedlm(formula='SDMT ~ 1 + TIME',
                  data=df_train.reset_index(),
                  groups="ID", re_formula="~ 1 + TIME").fit()
lmm.summary()



def rename_lmm_params(lmm_params: pd.DataFrame) -> pd.DataFrame:
    """Helper function to cope with different parameters names in LMM results."""
    return lmm_params.add_prefix('SLOPE_').rename({
        'SLOPE_ID': 'INTERCEPT', # for random effect
        'SLOPE_Intercept': 'INTERCEPT', # for fixed effect
    }, axis=1)

def fixed_parameters_of_mixed_model(lmm) -> pd.Series:
    return pd.concat([
        rename_lmm_params(lmm.fe_params),
        pd.Series({'noise_std': lmm.scale ** .5})
    ]).rename('Fixed parameters')

fe = fixed_parameters_of_mixed_model(lmm)
fe.to_frame()

df_random_effects = rename_lmm_params(pd.DataFrame.from_dict(lmm.random_effects, orient='index'))


def compute_individual_parameters_from_lmm(lmm):

    df_fixed_effects = rename_lmm_params(lmm.fe_params)
    df_random_effects = rename_lmm_params(pd.DataFrame.from_dict(lmm.random_effects, orient='index'))

    # INTERCEPT and SLOPE_XXX
    df_individual_parameters = df_fixed_effects + df_random_effects

    return df_individual_parameters

df_individual_params = compute_individual_parameters_from_lmm(lmm)
df_individual_params.head()

df['Linear Mixed Effect Model'] = compute_individual_reconstruction(df, df_individual_params) 
store_rmse_train_test(df, overall_results, 'Linear Mixed Effect Model')
overall_results

plot_individuals(df, overall_results, 'Linear Mixed Effect Model', sublist=sublist)
plt.show()

#############################################################

# ---- We use only train data to compute mean & std
ages = df.index.to_frame()['TIME']
ages_train = ages.xs('train', level='SPLIT')
age_mean = ages_train.mean()
age_std = ages_train.std()

df['TIME_NORMALIZED'] = (ages -  age_mean) / age_std

# ---- Re-create again df_train & df_test to get the new column
df_train = df.xs('train',level='SPLIT')
df_test = df.xs('test',level='SPLIT') 
lmm_2 = smf.mixedlm(formula='SDMT ~ TIME_NORMALIZED + 1', data=df_train.reset_index(),
                    groups="ID", re_formula="~TIME_NORMALIZED + 1").fit()

lmm_2.summary()

df_individual_params_2 = compute_individual_parameters_from_lmm(lmm_2)
df_individual_params_2.head() 

fe_2 = fixed_parameters_of_mixed_model(lmm_2).rename({'INTERCEPT': 'INTERCEPT_NORMALIZED'}, axis=1)

fe_2['INTERCEPT'] = fe_2['INTERCEPT_NORMALIZED'] - fe_2['SLOPE_TIME_NORMALIZED'] * age_mean / age_std
fe_2['SLOPE_TIME'] = fe_2['SLOPE_TIME_NORMALIZED'] / age_std

# Noise std-dev is decreased in v2 and mean slope of age is increased
pd.concat([fe.rename('v1'), fe_2.rename('v2')], axis=1).reindex(fe.index)

df['Linear Mixed Effect Model - V2'] = compute_individual_reconstruction(df, df_individual_params_2) 
store_rmse_train_test(df, overall_results, 'Linear Mixed Effect Model - V2')
overall_results


plot_individuals(df, overall_results, 'Linear Mixed Effect Model - V2', sublist=sublist)
plt.show()


models = ['Individual Linear Regression', 'Linear Mixed Effect Model - V2']

def plot_rmse_by_number_of_visit(models, df=df, *, group_nbr_of_visits: bool = True):
    rmse = df[models].copy()
    rmse.columns.names = ['MODEL']
    rmse = rmse.stack()
    rmse = (rmse - df['SDMT']) ** 2

    nbr_of_visits = rmse.xs('train', level='SPLIT').groupby(['MODEL','ID']).size()
    # we group number of visits together to form ~ balanced groups
    if group_nbr_of_visits:
        nbr_of_visits = pd.qcut(nbr_of_visits, q=13, duplicates='drop')

    rmse = rmse.xs('test',level='SPLIT')
    rmse = rmse.reset_index().rename(columns={0: 'RMSE'})
    rmse = rmse.reset_index()[['ID', 'MODEL', 'RMSE']]
    rmse = rmse.groupby(['MODEL','ID'], sort=False).mean() ** .5
    rmse['Visits number in train'] = nbr_of_visits

    plt.figure(figsize=(12, 5))
    ax = sns.boxplot(data=rmse.reset_index(), y='RMSE',
                     x='Visits number in train', hue='MODEL',
                     showmeans=True, whis=[5,95], showfliers=False)
    ax = sns.stripplot(data=rmse.reset_index(), y='RMSE',
                       x='Visits number in train', hue='MODEL',
                       dodge=True, alpha=.5, linewidth=1, ax=ax)

    plt.grid(True, axis='both')
    plt.ylim(0,None)
    plt.legend()
    return ax

plot_rmse_by_number_of_visit(models, df)
plt.show()


##########################################################

# --- Import methods
from leaspy import Leaspy, Data, AlgorithmSettings
from leaspy import Leaspy, Data, AlgorithmSettings, IndividualParameters

# --- Format the data
data = Data.from_dataframe(df_train[['SDMT']])

# --- Fit a model
leaspy_univariate = Leaspy('univariate_linear')
settings_fit = AlgorithmSettings('mcmc_saem',n_iter=3000, seed=0)
leaspy_univariate.fit(data, settings_fit)

settings_personalize = AlgorithmSettings('scipy_minimize', seed=0)
individual_parameters = leaspy_univariate.personalize(data, settings_personalize)

timepoints = df_train.index
df['Leaspy linear'] = leaspy_univariate.estimate(timepoints, individual_parameters)

store_rmse_train_test(df, overall_results, 'Leaspy linear')
overall_results

plot_individuals(df, overall_results, 'Leaspy linear', sublist=sublist)
plt.show()


###################################################################
leaspy_univariate = Leaspy('univariate_logistic')
settings_fit = AlgorithmSettings('mcmc_saem', n_iter=3000, seed=0)
leaspy_univariate.fit(data, settings_fit)

settings_personalize = AlgorithmSettings('scipy_minimize', seed=0)
individual_parameters = leaspy_univariate.personalize(data, settings_personalize)

timepoints = df_train.index
df['Leaspy logistic'] = leaspy_univariate.estimate(timepoints, individual_parameters)

store_rmse_train_test(df, overall_results, 'Leaspy logistic')
overall_results


plot_individuals(df, overall_results, 'Leaspy logistic', sublist=sublist)
plt.show()

models = ['Individual Linear Regression',
          'Linear Mixed Effect Model - V2',
          'Leaspy linear',
          'Leaspy logistic']

df_train = df.xs('train',level='SPLIT')
plot_rmse_by_number_of_visit(models, df_train)
plt.show()

#################################################################
# Bonus: have a look at the LMM log-likelihood landscape
import statsmodels.formula.api as smf
from tqdm import tqdm

import site
site.addsitedir('./utils')

from loglikelihood_landscape_lmm import plt_ll_landscape, plt_ll_landscape_

# new model (TIME_norm)
lmm_model = smf.mixedlm(formula="SDMT ~ 1 + TIME_NORMALIZED",
                        data=df.xs('train',level='SPLIT').reset_index(),
                        groups="ID", re_formula="~ 1 + TIME_NORMALIZED")

lmm_model.cov_pen = None

view = dict(
    vmin=2200,
    vmid=2500,
    vmax=2630,
    levels=[2500, 2600, 2620]
)

Rs = {} # placeholder for results

compute = True # uncomment to perform the computations ; disabled by default since computationally expensive

for fix_k, fix_v in tqdm([
    # corr fixed
    ('corr', -.99),
    ('corr', -.9),
    ('corr', -.75),
    ('corr', -.5),
    ('corr', -.25), # ~best
    ('corr', 0),
    ('corr', .5),
    ('corr', .75),
    ('corr', .99),
    # sx fixed
    ('sx', 1.),
    ('sx', 5),
    ('sx', 10), # ~best
    ('sx', 20),
    ('sx', 40),
    # sy fixed
    ('sy', 1),
    ('sy', 6), # ~best
    ('sy', 10),
    ('sy', 20),
]):
    if (fix_k, fix_v) in Rs:
        plt_ll_landscape_(Rs[(fix_k, fix_v)], **view)
    elif compute:
        Rs[(fix_k, fix_v)] = plt_ll_landscape(lmm_model, N=50, **{fix_k: fix_v}, **view)