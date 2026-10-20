import os
os.chdir(r"C:\Users\jorge\OneDrive\Documentos\Universidad\TFG neuro\demo\multileaspy")

import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import seaborn as sns
from leaspy import Leaspy, Data, AlgorithmSettings, IndividualParameters, Plotter#, __watermark__
from sklearn.linear_model import LinearRegression
from scipy import stats
import random
from percentile_interp import balanceSamples, percentile_interp, ecdf
from scipy.interpolate import interp1d



# =============================================================================
# df_train = pd.read_csv("./data/LACT_train_ID.csv")
# df_test = pd.read_csv("./data/LACT_test_ID.csv")
# df_grace = pd.read_csv("./data/grace_data.csv")
# df_grace["ID"] = df_grace["ID"].astype(str)
# 
# id_train = df_train["ID"]
# id_test = df_test["ID"]
# 
# df_train1 = pd.read_csv("./data/leaspy_train.csv")
# df_test1 = pd.read_csv("./data/leaspy_test.csv")
# 
# df_leaspy = pd.concat([df_train1, df_test1], axis= 0)
# 
# mask_train = df_leaspy["ID"].isin(id_train)
# mask_test = df_leaspy["ID"].isin(id_test)
# df_train2 = df_leaspy.loc[mask_train]
# df_test2 = df_leaspy.loc[mask_test]
# =============================================================================


def leaspyFit(df_train, markers, save= True, save_path= "./data/leaspyFit.json"):
    #Convert ID column to string
    df_train["ID"] = df_train["ID"].astype(str)
    #Set index for train and test data
    df_train = df_train.set_index(['ID', 'TIME'], verify_integrity=True).sort_index()
    #Leaspy model parameters
    leaspy_model = 'logistic'
    algo_settings = AlgorithmSettings('mcmc_saem', 
                                      n_iter=3000, 
                                      seed=0,
                                      progress_bar=True)
    #Create data structures for leaspy code
    data_train = Data.from_dataframe(df_train[markers])
    #Fit the data
    leaspy = Leaspy(leaspy_model, source_dimension=2, noise_model="gaussian_diagonal")
    leaspy.fit(data_train, algo_settings)
    #Save the model
    if(save):
        leaspy.save(save_path, indent= 2)
    #Return the model
    return leaspy
   
def leaspyPersonalize(df_test, markers, leaspy_model):
    #Convert ID column to string
    df_test["ID"] = df_test["ID"].astype(str)
    #Set index for train and test data
    df_test = df_test.set_index(['ID', 'TIME'], verify_integrity=True).sort_index()
    #Create data structures for leaspy code
    data_test = Data.from_dataframe(df_test[markers])
    #Personalize the data
    settings_personalization = AlgorithmSettings('scipy_minimize', progress_bar=True, use_jacobian=True)
    ip = leaspy_model.personalize(data_test, settings_personalization)
    #Convert ip to dataframe
    df_ip = ip.to_dataframe()
    #Resetea el índice
    df_test = df_test.reset_index()
    #Merge individual parameters to df_test
    df_test = df_test.merge(df_ip[["tau", "xi", "sources_0"]], on= "ID", how="left")
    #Leaspy estimated conversion time
    df_test["leaspy_estimation"] = df_test["AGE"] - df_test["tau"]
    #Return individual parameters
    return df_test, ip

def leaspyEstimate(df_test, markers, ascending, model, ip):
    ids = df_test["ID"].unique()
    estimation = np.full(df_test[markers].shape, np.nan)
    markers_orig = [m + "_orig" for m in markers]
    
    for id in ids: 
        pos = np.sum(~np.isnan(estimation[:, 0]))
        time_subj = df_test[df_test["ID"] == id]["TIME"].to_numpy()
        estimation[pos: pos + len(time_subj), :] = model.estimate({str(id):time_subj}, 
                                                       ip._individual_parameters)[str(id)]
        
    estimation_orig = np.zeros(estimation.shape)
    
    # scu - 0, mci - 1, pcu - 2
    group = pd.DataFrame([0 if (mci, conv) == (0,0) else 
             1 if mci == 1 and conv == 0 else 
             2 if mci == 0 and conv > 0 else None
                 for mci,conv in zip(df_test["group"], df_test["CONVERSION_TIME_DAYS"])])

    _, percent, point_feat = percentile_interp(
        (lambda df_test: 
             df_test[markers_orig].rename(columns=
                                  {col: i for i, col in
                                   enumerate(df_test[markers_orig].columns)}))(df_test), group)
    
    #Convert to original scale
    for j in range(estimation.shape[1]):
        interp_func = interp1d(percent[j], point_feat[j], kind='linear', fill_value='extrapolate')
        estimation_orig[:,j] = interp_func(estimation[:, j])
        
    return estimation, estimation_orig


#df = pd.read_excel("./data_multi/preclinical_CAM_jul22_1_2/dd.xlsx")
# =============================================================================
# markers = ["RAVLT_learning", "ADAS13", "CDRSB", "TRABSCOR"]
# ascending = [True, False]
# #model = Leaspy.load("model_parameters_LACT.json")
# df_trainp = pd.read_csv("LACT_train_ID.csv")
# df_trainp["ID"] = df_trainp["ID"].astype(str)
# model = leaspyFit(df,markers,True, "fit.json")
# _, ip = leaspyPersonalize(df, markers, model)
# ip_df = pd.DataFrame(ip._individual_parameters).T
# ip_df.to_excel("ip.xlsx")
# e, eo = leaspyEstimate(df, markers, ascending, model, ip)
# =============================================================================

