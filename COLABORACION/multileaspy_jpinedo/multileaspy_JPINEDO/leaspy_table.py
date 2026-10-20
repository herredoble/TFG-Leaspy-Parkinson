import os
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import seaborn as sns
from leaspy import Leaspy, Data, AlgorithmSettings, IndividualParameters, Plotter#, __watermark__
from sklearn.linear_model import LinearRegression
from scipy import stats
import random
from scipy.interpolate import interp1d
from leaspy_alg import leaspyFit, leaspyPersonalize
from percentile_interp import balanceSamples, percentile_interp, ecdf

def leaspy_table(df_test, btstrp, timepoints, percent, point_feat):
    ids = df_test["ID"].unique()
    
    leaspy_matrix = pd.DataFrame(columns = [
                                            ["RID", "Age"] + 
                                            [col for col_group in zip(
                                                    [f"dps_{i}" for i in range(1,btstrp +1)],
                                                    [f"ADAS13_{i}" for i in range(1,btstrp +1)],
                                                    [f"CDRSB_{i}" for i in range(1,btstrp +1)],
                                                    [f"TRABSCOR_{i}" for i in range(1,btstrp +1)],
                                                    [f"RAVLT_{i}" for i in range(1,btstrp +1)]
                                                )
                                                for col in col_group]
                                        ],
            index = range(len(ids) * len(timepoints))
        )    
    
    for n in range(1, btstrp +1):
        leaspy_model = Leaspy.load("./data/leaspy_btstrp_" + str(n) + ".json")
        _, ip = leaspyPersonalize(df_test, leaspy_model) # Individual parameters
        
        for i, id in enumerate(ids):
            estimation = leaspy_model.estimate({str(id):timepoints}, ip._individual_parameters)[str(id)]
            
            for j in range(estimation.shape[1]):
                interp_func = interp1d(percent[j], point_feat[j], kind='linear', fill_value='extrapolate')
                estimation[:,j] = interp_func(estimation[:, j])
                
            estimation[:, -1] = -estimation[:,-1] #Ravlt learning
            
            tau = ip._individual_parameters[str(id)]['tau']
            xi = ip._individual_parameters[str(id)]['xi']
            DPS = getDPS(timepoints, tau, xi)
            
            leaspy_matrix.loc[len(timepoints) * i : len(timepoints) * (i+1) - 1, 
                              ['dps_' + str(n),'ADAS13_' + str(n), 'CDRSB_' + str(n),
                               'TRABSCOR_' + str(n), 'RAVLT_' + str(n)]] = \
                np.vstack([DPS, np.array([estimation[:, i] for i in range(0,4)])]).T 

    leaspy_matrix.loc[:, ['RID', 'Age']] = np.vstack([
        np.array([id for id in ids for _ in range(len(timepoints))]),
        np.array([point for _ in range(len(ids)) for point in timepoints])
        ]).T
    
    return leaspy_matrix
        
            

def getDPS(tij, tau, xi):
    return math.exp(xi)*(tij-tau)
        
df_train = pd.read_csv("./data/train.csv")
df_test = pd.read_csv("./data/test.csv")
df_test = df_test[df_test["ID"].isin([6,8])]

markers = ["ADAS13_orig", "CDRSB_orig", "TRABSCOR_orig", "RAVLT_learning_orig"]

group = pd.DataFrame([0 if (mci, conv) == (0,0) else 
         1 if mci == 1 and conv == 0 else 
         2 if mci == 0 and conv > 0 else None
             for mci,conv in zip(df_train["group"], df_train["CONVERSION_TIME_DAYS"])])

_, percent, point_feat = percentile_interp(
    (lambda df_train: 
         df_train[markers].rename(columns=
                              {col: i for i, col in enumerate(df_train[markers].columns)}))(df_train), 
    group)

leaspy_matrix=leaspy_table(df_test, 1, np.linspace(60,100,100), percent, point_feat)

leaspy_model = Leaspy.load("./data/leaspy_btstrp_" + str(1) + ".json")
_,ip = leaspyPersonalize(df_train, leaspy_model)
ip = pd.DataFrame(ip._individual_parameters).T
