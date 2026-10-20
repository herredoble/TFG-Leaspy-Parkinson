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
from leaspy_alg import leaspyFit, leaspyPersonalize


def buildModelLeaspy(df_train, markers, btstrp, diagnose, classes, save_model_path):
    ids = df_train["ID"].unique()
    n_train = len(ids)
    num_points = np.array([sum(df_train["ID"] == id) for id in ids])
    max_visits = num_points.max()
    Idx_train = np.zeros((btstrp, n_train), dtype=int)
    Id_train = np.zeros((btstrp, n_train), dtype=int)
    Xi = np.zeros((btstrp, n_train), dtype=float)
    Tau = np.zeros((btstrp, n_train), dtype=float)
    DPS_Train = np.full((btstrp, n_train, max_visits), np.nan) 
    
    for n in range(1, btstrp + 1):
        idx_btstrp = []
        np.random.seed(n)
        
        for u in range(len(classes)):
            for uu in range(len(classes)):
                idx_u = np.where((diagnose[0] == classes[u]) & (diagnose[1] == classes[uu]))[0]
                med_u = np.median(num_points[idx_u])
                idx_many = idx_u[num_points[idx_u] > med_u]
                
                if len(idx_many) > 0:
                    idx_rnd = np.random.choice(idx_many, size=len(idx_many), replace=True)
                    idx_btstrp.extend(idx_rnd)
                    
                idx_few = idx_u[num_points[idx_u] <= med_u]
                
                if len(idx_few) > 0:
                    idx_rnd = np.random.choice(idx_few, size=len(idx_few), replace=True)
                    idx_btstrp.extend(idx_rnd)
        
        Idx_train[n-1, :] = idx_btstrp
        Id_train[n-1, :] = ids[idx_btstrp]
        
        #Train
        path = save_model_path + str(n) + ".json"
        df_train_btstrp = df_train[df_train["ID"].isin(Id_train[n-1,:].astype(str))]
        
        if not os.path.exists(path):
            model = leaspyFit(df_train_btstrp, markers, save=True, save_path=path)
            Leaspy.save(model, path)
        else:
            model = Leaspy.load(path)
        
        #Personalize
        _, ip = leaspyPersonalize(df_train_btstrp, markers, model)
        ip = pd.DataFrame(ip._individual_parameters).T
        Tau[n-1, :] = ip.loc[Id_train[n-1, :].astype(str)]["tau"]
        Xi[n-1, :] = ip.loc[Id_train[n-1, :].astype(str)]["xi"]
        
        #DPS train
        for i, (idx, id) in enumerate(zip(Idx_train[n-1, :], Id_train[n-1, :])): 
            tij = df_train[df_train["ID"] == id.astype(str)]["TIME"]
            tau = Tau[n-1, i]
            xi = Xi[n-1, i]
            dps = getDPS(tij, tau, xi)
            DPS_Train[n-1, i, : num_points[idx]] = dps
        
    return Id_train, Idx_train, Tau, Xi, DPS_Train

    


def getDiagnose(df):
    ids = df["ID"].unique()
    diagnose = np.zeros((2, len(ids)), dtype= 'O')
    
    for idx, id in enumerate(ids):
        df_sub = df[df["ID"] == id]
        df_sub_dx = df_sub["DX"]
        diagnose[0,idx] = df_sub_dx.iloc[0]
        diagnose[1,idx] = df_sub_dx.iloc[-1]
        
    return diagnose

def getDPS(tij, tau, xi):
    return math.exp(xi)*(tij-tau)

def getSubVisitMatrix(df, col):
    ids = df["ID"].unique()
    dx_per_sub = df.groupby("ID")[col]
    ncol = dx_per_sub.count().max()
    
    col = np.full((len(ids), ncol), np.nan, dtype='O')
    
    for idx, id in enumerate(ids):
        dx_sub = dx_per_sub.get_group(id)
        col[idx, :len(dx_sub)] = dx_sub
        
    return col

# =============================================================================
# df_train = pd.read_csv("./data/prepared_data.csv")
# df_train["ID"] = df_train["ID"].astype(str)
# diagnose = getDiagnose(df_train)
# classes = ["CN", "MCI", "Dementia"]
# markers = ["ADAS13", "CDRSB", "TRABSCOR", "RAVLT_learning"]
# getLabels(df_train)
# buildModelLeaspy(df_train, markers, 2, diagnose, classes, "./data_multi/preclinical_CAM_jul22_1_2/model")
# 
# =============================================================================
