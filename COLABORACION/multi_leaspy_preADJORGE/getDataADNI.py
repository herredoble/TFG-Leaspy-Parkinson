import os
#os.chdir(r"C:\Users\jorge\OneDrive\Documentos\Universidad\TFG neuro\demo\multileaspy")

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

def daysFromBaseline(df):  
    df['D'] = 0
    for id in df['ID'].unique():
        subject = df[df['ID'] == id]
        baseline_date = pd.to_datetime(subject['EXAMDATE'].iloc[0], format="%Y/%m/%d")
        df.loc[df['ID'] == id, 'D'] = (pd.to_datetime(subject['EXAMDATE'], format="%Y/%m/%d") - baseline_date).dt.days
    
    return df

def conversionTimeDays(df):
    def computeTime(subject):
        dx_bl = subject['DX'].iloc[0]
        for ivisit in range(1, len(subject)):
            visit = subject.iloc[ivisit]
            if dx_bl != visit['DX']:
                return ((visit['D'] + subject.iloc[ivisit - 1]['D']) / 2) / 365
        return np.nan
    
    df['CONVERSION_TIME_DAYS'] = 0
    df['CONVERT'] = 0
    
    for id in df['ID'].unique():
        subject = df[df['ID'] == id]
        if len(subject) > 1 and (subject['DX'].iloc[0] != subject['DX'].iloc[-1]):
            df.loc[df['ID'] == id, 'CONVERSION_TIME_DAYS'] = computeTime(subject)
            df.loc[df['ID'] == id, 'CONVERT'] = 1
        else:
            df.loc[df['ID'] == id, 'CONVERT'] = 0
    
    return df


def centeredTimeDays(df):
    df['onset'] = 0
    for id in df['ID'].unique():
        subject = df[df['ID'] == id]
        dx_bl = subject['DX'].iloc[0]
        convert_bl = subject['CONVERT'].iloc[0]
        
        if dx_bl == "CN" and convert_bl == 1:  # pCU
            conversion_time_days = subject['CONVERSION_TIME_DAYS'].iloc[0] * 365
            df.loc[df['ID'] == id, 'onset'] = (df.loc[df['ID'] == id, 'D'] - conversion_time_days) / 30.4167
        
        elif dx_bl == "CN" and convert_bl == 0:  # sCU
            last_d = subject['D'].iloc[-1]
            df.loc[df['ID'] == id, 'onset'] = (df.loc[df['ID'] == id, 'D'] - last_d) / 30.4167
        
        elif (dx_bl in ["MCI", "Dementia"]) and convert_bl == 0:  # sMCI
            df.loc[df['ID'] == id, 'onset'] = df.loc[df['ID'] == id, 'D'] / 30.4167
    
    return df

def getDataADNI(path, feats, asc):
    df = pd.read_csv(path)
    df.rename(columns={"RID" : "ID"}, inplace = True)
    df = df.dropna(subset=["AGE", "DX"])
# =============================================================================
#     ids_ad = df[(df["DX_bl"] == "AD") | (df["DX"] == "Dementia")]["ID"]
# =============================================================================
    ids_ad = df[(df["DX_bl"] == "AD") | ((df["VISCODE"] == "bl") & (df["DX"] == "Dementia"))]["ID"]
    df = df[~df["ID"].isin(ids_ad)]
# =============================================================================
#     ids_ad = df[df["DX"] == "Dementia"]["ID"]
#     df = df[~df["ID"].isin(ids_ad)]
# =============================================================================
    df["TIME"] = df["AGE"] + df["M"]/12
    
# =============================================================================
#     df["CONVERT_TIME"] = 0
#     df["CENSURE_TIME"] = 0
#     df["DX_bl"].replace("SMC","CN",inplace=True)
#     df["DX_bl"].replace(["LMCI","EMCI"], "MCI", inplace=True)
#     _,_ = estableProgressive(df)
# =============================================================================
    
    df = daysFromBaseline(df)
    df = conversionTimeDays(df)
    df = centeredTimeDays(df)
    
    # cn --> 0, mci --> 1
    df["group"] = 1
    ids_cn = df.loc[(df["DX"] == "CN") & (df["M"] == 0), "ID"]
    df.loc[df["ID"].isin(ids_cn), "group"] = 0
    
    ids_pmci = df[(df["group"] == 1) & (df["CONVERT"] == 1)]["ID"]
    df = df[~df["ID"].isin(ids_pmci)]
    
    for i, feat in enumerate(feats):
        if asc[i] == False:
            df[feat] = -df[feat]
    
    group = pd.DataFrame([0 if (mci, conv) == (0,0) else 
             1 if mci == 1 and conv == 0 else 
             2 if mci == 0 and conv > 0 else None
                 for mci,conv in zip(df["group"], df["CONVERSION_TIME_DAYS"])])
    
    feat2per, percent, point_feat = percentile_interp(
        (lambda df: 
             df[feats].rename(columns=
                                  {col: i for i, col in enumerate(df[feats].columns)}))(df), 
        group)
    
    cols = ["ID", "TIME", "AGE", "DX", "EXAMDATE", "M", "D", "group", 
            "CONVERT", "CONVERSION_TIME_DAYS", "onset"] +  feats
    feats_orig = [feat + "_orig" for feat in feats]
    feats_dict = dict(zip(feats, feats_orig))
    df = df[cols]
    df.rename(columns = feats_dict, inplace = True)
    df = df.reset_index(drop=True)
    df[feats] = feat2per
    
    n=2
    num_visits = df.groupby("ID").size().reset_index(name="num_visits")
    ids_few = num_visits[num_visits["num_visits"] < n]["ID"].tolist()
    
    df = df[~df["ID"].isin(ids_few)]
    
    return df
    
    
    
