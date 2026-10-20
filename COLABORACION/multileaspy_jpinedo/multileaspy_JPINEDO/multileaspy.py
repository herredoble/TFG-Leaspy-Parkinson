import os
os.chdir(r"C:\Users\jorge\OneDrive\Documentos\Universidad\TFG neuro\demo\multileaspy")

from leaspy_alg import leaspyFit, leaspyPersonalize
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import seaborn as sns
from leaspy import Leaspy, Data, AlgorithmSettings, IndividualParameters, Plotter#, __watermark__
from sklearn.linear_model import LinearRegression
from scipy import stats
import random
from itertools import combinations
from build_model_leaspy import buildModelLeaspy, getDiagnose
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde
from scipy.io import loadmat
from collections import defaultdict
import matplotlib.pyplot as plt
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import train_test_split
from build_model_leaspy import buildModelLeaspy, getDiagnose, getDPS, getSubVisitMatrix
from leaspy import Leaspy, Data, AlgorithmSettings, IndividualParameters, Plotter#, __watermark__
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
import scipy.io
from build_leaspy import buildLeaspy, bayesTrain, getAUC


def multiLeaspy(df_train, btstrp):
    header = "./data_multi/preclinical_CAM_jul22_"
    vector_range = range(4, 4 + 1)
    features_set = [3,5,8,12]
    
    classes = ["CN", "MCI", "Dementia"]
    classes_auc = ["CN", "MCI"]
    
    markers_df = pd.DataFrame({
    'markers': ["RAVLT_forgetting", "RAVLT_immediate", "RAVLT_learning", 
                "RAVLT_perc_forgetting", "ADAS13", "FAQ", "MMSE", "CDRSB", 
                "LDELTOTAL", "mPACCdigit", "mPACCtrailsB", "TRABSCOR"],
    'ascending': [True, False, False, 
                  True, True, True, False, True, 
                  False, False, False, True],
    'id': [1, 2, 3, 
           4, 5, 6, 7, 8, 
           9, 10, 11, 12]
    })
    
    feat_sub = pd.DataFrame(columns = [
                                        [f"feat_{i}" for i in range(1,max(vector_range)+1)] + 
                                        [f"id_{i}" for i in range(1,max(vector_range)+1)] +
                                        [f"ascending_{i}" for i in range(1,max(vector_range)+1)] +
                                        ["dir"]
                                    ])
    
    #Results table     
    results_table = pd.DataFrame()
    
    columns = [f"feat_{i}" for i in range(1, max(vector_range) +1)] + \
                [f"id_{i}" for i in range(1, max(vector_range) +1)] + \
                ["Corr_MCI_age", "Corr_MCI_reserve", "acc_scu",
                 "acc_smci", "acc_pcu", "AUC_train", "AUC_test"]
    
    for col in columns:
        results_table[col] = []
        
    
    
    #Train-test split
    df["ID"] = df["ID"].astype(str)
    ids = df["ID"].unique()
    id_train, id_test = train_test_split(ids, test_size=0.2, random_state=42)
    id_train = np.sort(id_train.astype(int)).astype(str)
    id_test = np.sort(id_test.astype(int)).astype(str)
    df_train = df[df["ID"].isin(id_train)]
    df_test = df[df["ID"].isin(id_test)]
    id_train = df_train["ID"].unique()
    id_test = df_test["ID"].unique()
    ids = np.concatenate((id_train, id_test), axis=0)
    #Diagnose
    diagnose_train = getDiagnose(df_train)
    diagnose_test = getDiagnose(df_test)
    #Labels
    labels_train = getSubVisitMatrix(df_train, "DX")
    labels_test = getSubVisitMatrix(df_test, "DX")
    #Ages
    ages_train = getSubVisitMatrix(df_train, "TIME")
    ages_test = getSubVisitMatrix(df_test, "TIME")
    
    # Create table with feats and directories to save data
    for n in vector_range:
        ncomb = list(combinations(markers_df.loc[markers_df["id"].isin(features_set), "markers"],n))
        
        for comb in ncomb:
            feat = [item for item in comb]
            id = list(map(lambda x: markers_df[markers_df["markers"] == x]["id"].iloc[0],feat))
            ascending = list(map(lambda x: markers_df[markers_df["markers"] == x]["ascending"].iloc[0],feat))
            
            dir = [header + "_".join(map(str, id))]
            feat = feat + [np.nan] * (max(vector_range) - n)
            id = id + [np.nan] * (max(vector_range) - n)
            ascending = ascending + [np.nan] * (max(vector_range) - n)
            
            row = feat + id + ascending + dir

            feat_sub = feat_sub.append(pd.Series(row, index = feat_sub.columns), ignore_index=True)
    
    for idx, row in feat_sub.iterrows():
        dir = row["dir"]
        file = dir + "/dd.xlsx"
        
        if os.path.exists(file):
            continue
        
        else:
            if not os.path.exists(dir):
                os.makedirs(dir)
                
            feat = row.filter(regex="feat*")
            feat = list(feat.dropna().values)
            feat_orig = [f + "_orig" for f in feat]
            
            cols = ["ID", "AGE", "TIME", "EXAMDATE", "M", "DX", "D", 
                    "CONVERSION_TIME_DAYS", "CONVERT", "onset", "group"] + feat + feat_orig
            
            df_subs = df_train[cols]
            df_subs.to_excel(file, index = False)
        
    for idx, row in feat_sub.iterrows():
        dir = row["dir"]
        file_in = dir + "/dd.xlsx"
        save_model_path = dir + "/leaspy_btstrp"
        
        if os.path.exists(save_model_path + "_data.mat"):
            continue
        
        else:
            feat = row.filter(regex="feat*")
            feat = list(feat.dropna().values)
            asc = row.filter(regex="ascending*")
            asc = list(asc.dropna().values)
            ids = row.filter(regex="id_")
            ids = list(ids.dropna().values)
            dd = pd.read_excel(file_in)
            dd["ID"] = dd["ID"].astype(str)
            AUC_test, AUC_train, acc_scu, acc_smci, acc_pcu, corr_MCI_age, corr_MCI_reserve = \
            buildLeaspy(dd, df_test, ids, id_train, id_test, 
                        diagnose_train, diagnose_test, labels_train,
                        labels_test, ages_train, ages_test, feat, asc, btstrp, 
                        classes, classes_auc, save_model_path)
            
            AUC_test_mean = np.mean(AUC_test)
            AUC_test_std = np.std(AUC_test)
            AUC_train_mean = np.mean(AUC_train)
            AUC_train_std = np.std(AUC_train)
            acc_scu_mean = np.mean(acc_scu)
            acc_scu_std = np.std(acc_scu)
            acc_smci_mean = np.mean(acc_smci)
            acc_smci_std = np.std(acc_smci)
            acc_pcu_mean = np.mean(acc_pcu)
            acc_pcu_std = np.std(acc_pcu)
            corr_MCI_age_mean = np.mean(corr_MCI_age)
            corr_MCI_age_std = np.std(corr_MCI_age)
            corr_MCI_reserve_mean = np.mean(corr_MCI_reserve)
            corr_MCI_reserve_std = np.std(corr_MCI_reserve)
            
            res = feat + ids + \
                [corr_MCI_age_mean, corr_MCI_reserve_mean, 
                 acc_scu_mean, acc_smci_mean, acc_pcu_mean,
                 AUC_train_mean, AUC_test_mean]
                
            res = dict(zip(results_table.columns, res))
            results_table = results_table.append(res, ignore_index=True)
            
    results_xlsx = "preclinical_CAM_jul22_" + "_".join(map(str, markers_df["id"])) + ".xlsx"
    results_table.to_excel(results_xlsx, index = False)
            
# =============================================================================
# 
#  Integrar los accuracies a la tabla de resultados final
#
# Tareas:
    # Mae = media del error absoluto entre y_test_real e y_test_fitted
    #   (creo que y es el valor de los marcadores no normalizados)
    # Nmae = media del error absoluto entre la desviación típica?
    # Hacer uno por cada bootstrap => media y desviación típica de ambos
    # Onset: calcular a partir de las predicciones (MCI si probabilidad por Bayes
    # es mayor que 0,5) y calcular el tiempo de onset según la visita en la cual
    # convierta (tengo la función). Comparar con la conversión real. 
    # Accuracies: a partir del onset, comparar si realmente pertenece al grupo
    # estimado.
# =============================================================================
def getResults():
    header = "./data_multi/preclinical_CAM_jul22_"
    nfeat = 2
    
    markers_df = pd.DataFrame({
    'markers': ["RAVLT_forgetting", "RAVLT_immediate", "RAVLT_learning", 
                "RAVLT_perc_forgetting", "ADAS13", "FAQ", "MMSE", "CDRSB", 
                "LDELTOTAL", "mPACCdigit", "mPACCtrailsB", "TRABSCOR"],
    'ascending': [True, False, False, 
                  True, True, True, False, True, 
                  False, False, False, True],
    'id': [1, 2, 3, 
           4, 5, 6, 7, 8, 
           9, 10, 11, 12]
    })
    
    feat_sub = pd.DataFrame(columns = [
                                        [f"feat_{i}" for i in range(1,nfeat+1)] + 
                                        [f"id_{i}" for i in range(1,nfeat+1)] +
                                        [f"ascending_{i}" for i in range(1,nfeat+1)] +
                                        ["dir"]
                                    ])
    
    # Create table with feats and directories to save data
    for n in range(2, nfeat + 1):
        ncomb = list(combinations(markers_df["markers"],n))
        
        for comb in ncomb:
            feat = [item for item in comb]
            id = list(map(lambda x: markers_df[markers_df["markers"] == x]["id"].iloc[0],feat))
            ascending = list(map(lambda x: markers_df[markers_df["markers"] == x]["ascending"].iloc[0],feat))
            
            dir = [header + "_".join(map(str, id))]
            feat = feat + [np.nan] * (nfeat - n)
            id = id + [np.nan] * (nfeat - n)
            ascending = ascending + [np.nan] * (nfeat - n)
            
            row = feat + id + ascending + dir

            feat_sub = feat_sub.append(pd.Series(row, index = feat_sub.columns), ignore_index=True)
            
    results_table = pd.DataFrame()
    
    columns = [f"feat_{i}" for i in range(1, nfeat +1)] + \
                [f"id_{i}" for i in range(1, nfeat +1)] + \
                ["Correlation", "tonset", "acc_scu",
                 "acc_smci", "acc_pcu_bl", "acc_pcu_end"]
    
    for col in columns:
        results_table[col] = []
    
    for idx, row in feat_sub.iterrows():
        file_in = row["dir"] + "/dd1.xlsx"
        dd1 = pd.read_excel(file_in)
        feats = row.filter(regex="feat.*")
        feats = list(feats.dropna().values)
        ids = row.filter(regex="id.*")
        ids = list(ids.dropna().values)
        
        times_table = pd.DataFrame(columns = ["ID", "min_time", "max_time", 
                                              "conv_time", "group"])
        
        for id in dd1["ID"].unique():
            filtered_dd1 = dd1[dd1["ID"] == id]
            min_time = (filtered_dd1["M"]/12 + filtered_dd1["leaspy_estimation"]).min()
            max_time = (filtered_dd1["M"]/12 + filtered_dd1["leaspy_estimation"]).max()
            conv_time = filtered_dd1["CONVERSION_TIME_DAYS"].iloc[0]
            group = "smci" if filtered_dd1["group"].iloc[0] == 1 else (
                        "pcu" if filtered_dd1["CONVERT"].iloc[0] == 1 else
                        "scu"
                    )
            
            times_table = times_table.append(pd.Series(
                [id, min_time, max_time, conv_time, group],
                index = times_table.columns), ignore_index=True)
            
        time_thresh = np.arange(-2, 8.2, 0.2)
        acc_table = pd.DataFrame(columns = ["tonset", "acc_scu", "acc_smci", 
                                            "acc_pcu_bl", "acc_pcu_end", "acc_sum"])
        scu_total = (times_table["group"] == "scu").sum()
        smci_total = (times_table["group"] == "smci").sum()
        pcu_total = (times_table["group"] == "pcu").sum()
        
        for thresh in time_thresh:
            acc_scu = ((times_table["group"] == "scu") & 
                       (times_table["max_time"] < thresh)).sum() / scu_total
            
            acc_smci = ((times_table["group"] == "smci") & 
                        (times_table["min_time"] > thresh)).sum() / smci_total
            
            acc_pcu_bl = ((times_table["group"] == "pcu") &
                          (times_table["min_time"] < thresh)).sum() / pcu_total
            
            acc_pcu_end = ((times_table["group"] == "pcu") &
                           (times_table["max_time"] > thresh)).sum() / pcu_total
            
            acc_sum = sum([acc_scu, acc_smci, acc_pcu_bl, acc_pcu_end])
            acc_tonset = thresh
            
            acc_table = acc_table.append(pd.Series(
                [acc_tonset, acc_scu, acc_smci, acc_pcu_bl, acc_pcu_end, acc_sum],
                index = acc_table.columns), ignore_index=True)
            
        #Optimum time onset    
        idx_opt = acc_table["acc_sum"].idxmax()
        #Time onset
        tonset = acc_table["tonset"].iloc[idx_opt]
        #Correlation
        pcu_ids = times_table[times_table["group"] == "pcu"]["ID"].values
        r = np.corrcoef(tonset - times_table[times_table["ID"].isin(pcu_ids)]["min_time"].astype(float),
                        times_table[times_table["ID"].isin(pcu_ids)]["conv_time"].astype(float))[0,1]
        #Results table
        res =  feats + ids + \
                [r, tonset,
                 acc_table["acc_scu"].iloc[idx_opt],
                 acc_table["acc_smci"].iloc[idx_opt],
                 acc_table["acc_pcu_bl"].iloc[idx_opt],
                 acc_table["acc_pcu_end"].iloc[idx_opt]]
                
        res = dict(zip(results_table.columns, res))
        
        results_table = results_table.append(res, ignore_index=True)
        
    results_xlsx = "preclinical_CAM_jul22_" + "_".join(map(str, markers_df["id"])) + ".xlsx"
    results_table.to_excel(results_xlsx, index = False)
    
    return results_table
    
    
                    

df = pd.read_csv("./data/prepared_data.csv")
# df["ID"] = df["ID"].astype(float)
# =============================================================================
# ids_ad = df[df["DX"] == "Dementia"]["ID"]
# df = df[~df["ID"].isin(ids_ad)]
# 
# df_trainp = pd.read_excel("./data/df_train.xlsx")
# df_testp = pd.read_excel("./data/df_test.xlsx")
# df_trainp["ID"] = df_trainp["ID"].astype(float)
# df_testp["ID"] = df_testp["ID"].astype(float)
# id_train = df_trainp["ID"].unique()
# id_test = df_testp["ID"].unique()
# 
# df_train = df[df["ID"].isin(id_train.astype(float))]
# df_test = df[df["ID"].isin(id_test)]
# =============================================================================

#df_train = df_train[df_train["ID"]!='6780']
# =============================================================================
# ids = df_train["ID"].unique()[:200]
# df_train = df_train[df_train["ID"].isin(ids)]
# =============================================================================
multiLeaspy(df, 1)
#tt = getResults()    

# =============================================================================
# #Train-test split
# df["ID"] = df["ID"].astype(str)
# ids = df["ID"].unique()
# id_train, id_test = train_test_split(ids, test_size=0.2, random_state=42)
# id_train = np.sort(id_train.astype(int)).astype(str)
# id_test = np.sort(id_test.astype(int)).astype(str)
# df_train = df[df["ID"].isin(id_train)]
# df_test = df[df["ID"].isin(id_test)]
# 
# dfp = pd.concat([df_train, df_test], axis=0)
# leaspy = Leaspy.load("./data_multi/preclinical_CAM_jul22_3_5_8_12/leaspy_btstrp1.json")
# _,ip = leaspyPersonalize(dfp, ["RAVLT_learning", "ADAS13", "CDRSB", "TRABSCOR"], leaspy)
# ip_df = pd.DataFrame(ip._individual_parameters).T
# ip_df.to_excel("ip.xlsx")
# 
# =============================================================================
