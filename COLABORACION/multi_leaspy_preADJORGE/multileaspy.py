import os
os.chdir(r"C:\Users\rroca\Desktop\TFG\COLABORACION\multi_leaspy_preADJORGE")
#os.chdir(os.path.abspath(__file__))

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
from getDataADNI import getDataADNI


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
    
    # MAE table
    MAE_table = pd.DataFrame()
    
    MAE_columns = ["feat_id"] + \
              [m for m in markers_df[markers_df["id"].isin(features_set)]["markers"]]
        
    for col in MAE_columns:
        MAE_table[col] = []
    
    #Train-test split
    df["ID"] = df["ID"].astype(str)
    ids = df["ID"].unique()
    id_train, id_test = train_test_split(ids, test_size=0.2, random_state=42)
    id_train = np.sort(id_train.astype(int)).astype(str)
    id_test = np.sort(id_test.astype(int)).astype(str)
    df_train = df[df["ID"].isin(id_train)]
    df_test = df[df["ID"].isin(id_test)]
    id_test = df_test["ID"].unique()
    
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
            all_markers_nan_idx = df_subs[feat_orig].isna().all(axis=1)
            df_subs = df_subs[~all_markers_nan_idx]
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
            
            id_train = dd["ID"].unique()
            #Diagnose
            diagnose_train = getDiagnose(dd)
            diagnose_test = getDiagnose(df_test)
            #Labels
            labels_train = getSubVisitMatrix(dd, "DX")
            labels_test = getSubVisitMatrix(df_test, "DX")
            #Ages
            ages_train = getSubVisitMatrix(dd, "TIME")
            ages_test = getSubVisitMatrix(df_test, "TIME")
            
            AUC_test, AUC_train, acc_scu, acc_smci, acc_pcu, corr_MCI_age, corr_MCI_reserve, MAE_test = \
            buildLeaspy(dd, df_test, id_train, id_test, 
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
            MAE_test_mean = MAE_test.mean(axis=0)
            MAE_test_std = MAE_test.std(axis=0)
            
            res = [feat[i] if (len(feat) > i) else "" for i in range(0,max(vector_range))] + \
                  [ids[i] if (len(feat) > i) else "" for i in range(0,max(vector_range))] + \
                  [corr_MCI_age_mean, corr_MCI_reserve_mean, 
                   acc_scu_mean, acc_smci_mean, acc_pcu_mean,
                   AUC_train_mean, AUC_test_mean]
                
            res = dict(zip(results_table.columns, res))
            results_table = results_table.append(res, ignore_index=True)
            
            str_id = "_".join(map(str, ids))
            MAE_row = [str_id] + \
                      ["" if id not in ids else MAE_test_mean[ids.index(id)] for id in features_set]
            
            MAE_row = dict(zip(MAE_table.columns, MAE_row))
            MAE_table = MAE_table.append(MAE_row, ignore_index = True)
            
    results_xlsx = "preclinical_CAM_jul22_" + "_".join(map(str, markers_df["id"])) + ".xlsx"
    results_table.to_excel(results_xlsx, index = False)
    
    MAE_table.to_excel("MAE_test_preclinical_CAM_jul22.xlsx", index=False)
               
    
                    
feats = ["RAVLT_forgetting", "RAVLT_immediate", "RAVLT_learning", 
            "RAVLT_perc_forgetting", "ADAS13", "FAQ", "MMSE", "CDRSB", 
            "LDELTOTAL", "mPACCdigit", "mPACCtrailsB", "TRABSCOR"] 
asc = [True, False, False, 
              True, True, True, False, True, 
              False, False, False, True] 
df =  getDataADNI("./data/ADNIMERGE_220706.csv", feats, asc) 

#df = pd.read_csv("./data/prepared_data.csv")
multiLeaspy(df, 10)