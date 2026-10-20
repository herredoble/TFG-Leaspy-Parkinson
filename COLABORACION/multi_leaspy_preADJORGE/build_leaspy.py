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
from leaspy_alg import leaspyFit, leaspyPersonalize, leaspyEstimate
import scipy.io
from getMae import getMaeLeaspy
from getOnset import getOnset, getTrueOnset

def buildLeaspy(df_train, df_test, id_train, id_test, diagnose_train, diagnose_test, labels_train,
                labels_test, ages_train, ages_test, markers, ascending, btstrp, classes, classes_auc,
                save_model_path):
       
    #DPM train
    Id_train, Idx_train, Tau_train, Xi_train, DPS_train = buildModelLeaspy(
        df_train, markers, btstrp, diagnose_train, classes, save_model_path)
    
    #Test personalize
    n_test = len(id_test)
    num_points_test = np.array([sum(df_test["ID"] == id) for id in id_test])
    max_visits_test = num_points_test.max()
    Xi_test = np.zeros((btstrp, n_test), dtype=float)
    Tau_test = np.zeros((btstrp, n_test), dtype=float)
    DPS_test = np.full((btstrp, n_test, max_visits_test), np.nan) 
    MAE_test = np.full((btstrp, len(markers)), np.nan)
    
    for n in range(1, btstrp + 1):
        path = save_model_path + str(n) + ".json"
        model = Leaspy.load(path)
        _, ip_test_orig = leaspyPersonalize(df_test, markers, model)
        ip_test = pd.DataFrame(ip_test_orig._individual_parameters).T
        Tau_test[n-1, :] = ip_test.loc[id_test.astype(str)]["tau"]
        Xi_test[n-1, :] = ip_test.loc[id_test.astype(str)]["xi"]
        
        #MAE
        MAE_test[n-1, :] = getMaeLeaspy(df_test, markers, ascending, model, ip_test_orig)
    
        for i, id in enumerate(id_test): 
            tij = df_test[df_test["ID"] == id]["TIME"]
            tau = Tau_test[n-1, i]
            xi = Xi_test[n-1, i]
            dps = getDPS(tij, tau, xi)
            DPS_test[n-1, i, : num_points_test[i]] = dps 
    
    #Bayesian classifier train
    posterior_bayes = np.empty(btstrp, dtype='O')
    likelihood = np.empty(btstrp, dtype='O')
    prior_bayes = np.empty(btstrp, dtype='O')
    evidence = np.empty(btstrp, dtype='O')

    for n in range(1, btstrp + 1):
        posterior_bayes[n-1], likelihood[n-1], prior_bayes[n-1], evidence[n-1] = bayesTrain(
            DPS_train[n-1, :, :], labels_train[Idx_train[n-1, :], :], classes)
    
    #AUC test
    AUC_test = np.zeros(btstrp, dtype = float)
    for n in range(1, btstrp + 1):
         AUC_test[n-1] = getAUC(DPS_test[n-1, :, :], posterior_bayes[n-1], 
                           likelihood[n-1], prior_bayes[n-1], labels_test, classes_auc)
         
    #AUC train
    AUC_train = np.zeros(btstrp, dtype = float)
    for n in range(1, btstrp + 1):
         AUC_train[n-1] = getAUC(DPS_train[n-1, :, :], posterior_bayes[n-1], 
                           likelihood[n-1], prior_bayes[n-1], labels_train[Idx_train[n-1, :], :],
                           classes_auc)
    
    #Metrics
    acc_scu, acc_smci, acc_pcu, corr_MCI_age, corr_MCI_reserve = getOnset(DPS_test, posterior_bayes, prior_bayes, likelihood, 
                                                                          ages_test, labels_test, btstrp)
    
    #Save data
    data = {"AUC_train" : AUC_train,
            "Id_train" : Id_train,
            "Idx_train" : Idx_train,
            "Tau_train" : Tau_train,
            "Xi_train" : Xi_train,
            "DPS_train" : DPS_train,
            "labels_train" : labels_train,
            "ages_train" : ages_train,
            "AUC_test" : AUC_test,
            "id_test" : id_test,
            "Tau_test" : Tau_test,
            "Xi_test" : Xi_test,
            "DPS_test" : DPS_test,
            "labels_test" : labels_test,
            "ages_test" : ages_test,
            "btstrp" : btstrp, 
            "acc_scu" : acc_scu,
            "acc_smci" : acc_smci,
            "acc_pcu" : acc_pcu,
            "corr_MCI_age" : corr_MCI_age,
            "corr_MCI_reserve" : corr_MCI_reserve,
            "MAE_test" : MAE_test
            }
    
    scipy.io.savemat(save_model_path + "data.mat", data)
         
    return AUC_test, AUC_train, acc_scu, acc_smci, acc_pcu, corr_MCI_age, corr_MCI_reserve, MAE_test

def bayesTrain(dps_train, labels, unique_labels):
    def evidence (likelihood, prior, points):
        ev = np.zeros_like(points).astype(np.float64)
        for i, (_,model) in enumerate(likelihood.items()):
            ev += model(points) * prior[i]
        return ev

    def posterior(likelihood, prior, points):
        ev = evidence(likelihood, prior, points)
        post = dict()
        
        for i, (_,model) in enumerate(likelihood.items()):
            post[i] = (prior[i] * likelihood[i](points)) / ev
            
        return post
   
    labels = labels.ravel()
    mask = np.vectorize(lambda x : x != x)(labels) # Nan != Nan
    labels = labels[~mask]
    dps_train = dps_train.ravel()
    dps_train = dps_train[~np.isnan(dps_train)]
    
    
    
    C = len(unique_labels)
    prior = np.zeros(C)
    likelihood = {}

    for c in range(C):
        prior[c] = np.sum(labels == unique_labels[c]) / len(labels)
        model = gaussian_kde(dps_train[labels == unique_labels[c]].astype(float), bw_method='scott') #KDE model
        likelihood[c] = model #Probability in point x is likelihood[c](x)

    return posterior, likelihood, prior, evidence

def getAUC(dps_test, posterior_bayes, likelihood, prior, labels, classes):
    dps_test = dps_test.ravel()
    mask = ~np.isnan(dps_test)
    dps_test = dps_test[mask]
    
    #One vs one AUC
    labels = labels.ravel()
    mask = np.vectorize(lambda x : x != x)(labels) # Nan != Nan
    labels = labels[~mask]
    mask_class = np.isin(labels, classes)
    
    labels = labels[mask_class]
    dps_test = dps_test[mask_class]
    
    predict_scores = posterior_bayes(likelihood, prior, dps_test)
    predict_scores = np.column_stack(list(predict_scores.values()))
    
    #Convert labels to numeric
    mapping = {nlabel: i for i, nlabel in enumerate(classes)}
    labels = np.array([mapping[nlabel] for nlabel in labels])
    
    if len(classes) == 2:
        #Binary case: the score corresponds to the probability of the class 
        #with the greater label (i.e. 2 if labels are 1 and 2)
        labels = -labels
        auc = roc_auc_score(
            labels,
            predict_scores[:, 0],
            average="macro"
        )
    else:
        #Multiclass case: the order of the class scores corresponds to the 
        #order of the labels, numerically (i.e. 1,2,3 of labels are 1,2,3)
        auc = roc_auc_score(
            labels,
            predict_scores,
            multi_class="ovo", 
            average="macro",
        )
    return auc

# =============================================================================
# df_train = pd.read_csv("./data/prepared_data.csv")
# df_train["ID"] = df_train["ID"].astype(str)
# ids_ad = df_train[df_train["DX"] == "Dementia"]["ID"]
# df_train = df_train[~df_train["ID"].isin(ids_ad)]
# ids = df_train["ID"].unique()[:200]
# df_train = df_train[df_train["ID"].isin(ids)]
# diagnose = getDiagnose(df_train)
# classes = ["CN", "MCI"]
# markers = ["ADAS13", "CDRSB", "TRABSCOR", "RAVLT_learning"]
# 
# df = df_train 
# df["ID"] = df["ID"].astype(str)
# ids = df["ID"].unique()
# id_train, id_test = train_test_split(ids, test_size=0.2, random_state=42)
# id_train = np.sort(id_train.astype(int)).astype(str)
# id_test = np.sort(id_test.astype(int)).astype(str)
# df_train = df[df["ID"].isin(id_train)]
# df_test = df[df["ID"].isin(id_test)]
# #Diagnose
# diagnose_train = getDiagnose(df_train)
# diagnose_test = getDiagnose(df_test)
# #Labels
# labels_train = getLabels(df_train)
# labels_test = getLabels(df_test)
# 
# #buildLeaspy(df_train, markers, 2, "./data_multi/preclinical_CAM_jul22_1_2/")
# buildLeaspy(df_train, df_test, ids, id_train, id_test, 
#             diagnose_train, diagnose_test, labels_train,
#             labels_test, markers, 2, classes, "./data_multi/preclinical_CAM_jul22_1_2/")
# =============================================================================
