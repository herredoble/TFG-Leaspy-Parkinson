from scipy.stats import pearsonr
import pandas as pd
import numpy as np
from scipy.io import loadmat
from build_model_leaspy import buildModelLeaspy, getDiagnose, getDPS, getSubVisitMatrix

def getOnset(DPS_test, posterior_bayes, prior_bayes, likelihood_bayes, 
             ages_test, labels_test, btstrp):
    I = DPS_test.shape[1] #Number of subjects
    prob_MCI = np.full(DPS_test.shape, np.nan)
    ages_onset = np.full((btstrp, I), np.nan)
    
    #MCI probability
    for n in range(0, btstrp):
        dps_test = DPS_test[n, :, :] 
        for i in range(0, dps_test.shape[0]):
            dps_subj = dps_test[i, :]
            mask_values = ~np.isnan(dps_subj)
            prob_MCI[n, i, mask_values] = posterior_bayes[n](
                likelihood_bayes[n], prior_bayes[n], dps_subj[mask_values])[1]
            
    for n in range(1, btstrp + 1):
        dps_test_n = DPS_test[n-1, :, :]
        prob_MCI = np.vectorize(lambda x : np.nan if np.isnan(x) else 1.0 if x > 0.5 else 0.0)(prob_MCI)
        
        for i in range(0, I):
            if np.nansum(prob_MCI[n-1,i, :] > 0):
                mask_values = ~np.vectorize(lambda x : x != x)(ages_test[i,:]) # Nan != Nan
                age_subj = ages_test[i, mask_values]
                label_subj = prob_MCI[n-1, i, mask_values]
                idx = np.where(label_subj == 1)[0]
                
                if label_subj[-1]:
                    if idx[0] > 1:
                        ages_onset[n-1, i] = (age_subj[idx[0]] + age_subj[idx[0]-1])/2
                    else:
                        ages_onset[n-1, i] = age_subj[idx[0]] 
    
    true_onset = getTrueOnset(labels_test, ages_test)
    
    #Metrics
    ages_bsl = ages_test[:, 0]
    
    acc_scu = np.full((btstrp, 1), np.nan)
    acc_smci = np.full((btstrp, 1), np.nan)
    acc_pcu = np.full((btstrp, 1), np.nan)
    corr_MCI_age = np.full((btstrp, 1), np.nan)
    corr_MCI_reserve = np.full((btstrp, 1), np.nan)
    
    n_scu = np.sum(np.isnan(true_onset))
    n_smci = np.sum(true_onset == ages_bsl)
    
    for n in range(0, btstrp):
        acc_scu[n] = 100 * np.sum(np.isnan(true_onset) & np.isnan(ages_onset[n, :])) / n_scu
        acc_smci[n] = 100 * np.sum((ages_onset[n, :] == ages_bsl) & (true_onset == ages_bsl)) / n_smci
        
    n_pcu = np.sum(true_onset > ages_bsl)
    
    for n in range(0,btstrp):
        acc_pcu[n] = 100 * np.sum((ages_onset[n, :] > ages_bsl) & (true_onset > ages_bsl)) / n_pcu
        corr_data_age = pd.DataFrame(np.column_stack((ages_onset[n, :], true_onset[0, :]))).dropna()
        corr_MCI_age[n] = corr_data_age.corr().iloc[0,1]
        corr_data_reserve = np.array(pd.DataFrame(np.column_stack((ages_onset[n, :] - ages_bsl, 
                                                          true_onset[0, :] - ages_bsl))).dropna())
        corr_MCI_reserve[n], _ = pearsonr(corr_data_reserve[:, 0], corr_data_reserve[:, 1])
        
    return acc_scu, acc_smci, acc_pcu, corr_MCI_age, corr_MCI_reserve
        

def getTrueOnset(labels_test, ages_test):
    I = labels_test.shape[0]
    true_onset = np.full((1, I), np.nan)
    MCI = np.where(labels_test == "MCI", 1, 0)
    
    for i in range(0, I):
        if np.nansum(MCI[i,:] > 0):
            mask_values = ~np.vectorize(lambda x : x != x)(ages_test[i,:]) # Nan != Nan
            age_subj = ages_test[i, mask_values]
            label_subj = MCI[i, mask_values]
            idx = np.where(label_subj == 1)[0]
            
            if label_subj[-1]:
                if idx[0] > 1:
                    true_onset[0, i] = (age_subj[idx[0]] + age_subj[idx[0]-1])/2
                else:
                    true_onset[0, i] = age_subj[idx[0]] 
    
    return true_onset

# =============================================================================
# data = loadmat("./data_multi/preclinical_CAM_jul22_1_2/leaspy_btstrp_data.mat")
# DPS_test = data["DPS_test"]
# posterior_bayes = data["posterior_bayes"]
# prior_bayes = data["prior_bayes"]
# likelihood_bayes = data["likelihood_bayes"]
# btstrp = 2
# df_test = pd.read_csv("./data_multi/preclinical_CAM_jul22_1_2/df_test.csv")
# prob_MCI0 = pd.read_csv("./data_multi/preclinical_CAM_jul22_1_2/prob_MCI0.csv")
# prob_MCI1 = pd.read_csv("./data_multi/preclinical_CAM_jul22_1_2/prob_MCI1.csv")
# prob_MCI = np.full(DPS_test.shape, np.nan) 
# prob_MCI[0,:,:] = prob_MCI0
# prob_MCI[1,:,:] = prob_MCI1
# labels_test = getSubVisitMatrix(df_test, "DX")
# ages_test = getSubVisitMatrix(df_test, "TIME")
# getOnset(DPS_test, posterior_bayes, prior_bayes, likelihood_bayes, ages_test, labels_test, btstrp, prob_MCI)
# =============================================================================

