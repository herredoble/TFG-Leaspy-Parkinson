import numpy as np
import pandas as pd
from leaspy_alg import leaspyFit, leaspyPersonalize, leaspyEstimate
from leaspy import Leaspy, Data, AlgorithmSettings, IndividualParameters, Plotter#, __watermark__

def getMaeLeaspy(df_test, markers, ascending, model, ip):
    markers_orig = [m + "_orig" for m in markers]
    y_test = df_test[markers_orig]
    MAE = np.full(len(markers), np.nan)
    
    y_test_fit, y_test_fit_orig = leaspyEstimate(df_test, markers, ascending, model, ip)
    
    for i, m in enumerate(markers_orig): 
        MAE[i] = np.nanmean(abs(y_test_fit_orig[:, i] - y_test[m]))
        
    return MAE
        

# =============================================================================
# df = pd.read_excel("./data_multi/preclinical_CAM_jul22_1_2/dd.xlsx")
# markers = ["RAVLT_forgetting", "RAVLT_immediate"]
# ascending = [True, False]
# model = Leaspy.load("./data_multi/preclinical_CAM_jul22_1_2/leaspy_btstrp_1.json")
# _, ip = leaspyPersonalize(df, markers, model)
# MAE = getMaeLeaspy(df, markers, ascending, model, ip)
# =============================================================================
