# # -*- coding: utf-8 -*-
# """
# Created on Fri Apr 12 09:49:58 2024

# @author: cplatero
# """

import numpy as np
from sklearn.metrics import accuracy_score, roc_auc_score, confusion_matrix

def get_MAE(matriz1, matriz2):
    # Calcular la diferencia absoluta entre las dos matrices
    diferencia_absoluta = np.abs(matriz1 - matriz2)   
    # Calcular el valor medio de la diferencia absoluta por columnas
    valor_medio_diff_abs = np.nanmean(diferencia_absoluta, axis=0)   
    return valor_medio_diff_abs


def get_metrics_class_bin(y_true, y_pred):
    # Calcular la matriz de confusión
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred).ravel()
    # Calcular la sensibilidad
    sensibilidad = tp / (tp + fn)
    # Calcular la especificidad
    especificidad = tn / (tn + fp)
    # Calcular la exactitud
    exactitud = accuracy_score(y_true, y_pred)
    # Calcular el área bajo la curva (AUC)
    auc = roc_auc_score(y_true, y_pred)
    return sensibilidad, especificidad, exactitud, auc
