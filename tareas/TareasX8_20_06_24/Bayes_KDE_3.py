# -*- coding: utf-8 -*-
"""
Created on Thu Jun 20 19:57:40 2024

@author: cplatero
"""

import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde

# Construcción del DataFrame con algunos valores NaN en DPS
data = {
    'DPS': [1.2, 2.4, 1.5, 2.3, np.nan, 3.6, 2.8, 3.3, np.nan, 1.8,
            2.9, 3.4, 2.0, np.nan, 1.3, 2.1, 3.5, 3.8, 1.4, 2.2],
    'DX': ['CU', 'CU', 'CU', 'CU', 'CU', 'CU', 'CU', 'CU', 'CU', 'CU',
           'MCI', 'MCI', 'MCI', 'MCI', 'MCI', 'MCI', 'MCI', 'MCI', 'MCI', 'MCI']
}

df = pd.DataFrame(data)

# Manejo de valores NaN: eliminación de filas con NaN en la columna DPS
df = df.dropna(subset=['DPS'])

# Separar los datos según las clases
dps_cu = df[df['DX'] == 'CU']['DPS']
dps_mci = df[df['DX'] == 'MCI']['DPS']

# Estimación de densidad de kernel para cada clase
kde_cu = gaussian_kde(dps_cu)
kde_mci = gaussian_kde(dps_mci)

# Función para predecir la etiqueta
def predict(dps_value):
    p_cu = kde_cu(dps_value)
    p_mci = kde_mci(dps_value)
    
    if p_cu > p_mci:
        return 'CU'
    else:
        return 'MCI'

# Ejemplo de uso:
nuevo_dps = 2.5
prediccion = predict(nuevo_dps)
print(f"Para DPS = {nuevo_dps}, la predicción es: {prediccion}")
