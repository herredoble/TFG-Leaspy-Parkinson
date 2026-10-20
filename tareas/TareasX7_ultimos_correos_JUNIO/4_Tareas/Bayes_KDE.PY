import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.neighbors import KernelDensity
from sklearn.metrics import accuracy_score, classification_report

# Crear un conjunto de datos de ejemplo
data = {
    'DPS': [2.3, 1.9, 3.1, 2.5, 3.8, 2.2, 3.0, 2.7, 3.2, 1.8,
            3.6, 2.1, 3.4, 2.4, 3.9, 2.0, 3.5, 2.8, 3.3, 1.7],
    'label': ['MCI', 'MCI', 'dementia', 'MCI', 'dementia', 'MCI', 'dementia', 'MCI', 'dementia', 'MCI',
              'dementia', 'MCI', 'dementia', 'MCI', 'dementia', 'MCI', 'dementia', 'MCI', 'dementia', 'MCI']
}

# Convertir el diccionario a un DataFrame
df_test = pd.DataFrame(data_test)
df_train = pd.DataFrame(data_train)

# Dividir los datos en características (X) y etiquetas (y)
X = df[['DPS']].values
y = df['label'].values

# Dividir los datos en conjuntos de entrenamiento y prueba
# X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=42)

def train_Bayes(X_train,y_train):
    # Separar los datos de entrenamiento por clase
    X_train_MCI = X_train[y_train == 'MCI']
    X_train_dementia = X_train[y_train == 'dementia']

    # Estimar la densidad de probabilidad utilizando KDE
    kde_MCI = KernelDensity(kernel='gaussian', bandwidth=0.2).fit(X_train_MCI)
    kde_dementia = KernelDensity(kernel='gaussian', bandwidth=0.2).fit(X_train_dementia)


    # Calcular las probabilidades a priori
    prior_MCI = len(X_train_MCI) / len(X_train)
    prior_dementia = len(X_train_dementia) / len(X_train)
    return kde_MCI, kde_dementia, prior_MCI, prior_dementia


def predict(X_test, kde_MCI, kde_dementia, prior_MCI, prior_dementia):
    # Calcular las log densidades de probabilidad
    log_dens_MCI = kde_MCI.score_samples(X_test)
    log_dens_dementia = kde_dementia.score_samples(X_test)

    # Calcular las probabilidades posteriores
    post_prob_MCI = np.exp(log_dens_MCI) * prior_MCI
    post_prob_dementia = np.exp(log_dens_dementia) * prior_dementia

    # Asignar la clase con la mayor probabilidad posterior
    y_pred = np.where(post_prob_MCI > post_prob_dementia, 'MCI', 'dementia')
    return y_pred
  
# Bayes  
  kde_MCI, kde_dementia, prior_MCI, prior_dementia = train_Bayes(X_train,y_train)
#Scores
y_pred_test = predict(X_test, kde_MCI, kde_dementia, prior_MCI, prior_dementia)
auc_test, report_test = get_SEN(y_pred_test,label_test)

y_pred_train = predict(X_train, kde_MCI, kde_dementia, prior_MCI, prior_dementia)
auc_train, report_train = get_SEN(y_pred_train,label_train)



# Evaluar el modelo
accuracy = accuracy_score(y_test, y_pred)
report = classification_report(y_test, y_pred)

print(f'Accuracy: {accuracy}')
print('Classification Report:')
print(report)
