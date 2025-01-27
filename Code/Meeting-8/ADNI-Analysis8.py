import pandas as pd
import tensorflow as tf
import numpy as np
from tensorflow.keras.models import Model
from tensorflow.keras.layers import Input, Dense
from tensorflow.keras.optimizers import Adam
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler, LabelEncoder
import matplotlib.pyplot as plt

def transform_data(columns, data):
    adata = data.copy()
    adata['ABETA'] = pd.to_numeric(adata['ABETA'].str.replace(">1700", "2500"), errors='coerce')
    adata['TAU'] = pd.to_numeric(adata['TAU'].str.replace(">1700", "2500"), errors='coerce')
    adata['PTAU'] = pd.to_numeric(adata['PTAU'].str.replace(">1700", "2500"), errors='coerce')
    adata['PTGENDER'] = adata['PTGENDER'].replace({'Male': 1, 'Female': 0})
    X = adata[columns].copy()
    y = adata['cluster'].copy()

    if y.dtype == 'object' or y.dtype.name == 'category':
        label_encoder = LabelEncoder()
        y = label_encoder.fit_transform(y)

    
    X = np.nan_to_num(X, nan=0.0)
    X = np.where(np.isinf(X), 0.0, X)
    
    return (X, y)

def train_model(columns, data, baseline_data, file_name):
    X, y = transform_data(columns, data)
    X_train, X_test, _, y_test =  train_test_split(X, y, test_size=0.2, random_state=42)
    # Define the autoencoder architecture
    input_dim = X_train.shape[1]  # Number of features
    encoding_dim = 4  # Dimension of the latent space

    input_layer = Input(shape=(input_dim,))
    encoded = Dense(16, activation='relu')(input_layer)
    encoded = Dense(8, activation='relu')(encoded)
    encoded = Dense(4, activation='relu')(encoded)
    latent_space = Dense(encoding_dim, activation='relu', name='latent_space')(encoded)

    decoded = Dense(4, activation='relu')(latent_space)
    decoded = Dense(8, activation='relu')(decoded)
    decoded = Dense(16, activation='relu')(decoded)
    output_layer = Dense(input_dim, activation='sigmoid')(decoded)

    autoencoder = Model(inputs=input_layer, outputs=output_layer)

    # Encoder model to extract latent space
    encoder = Model(inputs=input_layer, outputs=latent_space)
    autoencoder.compile(optimizer=Adam(learning_rate=0.001), loss='mse')
    history = autoencoder.fit(
        X_train, X_train,
        epochs=50,
        batch_size=32,
        shuffle=True,
        validation_data=(X_test, X_test)
    )
    X_baseline, y_baseline = transform_data(columns, baseline_data)
    
    # Extract features from the latent space
    latent_features = encoder.predict(X_baseline)
    val_loss = history.history['val_loss']
    write_csv(latent_features, y_baseline, encoding_dim, file_name)
    plot_val_loss(val_loss)
    return autoencoder

def write_csv(latent_space, label, encoding_dim, file_name):
    latent_df = pd.DataFrame(latent_space, columns=[f'latent_{i+1}' for i in range(encoding_dim)])
    latent_df["cluster"] = label
    latent_df.to_csv(file_name, index=False)
    print(latent_space)

def plot_val_loss(val_loss):
    plt.plot(val_loss, label='Validation Loss')
    plt.xlabel('Epochs')
    plt.ylabel('Loss')
    plt.title('Training and Validation Loss')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    data = pd.read_csv('ADNIMERGE-Simple-dtw.csv')
    baseline_data = pd.read_csv('ADNIMERGE-dtw.csv')

    cg_columns = ["ADAS13", "FAQ", "MOCA", "TRABSCOR", "LDELTOTAL", "DIGITSCOR", "ADASQ4", "RAVLT_immediate", "ADAS11", "MMSE", "RAVLT_learning", "RAVLT_forgetting"]
    cg_model = train_model(cg_columns, data, baseline_data, 'cg_latent.csv')

    csf_columns = ["ABETA", "TAU", "PTAU"]
    csf_model = train_model(csf_columns, data, baseline_data, 'csf_latent.csv')

    pet_columns = ["FDG", "AV45"]
    pet_model = train_model(pet_columns, data, baseline_data, 'pet_latent.csv')

    risk_columns = ["AGE", "APOE4", "PTEDUCAT", "PTGENDER"]
    risk_model = train_model(risk_columns, data, baseline_data, 'risk_latent.csv')

    mri_columns = ["Hippocampus", "Entorhinal", "MidTemp", "Fusiform", "Ventricles", "WholeBrain"]
    mri_model = train_model(mri_columns, data, baseline_data, 'mri_latent.csv')

    mri_risk_columns = ["Hippocampus", "Entorhinal", "MidTemp", "Fusiform", "Ventricles", "WholeBrain", "AGE", "APOE4", "PTEDUCAT", "PTGENDER"]
    mri_risk_latent = train_model(mri_risk_columns, data, baseline_data, 'mri_risk_latent.csv')


