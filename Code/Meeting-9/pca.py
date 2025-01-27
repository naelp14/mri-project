import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

def read_data(csv):
    data = pd.read_csv(csv)
    return data

def pca_percentage(data):
    pca = PCA()
    pca.fit(data)

    explained_variance_ratio = pca.explained_variance_ratio_
    cumulative_variance_ratio = np.cumsum(explained_variance_ratio)

    # Print the variance explained by each component
    for i, variance in enumerate(explained_variance_ratio):
        print(f"Principal Component {i+1}: {variance:.2%} of the variance")
    print(f"Cumulative Variance Explained: {cumulative_variance_ratio[-1]:.2%}")

def init_pca(data):
    data_without_label = data.drop('cluster', axis=1)
    pca = PCA(n_components=2)
    latent_pca = pca.fit_transform(data_without_label)
    return latent_pca

def plot_pca(name, pca_result, labels):
    plt.figure(figsize=(8, 6))
    scatter = plt.scatter(pca_result[:, 0], pca_result[:, 1], c=labels, cmap='viridis', s=30, alpha=0.7)
    plt.colorbar(scatter, label="Label")
    plt.title(name)
    plt.xlabel("PCA Component 1")
    plt.ylabel("PCA Component 2")
    plt.grid()
    plt.show()

if __name__ == "__main__":
    # note: This code won't be used anymore. We can just plot the latent space, so we don't really need PCA
    cg_latent = read_data('cg_latent.csv')
    pca_percentage(cg_latent)
    cg_pca = init_pca(cg_latent)
    plot_pca("PCA Components CG", cg_pca, cg_latent['cluster'].to_numpy())
    
    csf_latent = read_data('csf_latent.csv')
    pca_percentage(csf_latent)
    csf_pca = init_pca(csf_latent)
    plot_pca("PCA Components CSF", csf_pca, csf_latent['cluster'].to_numpy())

    pet_latent = read_data('pet_latent.csv')
    pca_percentage(pet_latent)
    pet_pca = init_pca(pet_latent)
    plot_pca("PCA Components PET", pet_pca, pet_latent['cluster'].to_numpy())

    mri_latent = read_data('mri_latent.csv')
    pca_percentage(mri_latent)
    mri_pca = init_pca(mri_latent)
    plot_pca("PCA Components MRI", mri_pca, mri_latent['cluster'].to_numpy())

    risk_latent = read_data('risk_latent.csv')
    pca_percentage(risk_latent)
    risk_pca = init_pca(risk_latent)
    plot_pca("PCA Components Risk", risk_pca, risk_latent['cluster'].to_numpy())

    mri_risk_latent = read_data('mri_risk_latent.csv')
    pca_percentage(mri_risk_latent)
    mri_risk_pca = init_pca(mri_risk_latent)
    plot_pca("PCA Components MRI-Risk", mri_risk_pca, mri_risk_latent['cluster'].to_numpy())