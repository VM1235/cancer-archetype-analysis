"""PCA helpers. Fitting space and visualization space are kept separate."""

import numpy as np
from sklearn.decomposition import PCA


def fit_pca(X, n_components, random_state=0):
    """Fit PCA on samples × features.

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)
    n_components : int
    """
    pca = PCA(n_components=n_components, svd_solver="full", random_state=random_state)
    scores = pca.fit_transform(np.asarray(X, dtype=float))
    return pca, scores


def cumulative_variance(pca):
    return np.cumsum(pca.explained_variance_ratio_)


def align_pca_signs(new_scores, saved_scores, pca):
    """Flip PCA component signs so new_scores match saved_scores."""
    saved = np.asarray(saved_scores, dtype=float)
    new = np.asarray(new_scores, dtype=float)
    n = min(saved.shape[1], new.shape[1])
    signs = np.ones(pca.components_.shape[0])
    for j in range(n):
        if np.corrcoef(new[:, j], saved[:, j])[0, 1] < 0:
            signs[j] = -1
            pca.components_[j] *= -1
            new[:, j] *= -1
    return pca, new, signs


def inverse_transform_scores(pca, scores):
    """Map PC-space coordinates back to gene space."""
    return pca.inverse_transform(np.asarray(scores, dtype=float))
