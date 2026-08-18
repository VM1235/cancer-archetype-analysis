"""Parametric ComBat (Johnson 2007), matching sva::ComBat with mod = ~1.

Used for Groves Panel C: batch = cell-line vs tumor, optional reference batch.
Does not implement the sva/fsva surrogate-variable path; Groves saved the
ComBat (bc2) matrix, not the fsva matrix.
"""

from __future__ import annotations

import numpy as np


def _aprior(gamma_hat):
    m = float(np.mean(gamma_hat))
    s2 = float(np.var(gamma_hat, ddof=1))
    return (2.0 * s2 + m * m) / s2


def _bprior(gamma_hat):
    m = float(np.mean(gamma_hat))
    s2 = float(np.var(gamma_hat, ddof=1))
    return (m * s2 + m * m * m) / s2


def _postmean(g_hat, g_bar, n, d_star, t2):
    return (t2 * n * g_hat + d_star * g_bar) / (t2 * n + d_star)


def _postvar(sum2, n, a, b):
    return (0.5 * sum2 + b) / (n / 2.0 + a - 1.0)


def _it_sol(sdat, g_hat, d_hat, g_bar, t2, a, b, conv=1e-4):
    n = float(sdat.shape[1])
    g_old = g_hat.copy()
    d_old = d_hat.copy()
    change = 1.0
    count = 0
    while change > conv:
        g_new = _postmean(g_hat, g_bar, n, d_old, t2)
        sum2 = np.sum((sdat - g_new[:, None]) ** 2, axis=1)
        d_new = _postvar(sum2, n, a, b)
        change = max(
            float(np.max(np.abs(g_new - g_old) / g_old)),
            float(np.max(np.abs(d_new - d_old) / d_old)),
        )
        g_old = g_new
        d_old = d_new
        count += 1
        if count > 100:
            break
    return g_new, d_new


def combat(dat, batch, ref_batch=None, par_prior=True):
    """Adjust genes × samples matrix for a known batch factor.

    Parameters
    ----------
    dat : ndarray, shape (n_genes, n_samples)
    batch : array-like of length n_samples
    ref_batch : optional batch level left unadjusted (sva ref.batch)
    par_prior : if True, parametric empirical Bayes (sva default)
    """
    dat = np.asarray(dat, dtype=float)
    batch = np.asarray(batch)
    if dat.shape[1] != batch.shape[0]:
        raise ValueError("dat columns must match batch length")
    if np.isnan(dat).any():
        raise ValueError("ComBat input contains NaN")

    levels = np.array(sorted(set(batch.tolist()), key=str))
    if ref_batch is not None and ref_batch not in set(levels.tolist()):
        raise ValueError(f"ref_batch {ref_batch!r} not in {levels.tolist()}")
    n_batch = len(levels)
    batches = [np.where(batch == lev)[0] for lev in levels]
    n_batches = np.array([len(idx) for idx in batches], dtype=float)
    n_array = dat.shape[1]

    # intercept-only biological model (Groves: mod0 <- model.matrix(~1, phen))
    design_cov = np.ones((n_array, 1))
    batchmod = np.zeros((n_array, n_batch))
    for i, idx in enumerate(batches):
        batchmod[idx, i] = 1.0
    design = np.hstack([batchmod, design_cov])

    if ref_batch is not None:
        ref = int(np.where(levels == ref_batch)[0][0])
    else:
        ref = None

    # cannot estimate batch + intercept without a constraint; drop last batch
    # column unless using a reference batch (then drop the reference column).
    drop = ref if ref is not None else (n_batch - 1)
    keep = [j for j in range(n_batch) if j != drop] + [n_batch]
    design_est = design[:, keep]

    n_gene = dat.shape[0]
    B = np.linalg.lstsq(design_est, dat.T, rcond=None)[0]  # (p, genes)
    # map coefficients back onto full batch+intercept layout
    grand = np.zeros((n_batch + 1, n_gene))
    grand[keep, :] = B
    if ref is not None:
        # reference batch additive effect is 0 by construction
        grand[ref, :] = 0.0
    else:
        grand[drop, :] = 0.0

    intercept = grand[n_batch, :]
    gamma_hat = grand[:n_batch, :]  # batch additive, genes in columns
    # residual variance after batch+intercept fit
    fitted = design @ grand
    res = dat - fitted.T
    var_pooled = np.sum(res ** 2, axis=1) / float(n_array - 1)
    var_pooled = np.maximum(var_pooled, 1e-12)

    stand_mean = intercept[:, None] * np.ones((1, n_array))
    s_data = (dat - stand_mean) / np.sqrt(var_pooled)[:, None]

    gamma_hat_s = np.vstack(
        [s_data[:, idx].mean(axis=1) for idx in batches]
    )
    delta_hat = np.vstack(
        [s_data[:, idx].var(axis=1, ddof=1) for idx in batches]
    )
    delta_hat = np.maximum(delta_hat, 1e-12)

    gamma_star = np.empty_like(gamma_hat_s)
    delta_star = np.empty_like(delta_hat)
    if par_prior:
        for i in range(n_batch):
            if ref is not None and i == ref:
                gamma_star[i] = 0.0
                delta_star[i] = 1.0
                continue
            g = gamma_hat_s[i]
            d = delta_hat[i]
            g_bar = float(np.mean(g))
            t2 = float(np.var(g, ddof=1))
            a = _aprior(d)
            b = _bprior(d)
            gs, ds = _it_sol(s_data[:, batches[i]], g, d, g_bar, t2, a, b)
            gamma_star[i] = gs
            delta_star[i] = ds
    else:
        gamma_star = gamma_hat_s.copy()
        delta_star = delta_hat.copy()
        if ref is not None:
            gamma_star[ref] = 0.0
            delta_star[ref] = 1.0

    bayes = np.empty_like(s_data)
    for i, idx in enumerate(batches):
        bayes[:, idx] = (s_data[:, idx] - gamma_star[i][:, None]) / np.sqrt(
            delta_star[i][:, None]
        )
    return bayes * np.sqrt(var_pooled)[:, None] + stand_mean
