"""PCHA fitting, explained-sample-variance, and t-ratio permutation tests.

PCHA is run in PC space. Volumes for the t-ratio are computed after
projecting onto the first (k-1) PCs, matching the ParTI / Hausser convention:
a k-vertex simplex is (k-1)-dimensional.
"""

from __future__ import annotations

import math

import numpy as np
from py_pcha import PCHA
from scipy.spatial import ConvexHull, QhullError


def fit_pcha(pc_scores, k, delta=0.1, conv_crit=1e-6, maxiter=500, n_retries=25):
    """Fit PCHA to samples × PCs.

    Parameters
    ----------
    pc_scores : ndarray, shape (n_samples, n_pcs)
    k : int
        Number of archetypes.

    Returns
    -------
    archetypes : ndarray, shape (k, n_pcs)
    S : ndarray, shape (n_samples, k)
        Mixture weights (rows sum to 1).
    varexpl : float
        Fraction of variance explained (ESV).
    """
    X = np.asarray(pc_scores, dtype=float).T  # (n_pcs, n_samples)
    last_err = None
    for _ in range(int(n_retries)):
        try:
            XC, S, C, SSE, varexpl = PCHA(
                X, noc=int(k), delta=delta, conv_crit=conv_crit, maxiter=maxiter
            )
            archetypes = np.asarray(XC).T
            weights = np.asarray(S).T
            return archetypes, weights, float(varexpl)
        except Exception as err:
            # py_pcha furthest_sum can pick index == n_samples and raises
            # a locally defined InitializationException (not importable).
            if "Initialization does not converge" in str(err) or isinstance(err, IndexError):
                last_err = err
                continue
            raise
    raise RuntimeError(f"PCHA failed after {n_retries} init retries: {last_err}")


def simplex_volume(archetypes):
    """MATLAB ParTI simplex volume: |det(arcs[:, :-1] - arcs[:, -1])| / (k-1)!."""
    arcs = np.asarray(archetypes, dtype=float).T  # (dims, k) like MATLAB Archs
    k = arcs.shape[1]
    reduced = arcs[:, :-1] - arcs[:, [-1]]
    return float(abs(np.linalg.det(reduced[: k - 1, : k - 1])) / math.factorial(k - 1))


def fit_pcha_best(pc_scores, k, n_init=50, delta=0.0, conv_crit=1e-6, maxiter=500):
    """Match ParTI: fit PCHA n_init times on the first (k-1) PCs, keep max volume."""
    scores = np.asarray(pc_scores, dtype=float)[:, : int(k) - 1]
    best = None
    best_vol = -np.inf
    n_ok = 0
    for _ in range(int(n_init)):
        try:
            archetypes, weights, varexpl = fit_pcha(
                scores, k, delta=delta, conv_crit=conv_crit, maxiter=maxiter
            )
            vol = simplex_volume(archetypes)
            n_ok += 1
            if np.isfinite(vol) and vol > best_vol:
                best_vol = vol
                best = (archetypes, weights, varexpl, vol)
        except RuntimeError:
            continue
    if best is None:
        raise RuntimeError(f"PCHA failed all {n_init} inits for k={k}")
    return best[0], best[1], best[2], best[3], n_ok


def esv_curve(pc_scores, k_values, delta=0.1, seed=0):
    """Fit PCHA for each k and return explained variance."""
    rng = np.random.default_rng(seed)
    # seed is reserved so later random inits stay reproducible if PCHA is wrapped
    _ = rng
    rows = []
    for k in k_values:
        archetypes, weights, varexpl = fit_pcha(pc_scores, k, delta=delta)
        rows.append(
            {
                "k": int(k),
                "esv": varexpl,
                "archetypes": archetypes,
                "S": weights,
            }
        )
    return rows


def t_ratio(pc_scores, archetypes):
    """Volume(archetype simplex) / volume(data convex hull) in (k-1)-D PC space."""
    k = archetypes.shape[0]
    n_vol = k - 1
    if n_vol < 2:
        raise ValueError("t-ratio needs k >= 3 so the simplex is at least 2-D")
    data = np.asarray(pc_scores, dtype=float)[:, :n_vol]
    arcs = np.asarray(archetypes, dtype=float)[:, :n_vol]
    hull_vol = ConvexHull(data).volume
    arc_vol = ConvexHull(arcs).volume
    if hull_vol <= 0:
        return np.nan
    return float(arc_vol / hull_vol)


def _shuffle_columns(pc_scores, rng):
    shuffled = np.array(pc_scores, dtype=float, copy=True)
    for j in range(shuffled.shape[1]):
        rng.shuffle(shuffled[:, j])
    return shuffled


def permutation_t_ratio(
    pc_scores,
    k,
    n_perm=100,
    delta=0.1,
    seed=0,
    observed_archetypes=None,
    verbose=False,
):
    """Permute each PC independently, refit PCHA, and compare t-ratios.

    p = fraction of shuffled t-ratios >= observed t-ratio.
    """
    rng = np.random.default_rng(seed)
    if observed_archetypes is None:
        observed_archetypes, _, _ = fit_pcha(pc_scores, k, delta=delta)
    observed = t_ratio(pc_scores, observed_archetypes)

    null = []
    n_fail = 0
    n_perm = int(n_perm)
    for i in range(n_perm):
        shuffled = _shuffle_columns(pc_scores, rng)
        try:
            archetypes, _, _ = fit_pcha(shuffled, k, delta=delta)
            null.append(t_ratio(shuffled, archetypes))
        except (QhullError, ValueError, RuntimeError):
            n_fail += 1
        if verbose and (i + 1) % 100 == 0:
            print(f"    k={k}: {i + 1}/{n_perm} shuffles", flush=True)

    null = np.asarray(null, dtype=float)
    null = null[np.isfinite(null)]
    if len(null) == 0:
        p_value = np.nan
    else:
        p_value = float(np.mean(null >= observed))
    return {
        "k": int(k),
        "t_ratio": observed,
        "p_value": p_value,
        "n_perm": int(n_perm),
        "n_success": int(len(null)),
        "n_fail": int(n_fail),
        "null_t_ratios": null,
        "archetypes": observed_archetypes,
    }
