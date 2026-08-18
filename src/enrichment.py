"""Independent subtype labels, distance bins, and hypergeometric enrichment."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist
from scipy.stats import hypergeom, rankdata

MARKER_TFS = ("ASCL1", "NEUROD1", "POU2F3", "YAP1")
SUBTYPES = ("A", "A2", "N", "P", "Y")


def _bh_qvalues(p_values):
    p = np.asarray(p_values, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    q = ranked * n / (np.arange(n) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(q, 0, 1)
    return out


def spearman_average_clusters(expr, n_clusters=5):
    """Cluster samples (columns) with Spearman distance and average linkage.

    Parameters
    ----------
    expr : pandas.DataFrame
        genes × samples
    """
    X = expr.T.values.astype(float)
    ranks = np.apply_along_axis(rankdata, 1, X)
    dist = pdist(ranks, metric="correlation")
    Z = linkage(dist, method="average")
    cluster_ids = fcluster(Z, t=n_clusters, criterion="maxclust")
    return pd.Series(cluster_ids, index=expr.columns, name="cluster")


def _cluster_tf_means(expr, clusters):
    means = {}
    for cid in sorted(clusters.unique()):
        samples = clusters.index[clusters == cid]
        means[cid] = {tf: float(expr.loc[tf, samples].mean()) for tf in MARKER_TFS}
    return pd.DataFrame(means).T


def name_clusters_by_tfs(expr, clusters):
    """Map numeric clusters to {A, A2, N, P, Y} using marker TF means.

    N/P/Y are the clusters with the highest NEUROD1 / POU2F3 / YAP1.
    The remaining two (ASCL1-high) clusters become A (higher ASCL1) and A2.
    """
    tf_means = _cluster_tf_means(expr, clusters)
    remaining = set(tf_means.index)
    mapping = {}
    for label, tf in (("N", "NEUROD1"), ("P", "POU2F3"), ("Y", "YAP1")):
        cid = tf_means.loc[list(remaining), tf].idxmax()
        mapping[int(cid)] = label
        remaining.remove(cid)
    ascl1_order = tf_means.loc[list(remaining), "ASCL1"].sort_values(ascending=False)
    mapping[int(ascl1_order.index[0])] = "A"
    mapping[int(ascl1_order.index[1])] = "A2"
    named = clusters.map(mapping)
    named.name = "subtype"
    return named, mapping, tf_means


def name_clusters_by_overlap(clusters, author_labels):
    """Name our clusters by majority overlap with published labels."""
    aligned = author_labels.reindex(clusters.index)
    mapping = {}
    used = set()
    for cid in sorted(clusters.unique()):
        counts = aligned[clusters == cid].value_counts()
        for label in counts.index:
            if label not in used and label in SUBTYPES:
                mapping[int(cid)] = label
                used.add(label)
                break
    unnamed = [cid for cid in clusters.unique() if int(cid) not in mapping]
    leftover = [s for s in SUBTYPES if s not in used]
    for cid, label in zip(unnamed, leftover):
        mapping[int(cid)] = label
    named = clusters.map(mapping)
    named.name = "subtype"
    return named, mapping


def distance_bins(pc_scores, archetypes, n_bins=10):
    """Equal-count bins of Euclidean distance to each archetype.

    Returns arrays samples × archetypes with bin ids 0 (closest) .. n_bins-1.
    If n_samples is not divisible by n_bins, leftover samples go to the
    closest bins first (sizes differ by at most 1).
    """
    scores = np.asarray(pc_scores, dtype=float)
    arcs = np.asarray(archetypes, dtype=float)
    n_samples = scores.shape[0]
    if n_samples < n_bins:
        raise ValueError(f"n_samples={n_samples} < n_bins={n_bins}")
    base, rem = divmod(n_samples, n_bins)
    sizes = np.array([base + (1 if i < rem else 0) for i in range(n_bins)], dtype=int)
    edges = np.cumsum(sizes)
    bin_ids = np.empty((n_samples, arcs.shape[0]), dtype=int)
    distances = np.empty((n_samples, arcs.shape[0]), dtype=float)
    for j, arc in enumerate(arcs):
        d = np.linalg.norm(scores - arc, axis=1)
        distances[:, j] = d
        order = np.argsort(d, kind="mergesort")
        ranks = np.empty(n_samples, dtype=int)
        ranks[order] = np.arange(n_samples)
        bin_ids[:, j] = np.searchsorted(edges, ranks, side="right")
    return bin_ids, distances


def hypergeometric_enrichment(subtypes, bin_ids, fdr=0.1, subtype_levels=None):
    """Test subtype over-representation in each archetype × distance bin.

    Hypergeometric P(X >= k) with Benjamini-Hochberg FDR across all tests.
    """
    subtypes = pd.Series(subtypes)
    if subtype_levels is None:
        labels = [s for s in SUBTYPES if (subtypes == s).any()]
    else:
        labels = [s for s in subtype_levels if (subtypes == s).any()]
    n_samples, n_arcs = bin_ids.shape
    n_bins = int(bin_ids.max()) + 1
    rows = []
    for subtype in labels:
        K = int((subtypes == subtype).sum())
        for arc in range(n_arcs):
            for b in range(n_bins):
                in_bin = bin_ids[:, arc] == b
                n = int(in_bin.sum())
                k = int(((subtypes == subtype).values & in_bin).sum())
                # P(X >= k)
                p = float(hypergeom.sf(k - 1, n_samples, K, n)) if k > 0 else 1.0
                expected = n * (K / n_samples)
                fold = (k / n) / (K / n_samples) if K > 0 and n > 0 else np.nan
                rows.append(
                    {
                        "subtype": subtype,
                        "archetype": arc,
                        "bin": b,
                        "k_in_bin": k,
                        "bin_size": n,
                        "K_subtype": K,
                        "expected": expected,
                        "fold_enrichment": fold,
                        "p_value": p,
                    }
                )
    table = pd.DataFrame(rows)
    table["q_value"] = _bh_qvalues(table["p_value"].values)
    table["significant"] = table["q_value"] < fdr
    peak = (
        table.sort_values(["subtype", "archetype", "fold_enrichment"], ascending=[True, True, False])
        .groupby(["subtype", "archetype"], as_index=False)
        .first()[["subtype", "archetype", "bin"]]
        .rename(columns={"bin": "peak_bin"})
    )
    table = table.merge(peak, on=["subtype", "archetype"], how="left")
    table["sig_peak_at_bin0"] = (
        table["significant"] & (table["bin"] == 0) & (table["peak_bin"] == 0)
    )
    return table
