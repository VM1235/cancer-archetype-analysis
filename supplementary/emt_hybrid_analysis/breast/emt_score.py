"""EMT scoring model from George, Jolly, Xu, Somarelli & Levine,
*Cancer Res* 77(22):6415-6428 (2017), doi:10.1158/0008-5472.CAN-16-3521.

Reproduces the paper's two-predictor ordinal multinomial logistic model
(their Equations 2-5). The model takes exactly two inputs per sample -
log2 CLDN7 expression, and log2(VIM) - log2(CDH1) as a stand-in for the
paper's "VIM/CDH1" predictor - and returns a continuous EMT score mu in
[0, 2]: 0 = epithelial (E), 1 = maximally hybrid E/M, 2 = mesenchymal (M).

This module does not refit the model. It plugs the paper's own fitted
coefficients (their Results/Discussion: beta = [-7.87, 0.0413, 1.36, -1.96],
reported order (alpha1, alpha2, beta1, beta2)) into their published
equations. It is meant to be applied to *new* expression matrices (e.g.
DepMap breast cell lines here), the same way the paper applies its model
to cell lines outside the NCI-60 training set (Table 2).

SIGN CONVENTION CAVEAT (read before trusting mu blindly)
----------------------------------------------------------
The paper's Equation 2 is:
    log(pi_jk / (1 - pi_jk)) = alpha_k - (beta1*X1 + beta2*X2)
i.e. the linear predictor is *subtracted*. Combined with the reported
coefficient signs, this literal reading implies that *higher* CLDN7
(an epithelial / tight-junction marker) pushes a sample toward the
*mesenchymal* end of the score - which is biologically backwards for
CLDN7. This kind of sign flip is a known hazard when re-implementing a
cumulative-logit ("ordinal") model from a paper: different software
(here, MATLAB's mnrfit) and different write-ups use opposite sign
conventions for the same fitted numbers, and it is not resolvable from
the text alone without the original fitted object.

We do not silently guess. `score_cell_lines()` computes mu under BOTH
sign conventions ("paper" = literal Eq. 2, and "flipped" = negated
linear predictor) and reports both. `sanity_check_against_labels()`
lets you pick the convention that agrees with an independent, trusted
biological label (e.g. PAM50: Basal/triple-negative and claudin-low
breast lines are well-established to be more mesenchymal/EMT-active
than Luminal lines). Use that function before reporting mu anywhere.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.special import expit

# (alpha1, alpha2, beta1, beta2), George et al. 2017, Results:
# "The best-fit model we ultimately utilized is completely described by
#  beta = [-7.87, 0.0413, 1.36, -1.96]"
DEFAULT_COEFFICIENTS = {
    "alpha1": -7.87,
    "alpha2": 0.0413,
    "beta1": 1.36,   # coefficient on X1 = log2(CLDN7)
    "beta2": -1.96,  # coefficient on X2 = log2(VIM/CDH1)
}

EMT_PREDICTOR_GENES = ("CLDN7", "CDH1", "VIM")

# George et al. Eq. 5 category thresholds (Methods, "EMT metric"):
MU_EPITHELIAL_MAX = 0.5
MU_MESENCHYMAL_MIN = 1.5


def build_predictors(expr, genes=EMT_PREDICTOR_GENES):
    """Extract X1 (log2 CLDN7) and X2 (log2 VIM - log2 CDH1) per sample.

    Parameters
    ----------
    expr : pandas.DataFrame
        genes x samples, already log2(TPM+1) (DepMap convention, matches
        what this repo's Panel A matrices already contain - no extra log).
    genes : tuple
        (CLDN7, CDH1, VIM) row names to pull. Raises KeyError with a
        clear message if any is absent from expr.index, rather than
        silently dropping the sample set.

    Returns
    -------
    pandas.DataFrame indexed by sample, columns
        CLDN7, CDH1, VIM, X1, X2
    """
    cldn7_name, cdh1_name, vim_name = genes
    missing = [g for g in genes if g not in expr.index]
    if missing:
        raise KeyError(
            f"EMT predictor genes missing from expression matrix: {missing}. "
            "Cannot compute George et al. 2017 mu without CLDN7/CDH1/VIM."
        )
    cldn7 = expr.loc[cldn7_name].astype(float)
    cdh1 = expr.loc[cdh1_name].astype(float)
    vim = expr.loc[vim_name].astype(float)
    out = pd.DataFrame(
        {
            "CLDN7": cldn7,
            "CDH1": cdh1,
            "VIM": vim,
            "X1": cldn7,
            "X2": vim - cdh1,  # log2(VIM/CDH1) since inputs are already log2
        }
    )
    out.index.name = "sample"
    return out


def emt_probabilities(x1, x2, coefficients=DEFAULT_COEFFICIENTS, sign="paper"):
    """George et al. Eqs 2-4: cumulative + categorical probabilities.

    Parameters
    ----------
    x1, x2 : array-like
        Predictor values (log2 CLDN7, log2 VIM/CDH1).
    sign : {"paper", "flipped"}
        "paper": eta = beta1*x1 + beta2*x2, pi_k = expit(alpha_k - eta)
                 (literal reading of Eq. 2).
        "flipped": pi_k = expit(alpha_k + eta). Use if the "paper"
                 convention fails the biological sanity check (see
                 module docstring / sanity_check_against_labels()).

    Returns
    -------
    pE, pH, pM : numpy arrays, P(E), P(hybrid E/M), P(M) per sample.
    """
    if sign not in ("paper", "flipped"):
        raise ValueError("sign must be 'paper' or 'flipped'")
    a1 = coefficients["alpha1"]
    a2 = coefficients["alpha2"]
    b1 = coefficients["beta1"]
    b2 = coefficients["beta2"]
    x1 = np.asarray(x1, dtype=float)
    x2 = np.asarray(x2, dtype=float)
    eta = b1 * x1 + b2 * x2
    if sign == "paper":
        pi1 = expit(a1 - eta)  # P(Y <= 1) = P(E)
        pi2 = expit(a2 - eta)  # P(Y <= 2) = P(E or E/M)
    else:
        pi1 = expit(a1 + eta)
        pi2 = expit(a2 + eta)
    # Guard against the (rare, only near-degenerate-coefficient) case
    # pi1 > pi2 from floating point at extreme eta; clip so pH >= 0.
    pi2 = np.maximum(pi2, pi1)
    pE = pi1
    pH = pi2 - pi1
    pM = 1.0 - pi2
    return pE, pH, pM


def emt_mu(pE, pH, pM):
    """George et al. Eq. 5: mu in [0, 2], 0=E, 1=maximal hybrid, 2=M."""
    pE = np.asarray(pE, dtype=float)
    pH = np.asarray(pH, dtype=float)
    pM = np.asarray(pM, dtype=float)
    mu = np.where(
        pE > pM,
        pH,
        np.where(pE < pM, 2.0 - pH, 1.0),
    )
    return mu


def classify_mu(mu):
    """Bin mu into {E, E/M, M} using the paper's Methods thresholds."""
    mu = np.asarray(mu, dtype=float)
    labels = np.where(
        mu < MU_EPITHELIAL_MAX,
        "E",
        np.where(mu > MU_MESENCHYMAL_MIN, "M", "E/M"),
    )
    return labels


def score_cell_lines(expr, genes=EMT_PREDICTOR_GENES, coefficients=DEFAULT_COEFFICIENTS):
    """Full pipeline: expr (genes x samples) -> tidy per-sample EMT table.

    Computes mu under BOTH sign conventions ("paper" and "flipped") so
    the caller can pick one with sanity_check_against_labels() instead
    of the module silently assuming one is correct.
    """
    pred = build_predictors(expr, genes=genes)
    out = pred.copy()
    for sign in ("paper", "flipped"):
        pE, pH, pM = emt_probabilities(pred["X1"], pred["X2"], coefficients, sign=sign)
        mu = emt_mu(pE, pH, pM)
        cat = classify_mu(mu)
        out[f"pE_{sign}"] = pE
        out[f"pH_{sign}"] = pH
        out[f"pM_{sign}"] = pM
        out[f"mu_{sign}"] = mu
        out[f"category_{sign}"] = cat
    return out


def sanity_check_against_labels(scored, labels, mesenchymal_label, epithelial_labels):
    """Pick the sign convention whose mu agrees with an independent label.

    Parameters
    ----------
    scored : DataFrame from score_cell_lines(), indexed by sample.
    labels : Series indexed the same way, e.g. PAM50 subtype calls.
    mesenchymal_label : str
        A label expected, on independent biological grounds, to be more
        mesenchymal/EMT-active (e.g. "Basal" for breast).
    epithelial_labels : list[str]
        Labels expected to be more epithelial (e.g. ["LumB", "Her2"]).

    Returns
    -------
    dict with the chosen sign, both group means for both conventions,
    and the mean-difference used to decide. Does not mutate `scored`.
    """
    aligned = labels.reindex(scored.index)
    mes_mask = aligned == mesenchymal_label
    epi_mask = aligned.isin(epithelial_labels)
    report = {}
    for sign in ("paper", "flipped"):
        mu = scored[f"mu_{sign}"]
        mes_mean = float(mu[mes_mask].mean())
        epi_mean = float(mu[epi_mask].mean())
        report[sign] = {
            "mesenchymal_label_mean_mu": mes_mean,
            "epithelial_labels_mean_mu": epi_mean,
            "difference_mes_minus_epi": mes_mean - epi_mean,
            "agrees_with_biology": mes_mean > epi_mean,
        }
    if report["paper"]["agrees_with_biology"] and not report["flipped"]["agrees_with_biology"]:
        chosen = "paper"
    elif report["flipped"]["agrees_with_biology"] and not report["paper"]["agrees_with_biology"]:
        chosen = "flipped"
    else:
        # Both or neither agree: pick the larger, correctly-signed gap;
        # if neither is correctly signed, still report both, chosen=None.
        chosen = max(
            report, key=lambda s: report[s]["difference_mes_minus_epi"]
        )
        if not report[chosen]["agrees_with_biology"]:
            chosen = None
    report["chosen_sign"] = chosen
    return report
