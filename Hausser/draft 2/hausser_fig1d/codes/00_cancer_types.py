"""Registry of the 8 cancer types shown in Hausser et al. 2019 Fig. 1d.

Hausser et al. tested all 15 TCGA cancer types with >=250 primary tumor
samples (BLCA, BRCA, CESC, COAD, HNSC, KIRC, LGG, LIHC, LUAD, LUSC, OV,
PRAD, STAD, THCA, UCEC; see Methods / Data availability). Six were
significant at FDR<10% (breast, colon, thyroid, bladder, low-grade glioma,
liver) and two more were borderline (lung p=0.01, head&neck p=0.02).
Fig. 1d shows exactly these eight:

    thyroid, bladder, liver, colon, glioma, breast, lung, head & neck

Breast uses METABRIC (not TCGA-BRCA) per the Fig. 1d caption
("TCGA, breast cancer from Metabric").

This module is just a lookup table; it does not touch the network or the
filesystem. `xena_cohort` is the UCSC Xena "TCGA Hub" dataset name (the
folder name you'll see after downloading from
https://xenabrowser.net/datapages/ -> TCGA Hub -> <cohort> -> "IlluminaHiSeq
RNASeqV2 (unc.edu, gene RSEM log2 -- HiSeqV2)"). Confirm exact names on the
portal before downloading; Xena has renamed some cohorts over time.

Cross-checked against the paper's Methods section (2026-08-24): the 15
TCGA types tested (>=250 primary tumors) are listed there by exact TCGA
disease code, and LUAD/LUSC are listed as two *separate* entries -- so
Fig. 1d's "Lung" panel is definitely one specific one of the two, not a
pooled category. The main text still never states which one; that detail
is only in the Supplementary Methods PDF, which wasn't accessible from this
review. LUAD remains the default below pending that confirmation.
"""

from __future__ import annotations

CANCER_TYPES = {
    "THCA": {
        "label": "Thyroid",
        "source": "TCGA",
        "xena_cohort": "TCGA Thyroid Cancer (THCA)",
        "hausser_p": "p < 0.001",
        "fig1d_panel": True,
    },
    "BLCA": {
        "label": "Bladder",
        "source": "TCGA",
        "xena_cohort": "TCGA Bladder Cancer (BLCA)",
        "hausser_p": "p = 0.001",
        "fig1d_panel": True,
    },
    "LIHC": {
        "label": "Liver",
        "source": "TCGA",
        "xena_cohort": "TCGA Liver Cancer (LIHC)",
        "hausser_p": "p = 0.002",
        "fig1d_panel": True,
    },
    "COAD": {
        "label": "Colon",
        "source": "TCGA",
        "xena_cohort": "TCGA Colon Cancer (COAD)",
        "hausser_p": "p = 0.009",
        "fig1d_panel": True,
    },
    "LGG": {
        "label": "Glioma",
        "source": "TCGA",
        "xena_cohort": "TCGA Lower Grade Glioma (LGG)",
        "hausser_p": "p < 0.001",
        "fig1d_panel": True,
    },
    "BRCA_METABRIC": {
        "label": "Breast",
        "source": "METABRIC",
        "xena_cohort": None,  # cBioPortal, not Xena; see 01b_prepare_metabric_full.py
        "hausser_p": "p = 0.001",
        "fig1d_panel": True,
    },
    "LUAD_OR_LUSC": {
        "label": "Lung",
        "source": "TCGA",
        # Hausser lists both LUAD and LUSC among the 15 tested types but Fig. 1d
        # just says "Lung" - the paper's supplementary methods should say which
        # one (or a pooled LUAD+LUSC) was plotted. Default to LUAD; override
        # with --cancer-type LUSC on the command line if the SI says otherwise.
        "xena_cohort": "TCGA Lung Adenocarcinoma (LUAD)",
        "hausser_p": "borderline, p = 0.013",
        "fig1d_panel": True,
        "note": "Confirm LUAD vs LUSC against Hausser Supplementary Methods before treating this as final.",
    },
    "HNSC": {
        "label": "Head & Neck",
        "source": "TCGA",
        "xena_cohort": "TCGA Head and Neck Cancer (HNSC)",
        "hausser_p": "borderline, p = 0.022",
        "fig1d_panel": True,
    },
    # The other 7 of the 15 tested types (CESC, KIRC, OV, PRAD, STAD, UCEC,
    # and whichever of LUAD/LUSC isn't used above) appeared as "clouds"
    # without significant vertices and are NOT part of Fig. 1d. Left out of
    # this registry on purpose - add them back only if you want to also
    # reproduce Table 1 / Supplementary Fig. 1A in full.
}


def fig1d_types():
    """The 8 cancer-type codes that belong in the Fig. 1d reproduction."""
    return [code for code, info in CANCER_TYPES.items() if info.get("fig1d_panel")]


if __name__ == "__main__":
    for code in fig1d_types():
        info = CANCER_TYPES[code]
        print(f"{code:16s} {info['label']:14s} ({info['source']}, Hausser {info['hausser_p']})")
