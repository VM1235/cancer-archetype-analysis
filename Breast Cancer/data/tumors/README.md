# TCGA-BRCA tumors (Panel C)

Expected files from [UCSC Xena](https://xenabrowser.net/datapages/) (TCGA Hub, BRCA):

1. `HiSeqV2` — gene expression, already log2. **Gitignored** (~172 MB).  
2. `TCGA.BRCA.sampleMap-BRCA_clinicalMatrix` — clinical / IHC annotations.

`prepare_tcga_brca.py` keeps primary tumors and intersects genes with the Panel A matrix. IHC columns used for the figure:

- `breast_carcinoma_estrogen_receptor_status`
- `lab_proc_her2_neu_immunohistochemistry_receptor_status`
