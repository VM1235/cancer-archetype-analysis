"""Broad task families for corner names (Hallmark-based), shared by the summary figures. Assignment by programme, checked against
the cells at each corner (KRAS-down corners 96% epithelial; MYC corners median 77% S/G2M; EMT/matrix-type corners mostly mesenchymal)."""
FAM = {
 "Proliferation": ["cycling", "MYC / growth (translation)", "mTORC1 Signaling", "PI3K/AKT/mTOR  Signaling"],
 "Epithelial differentiation": ["KRAS-down (epithelial / differentiated)", "luminal / differentiation", "KRAS-up (squamous / differentiation)"],
 "Mesenchymal / matrix (EMT)": ["EMT / matrix", "angiogenesis / matrix", "coagulation / matrix-like", "muscle / matrix-like", "TGF-β response", "epithelial junction"],
 "Metabolism & hypoxia": ["metabolic", "hypoxia / glycolysis", "Xenobiotic Metabolism", "Fatty Acid Metabolism", "Adipogenesis", "Reactive Oxygen Species Pathway"],
 "Stress & immune signalling": ["inflammatory / stress response", "p53 / stress", "UV Response Dn", "Unfolded Protein Response", "IL-6/JAK/STAT3 Signaling", "IL-2/STAT5 Signaling", "interferon response", "stress response", "Allograft Rejection"],
}
FCOL = {"Proliferation": "#2b6cb0", "Epithelial differentiation": "#2f855a", "Mesenchymal / matrix (EMT)": "#c53030", "Metabolism & hypoxia": "#d69e2e", "Stress & immune signalling": "#7c3aed", "Unnamed": "#a0aec0", "Identity (one sample)": "#718096"}
FSHORT = {"Proliferation": "Prolif.", "Epithelial differentiation": "Epi-diff.", "Mesenchymal / matrix (EMT)": "Mes / ECM", "Metabolism & hypoxia": "Metab.", "Stress & immune signalling": "Stress", "Unnamed": "?"}
T2F = {t: f for f, ts in FAM.items() for t in ts}
def family(task): return T2F.get(task, "Unnamed")
