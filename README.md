# Pituitary Adenoma DEG Analysis

**Author:** Hafsah Shamsi · [github.com/HafsahShamsi](https://github.com/HafsahShamsi)

Differential gene expression analysis of pituitary adenoma vs normal pituitary tissue, with pathway enrichment and targeted analysis of HPA axis feedback genes.

---

## Background

Pituitary adenomas are benign tumours of the pituitary gland characterised by loss of negative feedback sensitivity. Normally, circulating hormones suppress their own production via the hypothalamus-pituitary axis (HPA). In adenoma cells, this feedback loop is disrupted — cells keep secreting hormones regardless of circulating levels.

This project maps the transcriptional signature of that feedback loss using publicly available microarray data.

---

## Dataset

**GSE26966** — GEO database  
23 samples · GPL570 (Affymetrix HG-U133 Plus 2.0)  
14 pituitary adenoma · 9 normal pituitary tissue

---

## Pipeline

| Script | What it does |
|---|---|
| `pituitary_deg.py` | Fetches GSE26966, builds expression matrix, computes DEGs (t-test + FDR correction), generates volcano plot, PCA, and HPA axis gene plot |
| `enrichment.py` | Runs Enrichr on upregulated and downregulated DEGs separately — KEGG 2021 and GO Biological Process 2021 |

---

## Key findings

**DEGs:** 1033 upregulated · 1556 downregulated (FDR < 0.05, |log2FC| > 1)

**PCA:** PC1 explains 61.9% of variance and cleanly separates tumor from normal — confirming a large-scale transcriptional shift.

**Most significant downregulated genes:**
- **POMC** (log2FC = −12.2) — proopiomelanocortin, ACTH precursor. Massive suppression of the core secretory gene even as downstream hormone output is dysregulated
- **POU1F1** (−6.9) — master pituitary transcription factor. Its loss reflects dedifferentiation of adenoma cells
- **NKX2-2** (−7.1) — developmental identity gene, lost in tumour
- **CDKN2A** (−3.5) — p16 tumour suppressor, loss drives unchecked proliferation

**HPA axis feedback genes — all downregulated in adenoma:**
- SSTR2, SSTR5 — somatostatin receptors (inhibitory feedback)
- NR3C1 — glucocorticoid receptor
- FKBP5 — cortisol feedback marker
- PRLR, GHR — hormone receptors

**Pathway enrichment (KEGG):**

Downregulated: PI3K-Akt (adj_p=3.25e−06), MAPK signalling (3.75e−05), ECM-receptor interaction — growth regulation and tissue architecture collapsed

Upregulated: Aldosterone synthesis, Dopaminergic/Glutamatergic/Cholinergic synapse, Circadian entrainment — neuroendocrine rewiring and loss of circadian HPA regulation

---

## Connection to neuroimmunology

Cortisol dysregulation from a malfunctioning pituitary has direct consequences for neuroinflammation. HPA axis dysfunction is documented in MS patients — chronically altered cortisol rhythms affect T-cell trafficking and CNS immune surveillance. This project sits at the intersection of endocrinology and neuroimmunology, extending the research direction established in [neuro-deg-scanner](https://github.com/HafsahShamsi/neuro-deg-scanner).

---

## Outputs

'''

outputs/
├── pituitary_deg_results.csv          # Full DEG table
├── volcano_pituitary.png              # Volcano plot
├── pca_pituitary.png                  # PCA — tumor vs normal
├── feedback_genes.png                 # HPA axis gene fold changes
├── enrichment_KEGG_2021_Human.png     # KEGG dot plot
└── enrichment_GO_Biological_Process_2021.png

'''

---

## Dependencies

pip install pandas numpy matplotlib scipy statsmodels GEOparse scikit-learn requests

---

## Context

Part of a broader computational neuroimmunology pipeline:

- [`neuro-deg-scanner`](https://github.com/HafsahShamsi/neuro-deg-scanner) — reusable DEG pipeline
- [`ms-pathway-enrichment`](https://github.com/HafsahShamsi/ms-pathway-enrichment) — KEGG/GO enrichment on MS DEGs
- [`ms-polygenic-risk-score`](https://github.com/HafsahShamsi/ms-polygenic-risk-score) — GWAS-based PRS pipeline
- `pituitary-adenoma-deg` — this repo
- `drug-target-analysis` — coming next
