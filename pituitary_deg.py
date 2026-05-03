"""
Pituitary Adenoma DEG Analysis
================================
Author: Hafsah Shamsi
GitHub: github.com/HafsahShamsi

Differential gene expression analysis of pituitary adenoma vs normal pituitary tissue.
Dataset: GSE26966 — 23 samples, GPL570 (Affymetrix HG-U133 Plus 2.0)

Pituitary adenomas lose negative feedback sensitivity — cells keep secreting
hormones regardless of circulating levels. This script identifies the transcriptional
signature of that feedback loss.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import GEOparse
import os
import warnings
warnings.filterwarnings("ignore")

from scipy import stats
from statsmodels.stats.multitest import multipletests

os.makedirs("./outputs", exist_ok=True)
os.makedirs("./geo_cache", exist_ok=True)

GOLD       = "#C9A84C"
DARK       = "#1A1A2E"
MUTED      = "#4A4A6A"
UP_COL     = "#E05C5C"
DOWN_COL   = "#5C9BE0"
NEUT_COL   = "#4A4A6A"

# ── FETCH + DEG PIPELINE ─────────────────────────────────────────────────────

def run_deg_pipeline(gse_id="GSE26966"):
    print(f"Fetching {gse_id}...")
    gse = GEOparse.get_GEO(geo=gse_id, destdir="./geo_cache", silent=True)
    gpl = gse.gpls[list(gse.gpls.keys())[0]]

    # Find gene symbol column
    sym_col = None
    for c in ["Gene Symbol", "SYMBOL", "gene_symbol"]:
        if c in gpl.table.columns:
            sym_col = c
            break
    if sym_col is None:
        for c in gpl.table.columns:
            if "symbol" in c.lower():
                sym_col = c
                break
    print(f"  Symbol column: {sym_col}")

    annot = gpl.table[["ID", sym_col]].dropna()
    annot.columns = ["probe_id", "gene"]
    annot = annot[annot.gene.str.strip() != ""].set_index("probe_id")

    # Build expression matrix
    tables = []
    for n, gsm in gse.gsms.items():
        try:
            t = gsm.table[["ID_REF", "VALUE"]].rename(
                columns={"ID_REF": "probe_id", "VALUE": n}
            ).set_index("probe_id")
            tables.append(t)
        except Exception:
            continue

    expr = pd.concat(tables, axis=1).apply(pd.to_numeric, errors="coerce").dropna()
    expr = expr.join(annot).groupby("gene").mean()
    print(f"  Expression matrix: {expr.shape}")

    # Sample classification
    tumor_s, normal_s = [], []
    for n, gsm in gse.gsms.items():
        title = gsm.metadata.get("title", [""])[0].lower()
        char  = " ".join(gsm.metadata.get("characteristics_ch1", [])).lower()
        combined = title + " " + char
        if "normal" in combined:
            normal_s.append(n)
        else:
            tumor_s.append(n)

    tumor_s  = [s for s in tumor_s  if s in expr.columns]
    normal_s = [s for s in normal_s if s in expr.columns]
    print(f"  Tumor={len(tumor_s)}, Normal={len(normal_s)}")

    # DEG computation
    fc  = expr[tumor_s].mean(axis=1) - expr[normal_s].mean(axis=1)
    pv  = expr.apply(
        lambda r: stats.ttest_ind(r[tumor_s].values, r[normal_s].values).pvalue, axis=1
    )
    _, padj, _, _ = multipletests(pv.fillna(1), method="fdr_bh")
    deg = pd.DataFrame({"log2FC": fc, "padj": padj}, index=expr.index)

    up   = deg[(deg.padj < 0.05) & (deg.log2FC > 1)]
    down = deg[(deg.padj < 0.05) & (deg.log2FC < -1)]
    print(f"  DEGs: {len(up)} up, {len(down)} down")

    deg.to_csv("./outputs/pituitary_deg_results.csv")
    return deg, up, down

# ── VOLCANO PLOT ──────────────────────────────────────────────────────────────

def plot_volcano(deg, up, down):
    fig, ax = plt.subplots(figsize=(9, 6))
    fig.patch.set_facecolor(DARK)
    ax.set_facecolor(DARK)

    neg_log_p = -np.log10(deg.padj.clip(1e-300))

    # Neutral
    mask_neut = ~deg.index.isin(up.index) & ~deg.index.isin(down.index)
    ax.scatter(deg.loc[mask_neut, "log2FC"], neg_log_p[mask_neut],
               c=NEUT_COL, s=6, alpha=0.4, linewidths=0)

    # Down
    ax.scatter(deg.loc[down.index, "log2FC"], neg_log_p[down.index],
               c=DOWN_COL, s=8, alpha=0.7, linewidths=0)

    # Up
    ax.scatter(deg.loc[up.index, "log2FC"], neg_log_p[up.index],
               c=UP_COL, s=8, alpha=0.7, linewidths=0)

    # Label top 10 by significance
    top = deg[deg.index.isin(up.index) | deg.index.isin(down.index)]
    top = top.nsmallest(10, "padj")
    for gene, row in top.iterrows():
        ax.annotate(gene,
                    xy=(row.log2FC, -np.log10(max(row.padj, 1e-300))),
                    fontsize=7, color="white", alpha=0.9,
                    xytext=(4, 2), textcoords="offset points")

    ax.axvline(1,  color=GOLD, linestyle="--", linewidth=0.8, alpha=0.6)
    ax.axvline(-1, color=GOLD, linestyle="--", linewidth=0.8, alpha=0.6)
    ax.axhline(-np.log10(0.05), color=GOLD, linestyle=":", linewidth=0.8, alpha=0.6)

    ax.set_xlabel("log2 Fold Change (Tumor vs Normal)", color=GOLD, fontsize=10)
    ax.set_ylabel("-log10 adjusted p-value", color=GOLD, fontsize=10)
    ax.set_title("Pituitary Adenoma vs Normal Pituitary\nDifferential Gene Expression (GSE26966)",
                 color=GOLD, fontsize=11, fontweight="bold")
    ax.tick_params(colors="white")
    for spine in ax.spines.values():
        spine.set_edgecolor(MUTED)

    up_patch   = mpatches.Patch(color=UP_COL,   label=f"Upregulated ({len(up)})")
    down_patch = mpatches.Patch(color=DOWN_COL, label=f"Downregulated ({len(down)})")
    ax.legend(handles=[up_patch, down_patch],
              facecolor=DARK, edgecolor=GOLD, labelcolor="white", fontsize=9)

    plt.tight_layout()
    plt.savefig("./outputs/volcano_pituitary.png", dpi=180, bbox_inches="tight", facecolor=DARK)
    plt.close()
    print("  Volcano saved -> ./outputs/volcano_pituitary.png")

# ── PCA ───────────────────────────────────────────────────────────────────────

def plot_pca(deg, gse_id="GSE26966"):
    from sklearn.preprocessing import StandardScaler
    from sklearn.decomposition import PCA

    gse = GEOparse.get_GEO(geo=gse_id, destdir="./geo_cache", silent=True)
    gpl = gse.gpls[list(gse.gpls.keys())[0]]

    sym_col = None
    for c in gpl.table.columns:
        if "symbol" in c.lower():
            sym_col = c
            break

    annot = gpl.table[["ID", sym_col]].dropna()
    annot.columns = ["probe_id", "gene"]
    annot = annot[annot.gene.str.strip() != ""].set_index("probe_id")

    tables = []
    labels = []
    for n, gsm in gse.gsms.items():
        try:
            t = gsm.table[["ID_REF", "VALUE"]].rename(
                columns={"ID_REF": "probe_id", "VALUE": n}
            ).set_index("probe_id")
            tables.append(t)
            title = gsm.metadata.get("title", [""])[0].lower()
            char  = " ".join(gsm.metadata.get("characteristics_ch1", [])).lower()
            labels.append("Normal" if "normal" in title + char else "Tumor")
        except Exception:
            continue

    expr = pd.concat(tables, axis=1).apply(pd.to_numeric, errors="coerce").dropna()
    expr = expr.join(annot).groupby("gene").mean()

    # Use top 500 most variable genes
    top500 = expr.var(axis=1).nlargest(500).index
    X = expr.loc[top500].T.values
    X = StandardScaler().fit_transform(X)

    pca = PCA(n_components=2)
    coords = pca.fit_transform(X)

    fig, ax = plt.subplots(figsize=(7, 5))
    fig.patch.set_facecolor(DARK)
    ax.set_facecolor(DARK)

    colors = [UP_COL if l == "Tumor" else DOWN_COL for l in labels]
    ax.scatter(coords[:, 0], coords[:, 1], c=colors, s=60, alpha=0.85, edgecolors=GOLD, linewidths=0.4)

    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)", color=GOLD, fontsize=10)
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)", color=GOLD, fontsize=10)
    ax.set_title("PCA — Pituitary Adenoma vs Normal\nTop 500 variable genes",
                 color=GOLD, fontsize=11, fontweight="bold")
    ax.tick_params(colors="white")
    for spine in ax.spines.values():
        spine.set_edgecolor(MUTED)

    tumor_patch  = mpatches.Patch(color=UP_COL,   label="Tumor")
    normal_patch = mpatches.Patch(color=DOWN_COL, label="Normal")
    ax.legend(handles=[tumor_patch, normal_patch],
              facecolor=DARK, edgecolor=GOLD, labelcolor="white", fontsize=9)

    plt.tight_layout()
    plt.savefig("./outputs/pca_pituitary.png", dpi=180, bbox_inches="tight", facecolor=DARK)
    plt.close()
    print("  PCA saved -> ./outputs/pca_pituitary.png")

# ── FEEDBACK GENE PLOT ────────────────────────────────────────────────────────

def plot_feedback_genes(deg):
    # Key HPA/feedback axis genes
    feedback_genes = [
        "NR3C1",   # glucocorticoid receptor
        "SSTR2",   # somatostatin receptor 2
        "SSTR5",   # somatostatin receptor 5
        "POMC",    # proopiomelanocortin
        "FKBP5",   # cortisol feedback marker
        "CRH",     # corticotropin releasing hormone
        "ACSL1",   # fatty acid metabolism
        "GHR",     # growth hormone receptor
        "PRLR",    # prolactin receptor
        "CDKN2A",  # cell cycle — p16
        "MKI67",   # proliferation marker
        "PCNA",    # proliferation
    ]

    present = [g for g in feedback_genes if g in deg.index]
    subset  = deg.loc[present].copy()
    subset  = subset.sort_values("log2FC")

    colors = [UP_COL if fc > 0 else DOWN_COL for fc in subset.log2FC]

    fig, ax = plt.subplots(figsize=(8, 6))
    fig.patch.set_facecolor(DARK)
    ax.set_facecolor(DARK)

    bars = ax.barh(range(len(subset)), subset.log2FC.values,
                   color=colors, edgecolor=GOLD, linewidth=0.5)
    ax.set_yticks(range(len(subset)))
    ax.set_yticklabels(subset.index.values, fontsize=9, color="white")
    ax.axvline(0, color=GOLD, linewidth=0.8)
    ax.set_xlabel("log2 Fold Change (Tumor vs Normal)", color=GOLD, fontsize=10)
    ax.set_title("HPA Axis & Feedback Gene Expression\nin Pituitary Adenoma",
                 color=GOLD, fontsize=11, fontweight="bold")
    ax.tick_params(colors="white")
    for spine in ax.spines.values():
        spine.set_edgecolor(MUTED)

    plt.tight_layout()
    plt.savefig("./outputs/feedback_genes.png", dpi=180, bbox_inches="tight", facecolor=DARK)
    plt.close()
    print("  Feedback genes saved -> ./outputs/feedback_genes.png")

# ── MAIN ──────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    print("="*55)
    print("Pituitary Adenoma DEG Analysis — GSE26966")
    print("="*55)

    deg, up, down = run_deg_pipeline()

    print("\nGenerating plots...")
    plot_volcano(deg, up, down)
    plot_pca(deg)
    plot_feedback_genes(deg)

    print("\nTop 10 upregulated:")
    print(up.nsmallest(10, "padj")[["log2FC","padj"]].to_string())
    print("\nTop 10 downregulated:")
    print(down.nsmallest(10, "padj")[["log2FC","padj"]].to_string())

    print(f"\nDone. Outputs -> ./outputs/")