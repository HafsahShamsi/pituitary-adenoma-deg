"""
Pathway Enrichment — Pituitary Adenoma DEGs
============================================
Author: Hafsah Shamsi
GitHub: github.com/HafsahShamsi
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import requests
import os

os.makedirs("./outputs", exist_ok=True)

GOLD  = "#C9A84C"
DARK  = "#1A1A2E"
MUTED = "#4A4A6A"
UP_COL   = "#E05C5C"
DOWN_COL = "#5C9BE0"

def enrichr_query(gene_list, gene_set, label):
    if not gene_list:
        return pd.DataFrame()
    r = requests.post(
        "https://maayanlab.cloud/Enrichr/addList",
        files={"list": (None, "\n".join(gene_list[:500])),
               "description": (None, label)}
    )
    lid = r.json()["userListId"]
    r2  = requests.get(
        f"https://maayanlab.cloud/Enrichr/enrich?userListId={lid}&backgroundType={gene_set}"
    )
    rows = r2.json().get(gene_set, [])
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows, columns=[
        "rank","term","pval","zscore","combined_score",
        "genes","adj_pval","old_pval","old_adj_pval"
    ])
    df["group"] = label
    return df

def plot_dotplot(up_df, down_df, gene_set_label):
    top_up   = up_df.nsmallest(10, "adj_pval").copy()
    top_down = down_df.nsmallest(10, "adj_pval").copy()

    top_up["direction"]   = "Upregulated"
    top_down["direction"] = "Downregulated"
    combined = pd.concat([top_up, top_down])

    combined["-log10p"] = -np.log10(combined["adj_pval"].clip(1e-30))
    combined["term_short"] = combined["term"].str[:55]

    fig, ax = plt.subplots(figsize=(10, 9))
    fig.patch.set_facecolor(DARK)
    ax.set_facecolor(DARK)

    colors = [UP_COL if d == "Upregulated" else DOWN_COL
              for d in combined["direction"]]

    scatter = ax.scatter(
        combined["-log10p"],
        range(len(combined)),
        c=colors,
        s=combined["-log10p"] * 12,
        alpha=0.85,
        edgecolors=GOLD,
        linewidths=0.4
    )

    ax.set_yticks(range(len(combined)))
    ax.set_yticklabels(combined["term_short"], fontsize=8, color="white")
    ax.set_xlabel("-log10 adjusted p-value", color=GOLD, fontsize=10)
    ax.set_title(f"Pathway Enrichment — Pituitary Adenoma\n{gene_set_label}",
                 color=GOLD, fontsize=11, fontweight="bold")
    ax.tick_params(colors="white")
    for spine in ax.spines.values():
        spine.set_edgecolor(MUTED)

    import matplotlib.patches as mpatches
    up_patch   = mpatches.Patch(color=UP_COL,   label="Upregulated DEGs")
    down_patch = mpatches.Patch(color=DOWN_COL, label="Downregulated DEGs")
    ax.legend(handles=[up_patch, down_patch],
              facecolor=DARK, edgecolor=GOLD, labelcolor="white", fontsize=9)

    ax.axhline(9.5, color=MUTED, linewidth=0.6, linestyle="--")

    plt.tight_layout()
    fname = f"./outputs/enrichment_{gene_set_label.replace(' ','_')}.png"
    plt.savefig(fname, dpi=180, bbox_inches="tight", facecolor=DARK)
    plt.close()
    print(f"  Saved -> {fname}")

if __name__ == "__main__":
    deg = pd.read_csv("./outputs/pituitary_deg_results.csv", index_col=0)

    up_genes   = deg[(deg.padj < 0.05) & (deg.log2FC > 1)].index.tolist()
    down_genes = deg[(deg.padj < 0.05) & (deg.log2FC < -1)].index.tolist()

    # Clean gene names — remove /// entries
    up_genes   = [g for g in up_genes   if "///" not in g]
    down_genes = [g for g in down_genes if "///" not in g]

    print(f"Upregulated genes for enrichment: {len(up_genes)}")
    print(f"Downregulated genes for enrichment: {len(down_genes)}")

    for gene_set in ["KEGG_2021_Human", "GO_Biological_Process_2021"]:
        print(f"\nQuerying {gene_set}...")
        up_df   = enrichr_query(up_genes,   gene_set, "upregulated")
        down_df = enrichr_query(down_genes, gene_set, "downregulated")

        if not up_df.empty:
            print(f"  Top upregulated pathways:")
            for _, r in up_df.nsmallest(5, "adj_pval").iterrows():
                print(f"    {r['term'][:60]}  adj_p={r['adj_pval']:.2e}")

        if not down_df.empty:
            print(f"  Top downregulated pathways:")
            for _, r in down_df.nsmallest(5, "adj_pval").iterrows():
                print(f"    {r['term'][:60]}  adj_p={r['adj_pval']:.2e}")

        if not up_df.empty and not down_df.empty:
            plot_dotplot(up_df, down_df, gene_set)

        up_df.to_csv(f"./outputs/enrichment_up_{gene_set}.csv", index=False)
        down_df.to_csv(f"./outputs/enrichment_down_{gene_set}.csv", index=False)

    print("\nDone.")