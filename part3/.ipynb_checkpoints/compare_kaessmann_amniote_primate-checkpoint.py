#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
compare_kaessmann_amniote_primate.py
====================================
Python adaptation of kaessmann.amniote.R and kaessmann.primate.R
Faithfully reproduces the R logic for RNA modification enzyme expression correlations
across amniote and primate orthologues.
"""

import os
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
from matplotlib.colors import LinearSegmentedColormap
import warnings
warnings.filterwarnings("ignore")

# === Paths ===
AMNIOTE_FILE = "NormalizedRPKM_ConstitutiveExons_Amniote1to1Orthologues.txt"
PRIMATE_FILE = "NormalizedRPKM_ConstitutiveExons_Primate1to1Orthologues.txt"
RNAMOD_FILE = "human_id_symbol_class.tsv"

OUTDIR = "results_kaessmann"
os.makedirs(OUTDIR, exist_ok=True)

# === 1. Load datasets ===
print("📥 Loading Kaessmann datasets...")

amniote = pd.read_csv(AMNIOTE_FILE, sep="\t")
primate = pd.read_csv(PRIMATE_FILE, sep="\t")
rnamod = pd.read_csv(RNAMOD_FILE, sep="\t")

rnamod.columns = [c.strip() for c in rnamod.columns]
print(f"✅ Amniote shape: {amniote.shape}")
print(f"✅ Primate shape: {primate.shape}")
print(f"✅ RNAmod annotation shape: {rnamod.shape}")

# === 2. Clean and normalize ===
def clean_expression(df):
    # Expecting first column as gene_id or symbol
    gene_col = df.columns[0]
    df = df.dropna(subset=[gene_col])
    df = df.set_index(gene_col)
    # Convert to numeric
    df = df.apply(pd.to_numeric, errors='coerce')
    # Replace 0s with small pseudocount
    df[df == 0] = np.nan
    df = np.log2(df + 1)
    return df

amniote_expr = clean_expression(amniote)
primate_expr = clean_expression(primate)

# === 3. Merge class annotation ===
# Harmonize IDs: handle both Ensembl IDs and gene symbols
amniote_expr.index = amniote_expr.index.str.split(".").str[0].str.upper().str.strip()
primate_expr.index = primate_expr.index.str.split(".").str[0].str.upper().str.strip()
rnamod["gene_id"] = rnamod["gene_id"].str.split(".").str[0].str.upper().str.strip()
rnamod["Symbol"] = rnamod["Symbol"].str.upper().str.strip()


# === Harmonize Ensembl IDs ===
def clean_ensembl_ids(x):
    return x.str.replace(r"\..*", "", regex=True).str.upper().str.strip()

amniote_expr.index = clean_ensembl_ids(pd.Series(amniote_expr.index))
primate_expr.index = clean_ensembl_ids(pd.Series(primate_expr.index))
rnamod["gene_id"] = clean_ensembl_ids(rnamod["gene_id"])
rnamod["Symbol"] = rnamod["Symbol"].str.upper().str.strip()

# Rebuild mapping
gene_class_map = dict(zip(rnamod["gene_id"], rnamod["Class"]))

# Determine shared genes
shared_genes_amniote = [g for g in amniote_expr.index if g in gene_class_map]
shared_genes_primate = [g for g in primate_expr.index if g in gene_class_map]
common_genes = sorted(set(shared_genes_amniote) & set(shared_genes_primate))

print(f"🧬 Common annotated genes across datasets: {len(common_genes)}")
if len(common_genes) == 0:
    print("⚠️ No matches found. Double-check Ensembl IDs in Kaessmann and RNAmod files.")
    

# Build mapping from Ensembl IDs to RNAmod class
gene_class_map = dict(zip(rnamod["gene_id"], rnamod["Class"]))

# Filter by Ensembl IDs shared across datasets
shared_genes_amniote = [g for g in amniote_expr.index if g in gene_class_map]
shared_genes_primate = [g for g in primate_expr.index if g in gene_class_map]
common_genes = sorted(set(shared_genes_amniote) & set(shared_genes_primate))

print(f"🧬 Common annotated genes across datasets: {len(common_genes)}")
if len(common_genes) < 5:
    print("⚠️ Very few or no overlapping genes — check Ensembl ID versions.")
print(f"🧬 Common annotated genes across datasets: {len(common_genes)}")

amniote_expr = amniote_expr.loc[common_genes]
primate_expr = primate_expr.loc[common_genes]
classes = [gene_class_map[g] for g in common_genes]

# === 4. Compute tissue-wise correlations ===
print("🔬 Calculating tissue correlation matrices...")
tissues_amniote = amniote_expr.columns
tissues_primate = primate_expr.columns
shared_tissues = sorted(set(tissues_amniote) & set(tissues_primate))

print(f"🧠 Shared tissues between datasets: {len(shared_tissues)}")
if len(shared_tissues) == 0:
    print("⚠️ No shared tissue names — please harmonize column labels!")
else:
    print(f"🔍 Matched tissues: {shared_tissues}")

def compute_correlations(df1, df2):
    """
    Compute gene-wise correlation between two expression matrices (same genes).
    """
    corrmat = np.zeros((df1.shape[1], df2.shape[1]))
    for i, col1 in enumerate(df1.columns):
        for j, col2 in enumerate(df2.columns):
            valid = df1[col1].notna() & df2[col2].notna()
            if valid.sum() > 5:
                corrmat[i, j] = pearsonr(df1[col1][valid], df2[col2][valid])[0]
            else:
                corrmat[i, j] = np.nan
    corr_df = pd.DataFrame(corrmat, index=df1.columns, columns=df2.columns)
    return corr_df

# === 5. Compute per-class mean correlations ===
summary = []
for cls in sorted(set(classes)):
    cls_genes = [g for g, c in zip(common_genes, classes) if c == cls]
    if len(cls_genes) < 3:
        continue
    a_df = amniote_expr.loc[cls_genes]
    p_df = primate_expr.loc[cls_genes]
    corr = compute_correlations(a_df, p_df)
    mean_corr = np.nanmean(corr.values)
    summary.append({"Class": cls, "N_genes": len(cls_genes), "MeanCorr": mean_corr})

summary_df = pd.DataFrame(summary)
summary_df.to_csv(os.path.join(OUTDIR, "kaessmann_class_correlations.csv"), index=False)
print("✅ Class-level summary saved!")

# === 6. Global correlation matrix ===
print("🌍 Building global correlation matrix...")
global_corr = compute_correlations(amniote_expr, primate_expr)
global_corr.to_csv(os.path.join(OUTDIR, "kaessmann_global_corrmatrix.csv"))

# === 7. Plotting ===
print("🎨 Generating clustered heatmaps...")

# Custom blue-white-red gradient
cmap = LinearSegmentedColormap.from_list("bwr_custom", ["#4575b4", "white", "#d73027"])

def plot_heatmap(corr_df, title, outname):
    import seaborn as sns
    import numpy as np
    import matplotlib.pyplot as plt
    import os

    # === Sanitize correlation matrix ===
    corr_df = corr_df.replace([np.inf, -np.inf], np.nan)
    corr_df = corr_df.fillna(0)  # Replace remaining NaNs
    corr_df = corr_df.astype(float)

    # Check for any remaining invalid entries
    if not np.isfinite(corr_df.to_numpy()).all():
        print("⚠️ Non-finite values still found, replacing with zeros.")
        corr_df = np.nan_to_num(corr_df, nan=0.0, posinf=0.0, neginf=0.0)

    # === Plot ===
    plt.figure(figsize=(10, 8))
    sns.clustermap(
        corr_df,
        cmap="coolwarm",
        vmin=-1,
        vmax=1,
        linewidths=0.2,
        figsize=(10, 10),
        cbar_pos=(0.02, 0.8, 0.03, 0.15),
        method="average",
        metric="euclidean"
    )
    plt.suptitle(title, fontsize=14, y=1.02)
    plt.savefig(os.path.join(OUTDIR, outname), dpi=300, bbox_inches="tight")
    plt.close()

# Plot global
plot_heatmap(global_corr, "Amniote vs Primate (All Genes)", "kaessmann_global_heatmap.png")

# Plot per class
for cls in sorted(set(classes)):
    cls_genes = [g for g, c in zip(common_genes, classes) if c == cls]
    if len(cls_genes) < 3:
        continue
    corr_cls = compute_correlations(amniote_expr.loc[cls_genes], primate_expr.loc[cls_genes])
    fname = f"kaessmann_{cls}_heatmap.png"
    plot_heatmap(corr_cls, f"Amniote vs Primate – {cls}", fname)

print("✅ All plots and tables saved in results_kaessmann/")