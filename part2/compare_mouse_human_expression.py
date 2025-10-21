#!/usr/bin/env python3
# ============================================================
# 🧬 GTEx (Human) vs ENCODE (Mouse) Expression Correlation
# Author: Oguzhan Begik | 2025
# Purpose: Replicate gtex_manipulation.R, encode_manipulation.R,
#          gtex_tissuewide.R, and gtex_vs_encode_similarity.R
# ============================================================

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path

# -------------------------------
# === Configuration ===
# -------------------------------
GTEX_PATH = "data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.gct.tsv"
ENCODE_PATH = "data/encode_TPM_renamed.tsv"
ANNOT_PATH = "data/mouse_id_symbol_class.tsv"
OUTDIR = Path("results_part2")
OUTDIR.mkdir(parents=True, exist_ok=True)

# -------------------------------
# === 1. Load and preprocess ===
# -------------------------------

def load_gtex():
    """Load GTEx gene-level TPM table and clean columns."""
    print("📥 Loading GTEx (human) dataset...")
    df = pd.read_csv(GTEX_PATH, sep="\t", comment="#")
    if "Description" in df.columns:
        df = df.rename(columns={"Name": "gene_id", "Description": "symbol"})
    else:
        df = df.rename(columns={"gene_id": "gene_id", "gene_name": "symbol"})
    # Drop duplicates and clean
    df = df.drop_duplicates(subset="gene_id")
    print(f"✅ GTEx shape: {df.shape}")
    return df


def load_encode():
    """Load ENCODE TPM (mouse) dataset."""
    print("📥 Loading ENCODE (mouse) dataset...")
    df = pd.read_csv(ENCODE_PATH, sep="\t")
    df.columns = df.columns.str.strip()
    # detect gene/symbol columns
    possible_cols = [c for c in df.columns if "gene" in c.lower() or "symbol" in c.lower()]
    gene_col = possible_cols[0]
    df = df.rename(columns={gene_col: "symbol"})
    print(f"✅ ENCODE shape: {df.shape}")
    return df


def load_annotations():
    """Load RNAmod annotation linking symbols to classes."""
    print("📥 Loading RNAmod annotations...")
    df = pd.read_csv(ANNOT_PATH, sep="\t")
    df.columns = df.columns.str.strip()
    df = df.rename(columns={"gene_id": "gene_id", "Symbol": "symbol", "Class": "class"})
    print(f"✅ Annotation shape: {df.shape}")
    return df

# -------------------------------
# === 2. Normalization ===
# -------------------------------

def normalize_tpm(df, log_transform=True):
    """Log-transform and z-score normalize."""
    df_expr = df.drop(columns=["gene_id", "symbol"], errors="ignore")
    df_expr = df_expr.apply(pd.to_numeric, errors="coerce").fillna(0)
    if log_transform:
        df_expr = np.log2(df_expr + 1)
    # z-score normalization per gene
    df_expr = (df_expr - df_expr.mean(axis=1).values.reshape(-1,1)) / df_expr.std(axis=1).values.reshape(-1,1)
    df_norm = pd.concat([df[["gene_id", "symbol"]], df_expr], axis=1)
    return df_norm

# -------------------------------
# === 3. Average per tissue ===
# -------------------------------

def summarize_tissues(df):
    """Compute average TPM per tissue (handles replicates)."""
    tissue_cols = [c for c in df.columns if c not in ["gene_id", "symbol"]]
    grouped = df[tissue_cols].groupby(level=0, axis=1).mean() if any("." in c for c in tissue_cols) else df[tissue_cols]
    df_sum = pd.concat([df[["gene_id", "symbol"]], grouped], axis=1)
    return df_sum

# -------------------------------
# === 4. Correlation Analysis ===
# -------------------------------

def compute_tissue_correlation(gtex, encode, class_name, out_prefix):
    """Compute Pearson correlation between tissues across matched genes."""
    common_genes = set(gtex["symbol"]) & set(encode["symbol"])
    if len(common_genes) < 3:
        print(f"⚠️ Skipping {class_name}: only {len(common_genes)} overlapping genes.")
        return None

    gtex_sub = gtex[gtex["symbol"].isin(common_genes)].set_index("symbol")
    encode_sub = encode[encode["symbol"].isin(common_genes)].set_index("symbol")

    gtex_vals = gtex_sub.drop(columns=["gene_id"], errors="ignore")
    encode_vals = encode_sub.drop(columns=["gene_id"], errors="ignore")

    corr = np.corrcoef(gtex_vals.T, encode_vals.T)[: len(gtex_vals.columns), len(gtex_vals.columns):]
    corr_df = pd.DataFrame(
        corr,
        index=gtex_vals.columns,
        columns=encode_vals.columns
    )

    out_csv = OUTDIR / f"{out_prefix}_correlation.csv"
    corr_df.to_csv(out_csv)
    print(f"💾 Saved correlation table → {out_csv}")

    # Plot heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(corr_df, cmap="coolwarm", center=0, linewidths=0.5)
    plt.title(f"Pearson correlation: {class_name}")
    plt.tight_layout()
    out_png = OUTDIR / f"{out_prefix}_heatmap.png"
    plt.savefig(out_png, dpi=300)
    plt.close()
    print(f"🧭 Saved heatmap → {out_png}")

    return corr_df

# -------------------------------
# === 5. Main workflow ===
# -------------------------------

def main():
    gtex = load_gtex()
    encode = load_encode()
    annot = load_annotations()

    # Normalize
    gtex_norm = normalize_tpm(gtex)
    encode_norm = normalize_tpm(encode)

    # Summarize per tissue
    gtex_sum = summarize_tissues(gtex_norm)
    encode_sum = summarize_tissues(encode_norm)

    # Merge annotations to ENCODE genes
    encode_annot = encode_sum.merge(annot, on="symbol", how="left")

    # Unique RNAmod classes
    classes = encode_annot["class"].dropna().unique()
    print(f"📊 Found {len(classes)} RNAmod classes: {list(classes)}")

    summary = []
    for cls in classes:
        sub_encode = encode_annot[encode_annot["class"] == cls]
        print(f"\n🧩 Processing class: {cls} ({len(sub_encode)} genes)")
        corr_df = compute_tissue_correlation(gtex_sum, sub_encode, cls, f"tissue_corr_{cls}")
        if corr_df is not None:
            mean_corr = corr_df.stack().mean()
            summary.append({"Class": cls, "MeanCorrelation": mean_corr})

    if summary:
        summary_df = pd.DataFrame(summary)
        summary_path = OUTDIR / "tissue_correlation_summary.csv"
        summary_df.to_csv(summary_path, index=False)
        print(f"\n✅ Summary saved → {summary_path}")
        print(summary_df)

if __name__ == "__main__":
    main()