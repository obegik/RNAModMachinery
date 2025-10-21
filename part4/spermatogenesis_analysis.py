#!/usr/bin/env python3
# ================================================================
# Spermatogenesis Expression Analysis (Green et al., 2019)
# Adapted from Oguzhan Begik (Genome Biology, 2020)
# ================================================================

import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path
import warnings
warnings.filterwarnings("ignore")

# ================================================================
# 1. Input arguments
# ================================================================

if len(sys.argv) < 3:
    print("Usage: python spermatogenesis_analysis.py <expression.tsv> <gene_map.tsv>")
    sys.exit(1)

expr_path = Path(sys.argv[1])
map_path = Path(sys.argv[2])
outdir = Path("results_spermatogenesis")
outdir.mkdir(exist_ok=True)

print(f"📥 Loading {expr_path.name} and {map_path.name} ...")

# ================================================================
# 2. Load and preprocess
# ================================================================

data = pd.read_csv(expr_path, sep="\t")
ensembl = pd.read_csv(map_path, sep="\t")

# Expect columns: Gene, ENSEMBL, HGNC
if "Gene" not in data.columns:
    data.rename(columns={data.columns[0]: "Gene"}, inplace=True)

joined = pd.merge(data, ensembl, on="Gene", how="left")
joined = joined.sort_values("ENSEMBL").dropna(subset=["ENSEMBL"])
joined = joined[["HGNC"] + [c for c in joined.columns if c not in ["HGNC"]]]

# ================================================================
# 3. Aggregate stages
# ================================================================

print("🧫 Aggregating spermatogenesis stages...")

clustered = pd.DataFrame({
    "HGNC": joined["HGNC"],
    "Spermatogonia": joined["GC1"],
    "Prelep_Spermatocyte": joined[["GC2", "GC3"]].mean(axis=1),
    "Spermatocytes": joined[["GC4", "GC5", "GC6", "GC7", "GC8"]].mean(axis=1),
    "Spermatids": joined[["GC9", "GC10", "GC11"]].mean(axis=1),
    "Elongating_Spermatids": joined["GC12"]
})

clustered.set_index("HGNC", inplace=True)
clustered = clustered.apply(pd.to_numeric)

# ================================================================
# 4. K-means clustering
# ================================================================

print("🔬 Performing k-means clustering (k=4)...")
kmeans = KMeans(n_clusters=4, n_init=25, random_state=42)
cluster_labels = kmeans.fit_predict(clustered)
clustered["Cluster"] = cluster_labels + 1  # 1–4

# ================================================================
# 5. PCA
# ================================================================

print("🧠 Performing PCA...")
scaler = StandardScaler()
X_scaled = scaler.fit_transform(clustered.drop(columns="Cluster"))
pca = PCA(n_components=2)
pca_data = pca.fit_transform(X_scaled)
pca_df = pd.DataFrame(pca_data, columns=["PC1", "PC2"], index=clustered.index)
pca_df["Cluster"] = clustered["Cluster"]

expl_var = np.round(pca.explained_variance_ratio_ * 100, 2)
pc_labels = [f"PC1 ({expl_var[0]}%)", f"PC2 ({expl_var[1]}%)"]

# --- PCA plots ---
with PdfPages(outdir / "pca_spermatogenesis_clustered.pdf") as pdf:
    plt.figure(figsize=(6, 5))
    sns.scatterplot(
        data=pca_df, x="PC1", y="PC2", hue="Cluster", palette="Set2", s=60, edgecolor="black"
    )
    plt.axhline(0, ls="--", color="gray")
    plt.axvline(0, ls="--", color="gray")
    plt.xlabel(pc_labels[0])
    plt.ylabel(pc_labels[1])
    plt.title("PCA of Spermatogenesis Expression")
    pdf.savefig(bbox_inches="tight")
    plt.close()

# ================================================================
# 6. PCA loadings (rotation)
# ================================================================

loadings = pd.DataFrame(
    pca.components_.T, index=clustered.columns[:-1], columns=["PC1", "PC2"]
)
with PdfPages(outdir / "pca_spermatogenesis_loadings.pdf") as pdf:
    plt.figure(figsize=(7, 5))
    sns.scatterplot(x="PC1", y="PC2", data=loadings, s=100)
    for i, txt in enumerate(loadings.index):
        plt.text(loadings.PC1[i], loadings.PC2[i], txt, fontsize=9)
    plt.axhline(0, ls="--", color="gray")
    plt.axvline(0, ls="--", color="gray")
    plt.title("PCA Stage Loadings")
    pdf.savefig(bbox_inches="tight")
    plt.close()

# ================================================================
# 7. Heatmap
# ================================================================

print("🔥 Generating heatmap...")
heatmap_data = clustered.drop(columns="Cluster").copy()
heatmap_data = np.log1p(heatmap_data)

plt.figure(figsize=(8, 12))
sns.clustermap(
    heatmap_data,
    cmap=sns.color_palette("coolwarm", as_cmap=True),
    row_cluster=True,
    col_cluster=False,
    yticklabels=False,
    figsize=(8, 12)
)
plt.savefig(outdir / "spermatogenesis_heatmap.png", dpi=300, bbox_inches="tight")
plt.close()

# ================================================================
# 8. Violin plots
# ================================================================

print("🎻 Plotting violin plots...")
long_df = (
    clustered.drop(columns="Cluster")
    .reset_index()
    .melt(id_vars="HGNC", var_name="Stage", value_name="Expression")
)
long_df["Cluster"] = np.repeat(clustered["Cluster"].values, len(clustered.columns) - 1)

plt.figure(figsize=(10, 5))
sns.violinplot(data=long_df, x="Stage", y="Expression", hue="Cluster", split=True, palette="Set2")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig(outdir / "spermatogenesis_violinplot.png", dpi=300)
plt.close()

# ================================================================
# 9. Barplots
# ================================================================

print("📊 Generating per-gene barplots...")
bar_long = (
    clustered.reset_index()
    .melt(id_vars=["HGNC", "Cluster"], var_name="Stage", value_name="Expression")
)
g = sns.FacetGrid(bar_long, col="HGNC", col_wrap=6, sharey=False, height=2)
g.map_dataframe(sns.barplot, x="Stage", y="Expression", hue="Stage", dodge=False)
g.set_titles("{col_name}")
g.set_xticklabels(rotation=45, ha="right")
plt.tight_layout()
g.savefig(outdir / "spermatogenesis_barplots.png", dpi=200)
plt.close()

# ================================================================
# 10. Export tables
# ================================================================

clustered.to_csv(outdir / "spermatogenesis_clustered.tsv", sep="\t")
pca_df.to_csv(outdir / "spermatogenesis_pca.tsv", sep="\t")
print(f"✅ Analysis complete! Results in {outdir.resolve()}")