#!/usr/bin/env python3
"""
analyze_trees_v3.py
Advanced phylogenetic + orthology analysis for RNAmodMachinery v2
------------------------------------------------------------------
Extracts:
- Branch statistics, duplication rates, entropy, and IQ-TREE info
- Domain × Species orthology/expansion matrix from HMMER tblouts

Author: Oguzhan Begik
"""

import os
import re
import json
import math
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ete3 import Tree
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed


# ============================================================
#                    Helper Functions
# ============================================================

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def shannon_entropy(values):
    total = sum(values)
    probs = [v / total for v in values if v > 0]
    return -sum(p * math.log2(p) for p in probs)


def parse_species_from_name(name):
    """Extract species code from UniProt-like name."""
    if "_" in name:
        return name.split("_")[-1].replace("|", "")
    return "Unknown"


def parse_tblout_counts(tbl_files, evalue_cutoff=1e-5):
    """Parse all tblout files to get counts per species."""
    counts = defaultdict(lambda: defaultdict(int))
    for tbl in tbl_files:
        domain = os.path.basename(os.path.dirname(tbl))
        species = os.path.splitext(os.path.basename(tbl))[0].replace(".tblout", "")
        try:
            with open(tbl) as f:
                for line in f:
                    if line.startswith("#"):
                        continue
                    parts = line.strip().split()
                    if len(parts) > 4 and float(parts[4]) < evalue_cutoff:
                        counts[domain][species] += 1
        except Exception as e:
            print(f"⚠️  Could not parse {tbl}: {e}")
    return counts


def analyze_tree(treefile, iqtree_file=None):
    """Analyze a single .treefile and return a feature dictionary."""
    t = Tree(treefile, format=1)
    leaves = t.get_leaves()
    n_leaves = len(leaves)
    branch_lengths = [n.dist for n in t.traverse() if not n.is_root()]
    mean_branch = np.mean(branch_lengths) if branch_lengths else 0
    max_branch = np.max(branch_lengths) if branch_lengths else 0

    # Species distribution
    species = [parse_species_from_name(l.name) for l in leaves]
    species_counts = Counter(species)
    n_species = len(species_counts)
    entropy = shannon_entropy(list(species_counts.values()))

    # Duplication detection
    duplication_nodes = 0
    for node in t.traverse("postorder"):
        if not node.is_leaf():
            child_species = [parse_species_from_name(l.name) for l in node.get_leaves()]
            if len(child_species) > len(set(child_species)):
                duplication_nodes += 1
    duplication_rate = duplication_nodes / len([n for n in t.traverse() if not n.is_leaf()])

    # Parse IQ-TREE info safely
    model, logl, mean_bs = None, None, None
    if iqtree_file and os.path.exists(iqtree_file):
        with open(iqtree_file) as f:
            for line in f:
                if "Best-fit model" in line:
                    model = line.split(":")[-1].strip()
                elif "Log-likelihood of the tree:" in line:
                    val = re.sub(r"[^\d\.\-]", "", line.split(":")[-1])
                    try:
                        logl = float(val)
                    except ValueError:
                        pass
                elif "Average support:" in line:
                    val = re.sub(r"[^\d\.]", "", line.split(":")[-1])
                    try:
                        mean_bs = float(val)
                    except ValueError:
                        pass

    return {
        "domain": os.path.basename(treefile).replace("_aligned.fasta.treefile", ""),
        "n_leaves": n_leaves,
        "n_species": n_species,
        "mean_branch_len": round(mean_branch, 4),
        "max_branch_len": round(max_branch, 4),
        "duplication_rate": round(duplication_rate, 4),
        "entropy_species": round(entropy, 4),
        "model": model,
        "log_likelihood": logl,
        "mean_support": mean_bs
    }


def process_domain(domain_dir):
    treefiles = glob.glob(os.path.join(domain_dir, "*aligned.fasta.treefile"))
    results = []
    for tf in treefiles:
        iqtree_file = tf.replace(".treefile", ".iqtree")
        try:
            feat = analyze_tree(tf, iqtree_file)
            results.append(feat)
        except Exception as e:
            print(f"⚠️  Failed to process {tf}: {e}")
    return results


# ============================================================
#                    Main Runner
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="Analyze phylogenetic trees and build orthology matrix.")
    parser.add_argument("--tree-dir", required=True, help="Directory containing domain subfolders with .treefile outputs")
    parser.add_argument("--outdir", required=True, help="Output directory for results")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for parallel processing")
    parser.add_argument("--evalue", type=float, default=1e-5, help="E-value cutoff for tblout parsing")
    args = parser.parse_args()

    ensure_dir(args.outdir)
    ensure_dir(os.path.join(args.outdir, "logs"))

    domain_dirs = [os.path.join(args.tree_dir, d) for d in os.listdir(args.tree_dir)
                   if os.path.isdir(os.path.join(args.tree_dir, d))]
    print(f"🧩 Found {len(domain_dirs)} domain folders to analyze")

    # ----- 1. Analyze Trees -----
    all_features = []
    with ThreadPoolExecutor(max_workers=args.threads) as ex:
        futures = {ex.submit(process_domain, d): d for d in domain_dirs}
        for fut in as_completed(futures):
            res = fut.result()
            if res:
                all_features.extend(res)

    if not all_features:
        print("❌ No trees processed successfully.")
        return

    df = pd.DataFrame(all_features)
    df.to_csv(os.path.join(args.outdir, "tree_features.tsv"), sep="\t", index=False)

    # ----- 2. Parse tblouts for Orthology Matrix -----
    all_tblouts = glob.glob(os.path.join(args.tree_dir, "*/*.tblout"))
    print(f"📂 Parsing {len(all_tblouts)} tblout files for orthology counts...")
    counts = parse_tblout_counts(all_tblouts, evalue_cutoff=args.evalue)

    df_matrix = pd.DataFrame(counts).fillna(0).astype(int).T
    df_matrix.to_csv(os.path.join(args.outdir, "orthology_matrix.tsv"), sep="\t")
    print(f"✅ Orthology matrix saved: {args.outdir}/orthology_matrix.tsv")
    # ----- 3. Labeled Heatmap (Bulletproof Save) -----
    try:
        import matplotlib
        matplotlib.use("Agg")  # Force non-interactive backend
        import matplotlib.pyplot as plt
        import seaborn as sns

        print("\n🧬 Starting heatmap generation...")

        # Check output path
        output_path = os.path.join(args.outdir, "orthology_matrix_heatmap_labeled.png")
        os.makedirs(args.outdir, exist_ok=True)
        print(f"📁 Output directory confirmed: {args.outdir}")

        # Check matrix validity
        if df_matrix.empty:
            print("❌ Matrix is empty — skipping heatmap.")
        else:
            df_plot = np.log1p(df_matrix)  # log(1+x) transform for visibility

            plt.figure(figsize=(max(12, len(df_matrix.columns) * 0.4),
                                max(8, len(df_matrix.index) * 0.35)))
            ax = sns.heatmap(
                df_plot,
                cmap="viridis",
                linewidths=0.1,
                cbar_kws={"label": "log(Number of Homologs + 1)"},
                xticklabels=True,
                yticklabels=True
            )

            plt.title("Domain × Species Expansion Matrix", fontsize=14, pad=15)
            plt.xlabel("Species", fontsize=12, labelpad=10)
            plt.ylabel("Pfam Domain", fontsize=12, labelpad=10)
            plt.xticks(rotation=90, fontsize=8)
            plt.yticks(fontsize=9)
            plt.tight_layout()

            plt.savefig(output_path, dpi=300, bbox_inches="tight")
            plt.close()
            print(f"✅ Heatmap saved successfully: {output_path}\n")

    except Exception as e:
        import traceback
        print("❌ Heatmap generation failed with error:")
        traceback.print_exc()
        

if __name__ == "__main__":
    main()