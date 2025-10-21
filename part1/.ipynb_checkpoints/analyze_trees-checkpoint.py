#!/usr/bin/env python3
"""
analyze_trees.py
Summarize domain-based phylogenetic results from RNAModMachinery v2:
Extract species presence, duplication counts, and orthology patterns.
Author: Oguzhan Begik (RNAModMachinery v2)
"""

import os, re, glob, pandas as pd
from collections import defaultdict
from Bio import Phylo

# ===============================================================
#                 Utility Functions
# ===============================================================

def parse_species_from_id(seqname):
    """Extract species name from UniProt-like ID, e.g. sp|Q9HC16|ABC3G_HUMAN → HUMAN"""
    m = re.search(r'_([A-Z0-9]+)$', seqname)
    if m:
        return m.group(1)
    else:
        # fallback for names like MyGene_Species
        parts = seqname.split('_')
        return parts[-1] if len(parts) > 1 else "Unknown"

def summarize_tree(treefile):
    """Extract per-species hit counts from a Newick tree"""
    species_counts = defaultdict(int)
    tree = Phylo.read(treefile, "newick")

    for leaf in tree.get_terminals():
        sp = parse_species_from_id(leaf.name)
        species_counts[sp] += 1

    return species_counts

# ===============================================================
#                 Main Analysis
# ===============================================================

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Summarize phylogenetic results (presence, paralogs) across Pfam domains.")
    parser.add_argument("--indir", required=True, help="Directory containing Pfam subfolders (e.g. results_part1_v2/)")
    parser.add_argument("--outdir", default="tree_analysis", help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    treefiles = sorted(glob.glob(os.path.join(args.indir, "*/*.treefile")))
    if not treefiles:
        print("❌ No .treefile files found — run the discovery pipeline first.")
        return

    print(f"🌳 Found {len(treefiles)} treefiles to analyze")

    # Collect data
    presence = defaultdict(dict)
    paralogs = defaultdict(dict)
    all_species = set()

    for tf in treefiles:
        domain = os.path.basename(os.path.dirname(tf))
        counts = summarize_tree(tf)
        for sp, n in counts.items():
            presence[sp][domain] = 1
            paralogs[sp][domain] = n
            all_species.add(sp)
        # Mark absences later

    # Fill missing species
    all_species = sorted(all_species)
    all_domains = sorted({os.path.basename(os.path.dirname(f)) for f in treefiles})

    for sp in all_species:
        for d in all_domains:
            presence[sp].setdefault(d, 0)
            paralogs[sp].setdefault(d, 0)

    # Convert to DataFrames
    df_presence = pd.DataFrame.from_dict(presence, orient="index")[all_domains]
    df_paralogs = pd.DataFrame.from_dict(paralogs, orient="index")[all_domains]

    # Save outputs
    df_presence.to_csv(os.path.join(args.outdir, "species_presence.tsv"), sep="\t")
    df_paralogs.to_csv(os.path.join(args.outdir, "paralog_summary.tsv"), sep="\t")

    print(f"✅ Saved species presence matrix → {args.outdir}/species_presence.tsv")
    print(f"✅ Saved paralog counts matrix → {args.outdir}/paralog_summary.tsv")

    # Optional summary stats
    summary = df_paralogs.astype(bool).sum().reset_index()
    summary.columns = ["Pfam_Domain", "Species_with_hits"]
    summary["Total_species"] = len(all_species)
    summary["Coverage_%"] = (summary["Species_with_hits"] / summary["Total_species"] * 100).round(2)
    summary.to_csv(os.path.join(args.outdir, "domain_coverage_summary.tsv"), sep="\t", index=False)

    print(f"📊 Domain-level summary written → {args.outdir}/domain_coverage_summary.tsv")

if __name__ == "__main__":
    main()