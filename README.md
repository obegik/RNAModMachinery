#  RNA Modification Machinery Analysis (Python Adaptation)

**This repository is a Python-based adaptation** of the original analysis scripts from:

> **Begik O**, Lucas MC, Liu H, Ramirez JM, Mattick JS, and Novoa EM.  
> *Integrative analyses of the RNA modification machinery reveal tissue- and cancer-specific signatures.*  
> **Genome Biology**, May 2020.  
> [🔗 Paper link](https://rdcu.be/b3Z5I) | [📄 DOI: 10.1186/s13059-020-02009-z](https://doi.org/10.1186/s13059-020-02009-z)

---

##  Data Availability

All accompanying datasets can be accessed from the public CRG server:  
 [https://public-docs.crg.es/enovoa/public/OguzhanBegik/begik_2020_RMP](https://public-docs.crg.es/enovoa/public/OguzhanBegik/begik_2020_RMP)

---

## Overview

This Python fork reproduces and extends the full *RNA Modification Machinery* (RMM) analysis pipeline originally written in R.  
The new implementation leverages **pandas**, **seaborn**, and **scikit-learn** to perform robust and reproducible analyses of RNA modification enzyme expression across:

-  Human and mouse tissues (GTEx, ENCODE)
-  Amniote and Primate orthologs
-  Tumor vs normal tissues (TCGA/GTEx)
-  Single-cell spermatogenesis datasets
-  Cross-species expression correlations

All analysis parts correspond to those in the original publication, modernized in Python.

---

##  PART 1 — RNA Modification Enzyme Discovery and Phylogenetic Analysis

###  Overview
This module identifies and characterizes RNA modification enzymes (writers, readers, erasers) across multiple species using Pfam HMM domain profiles.  
It replaces the legacy shell scripts (`find_homologs.sh`, `mafft.sh`, `iqtree.sh`) with a unified Python workflow.

##  1. Run HMM-based discovery
Search each proteome for proteins containing RNA modification domains.

```bash
python part1/rnamod_discover.py --pfam_dir pfam_list --proteome_dir proteomes --out_dir results_part1
```
##  2. Multiple Sequence Alignment & Tree Analysis

```bash
python part1/analyze_trees.py --input_dir results_part1 --out_dir results_part1/trees
```



##  PART 2 — Expression Analysis: Human (GTEx) and Mouse (ENCODE)

```bash
python compare_mouse_human_expression.py
```



##  PART 3 — Comparative Expression Across Amniote and Primate Species


##  Amniote/Primate Ortholog Expression

```bash
python compare_kaessmann_amniote_primate.py
  ```


##  PART 4 — Single-cell RNAseq Spermatogenesis Analysis (Green et al., 2019)

```bash
python scripts/spermatogenesis_analysis.py \
  spermatogenesis_scRNA_averageexpression.tsv \
  gene_hgnc_ensmus.tsv
 ```


##  PART 5 — Tumor vs Normal Tissue Analysis (TCGA / GTEx)

In preparation




Citation

If you use this repository or its derived Python analyses, please cite:

Begik O, Lucas MC, Liu H, Ramirez JM, Mattick JS, and Novoa EM.
Integrative analyses of the RNA modification machinery reveal tissue- and cancer-specific signatures.
Genome Biology, 2020.
DOI: 10.1186/s13059-020-02009-z



