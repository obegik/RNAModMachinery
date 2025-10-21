#!/usr/bin/env python3
"""
rnamod_discover_v3.py
Automated RNA modification protein discovery pipeline:
HMM search → filtering → MAFFT alignment → IQ-TREE (fast/robust) → JSON summary
Author: Oguzhan Begik (RNAModMachinery v3)
"""

import os, subprocess, glob, json, re
from Bio import SeqIO
from statistics import mean

# ===============================================================
#                    Helper functions
# ===============================================================

def run_cmd(cmd, dry=False):
    """Run a shell command safely."""
    print("▶️  Running:", " ".join(cmd))
    if not dry:
        subprocess.run(cmd, check=True)

def ensure_dir(path):
    os.makedirs(path, exist_ok=True)

def parse_tblout(tbl_file, evalue_cutoff=1e-5, coverage_cutoff=0.5):
    hits = []
    hmm_length = None
    with open(tbl_file) as f:
        for line in f:
            if line.startswith("#"):
                if "Query:" in line and "[M=" in line:
                    match = re.search(r"\[M=(\d+)\]", line)
                    if match:
                        hmm_length = int(match.group(1))
                continue
            parts = line.strip().split()
            if len(parts) < 18:
                continue
            seq_id = parts[0]
            evalue = float(parts[4])
            hmm_from, hmm_to = int(parts[15]), int(parts[16])
            coverage = abs(hmm_to - hmm_from + 1) / hmm_length if hmm_length else 1
            if evalue < evalue_cutoff and coverage >= coverage_cutoff:
                hits.append(seq_id)
    return hits

def extract_hits(tbl_file, fasta_file, output_fasta, evalue_cutoff=1e-5, coverage_cutoff=0.5):
    ids = parse_tblout(tbl_file, evalue_cutoff, coverage_cutoff)
    if not ids:
        print(f"  ⚠️  No high-confidence hits found in {tbl_file}")
        return 0

    records = [r for r in SeqIO.parse(fasta_file, "fasta") if any(i in r.id for i in ids)]
    if not records:
        print(f"  ⚠️  No sequences matched in {fasta_file}")
        return 0

    SeqIO.write(records, output_fasta, "fasta")
    print(f"✅ Extracted {len(records)} filtered sequences to {output_fasta}")
    return len(records)

def run_hmmsearch(hmm_file, fasta_file, output_tbl, threads=4):
    cmd = ["hmmsearch", "--cpu", str(threads), "--tblout", output_tbl, hmm_file, fasta_file]
    run_cmd(cmd)

def collapse_isoforms(fasta_in, fasta_out):
    """Keep only one sequence per UniProt base accession (e.g. Q9BUB4)."""
    seen = set()
    unique_records = []
    for record in SeqIO.parse(fasta_in, "fasta"):
        match = re.search(r"\|([A-Z0-9]+)(?:-\d+)?\|", record.id)
        if match:
            acc = match.group(1)
            if acc not in seen:
                seen.add(acc)
                unique_records.append(record)
        else:
            unique_records.append(record)
    SeqIO.write(unique_records, fasta_out, "fasta")
    print(f"🧹 Collapsed to {len(unique_records)} unique UniProt accessions → {fasta_out}")
    return len(unique_records)

def run_mafft(input_fasta, output_fasta, threads):
    cmd = ["mafft", "--thread", str(threads), "--auto", input_fasta]
    with open(output_fasta, "w") as out:
        subprocess.run(cmd, stdout=out, check=True)
    print(f"✅ MAFFT alignment completed: {output_fasta}")

def run_iqtree(aligned_fasta, threads, fast_mode=False):
    """Run IQ-TREE in fast or robust mode."""
    cmd = [
        "iqtree",
        "-s", aligned_fasta,
        "-st", "AA",
        "-m", "JTT+F+I+G4",
        "-nt", str(threads),
        "--safe",
        "--redo"
    ]

    if fast_mode:
        cmd.append("--fast")
        print("⚡ Running IQ-TREE in FAST mode (exploratory)")
    else:
        cmd += ["-bb", "100", "-alrt", "100"]
        print("🧬 Running IQ-TREE in ROBUST mode (100 bootstrap + 100 SH-aLRT)")

    run_cmd(cmd)
    print(f"✅ IQ-TREE completed for {aligned_fasta}")

def compute_alignment_stats(fasta_file):
    lengths = [len(rec.seq) for rec in SeqIO.parse(fasta_file, "fasta")]
    return {
        "nseq": len(lengths),
        "mean_length": round(mean(lengths), 2) if lengths else 0,
        "min_length": min(lengths) if lengths else 0,
        "max_length": max(lengths) if lengths else 0
    }

# ===============================================================
#                    Main pipeline
# ===============================================================

def main():
    import argparse
    parser = argparse.ArgumentParser(
        description="Automated RNA modification protein discovery: HMM search → filtering → alignment → tree + JSON summary"
    )
    parser.add_argument("--hmm-dir", required=True, help="Directory containing Pfam HMM files (.hmm)")
    parser.add_argument("--proteomes", required=True, help="Directory containing species proteomes (.fasta)")
    parser.add_argument("--threads", type=int, default=4, help="Number of CPU threads")
    parser.add_argument("--outdir", default="results_part1_v3", help="Output directory")
    parser.add_argument("--evalue", type=float, default=1e-5, help="E-value cutoff for filtering")
    parser.add_argument("--coverage", type=float, default=0.5, help="Coverage cutoff for filtering")
    parser.add_argument("--fast-tree", action="store_true", help="Use IQ-TREE fast mode (skip bootstraps)")
    args = parser.parse_args()

    ensure_dir(args.outdir)
    hmm_files = sorted(glob.glob(os.path.join(args.hmm_dir, "*.hmm")))
    proteomes = sorted(glob.glob(os.path.join(args.proteomes, "*.fasta")))

    print(f"🧩 Found {len(hmm_files)} Pfam domains")
    print(f"🧬 Found {len(proteomes)} proteomes")

    summary_data = {}

    for hmm_file in hmm_files:
        domain = os.path.splitext(os.path.basename(hmm_file))[0]
        domain_out = os.path.join(args.outdir, domain)
        ensure_dir(domain_out)

        print(f"\n=== Processing Pfam domain: {domain} ===")
        merged_hits = os.path.join(domain_out, f"{domain}_merged.fasta")
        collapsed_hits = os.path.join(domain_out, f"{domain}_collapsed.fasta")
        merged_records = []

        # HMMER search for each species
        for fasta_file in proteomes:
            species = os.path.splitext(os.path.basename(fasta_file))[0]
            tbl_file = os.path.join(domain_out, f"{species}.tblout")
            fasta_hits = os.path.join(domain_out, f"{species}_hits.fasta")

            try:
                run_hmmsearch(hmm_file, fasta_file, tbl_file, args.threads)
                nhits = extract_hits(tbl_file, fasta_file, fasta_hits, args.evalue, args.coverage)
                if nhits > 0:
                    merged_records += list(SeqIO.parse(fasta_hits, "fasta"))
            except subprocess.CalledProcessError:
                print(f"❌ HMMER failed for {species}")
                continue

        if not merged_records:
            print(f"⚠️  Skipping {domain}: no hits found across species")
            continue

        SeqIO.write(merged_records, merged_hits, "fasta")
        print(f"✅ Merged {len(merged_records)} total hits → {merged_hits}")

        # Collapse isoforms
        collapse_isoforms(merged_hits, collapsed_hits)

        # MAFFT alignment
        aligned_fasta = os.path.join(domain_out, f"{domain}_aligned.fasta")
        try:
            run_mafft(collapsed_hits, aligned_fasta, args.threads)
        except subprocess.CalledProcessError:
            print(f"❌ MAFFT failed for {domain}")
            continue

        # IQ-TREE (fast or robust)
        try:
            run_iqtree(aligned_fasta, args.threads, fast_mode=args.fast_tree)
        except subprocess.CalledProcessError:
            print(f"❌ IQ-TREE failed for {domain}")
            continue

        stats = compute_alignment_stats(aligned_fasta)
        stats.update({"domain": domain, "total_hits": len(merged_records)})
        summary_data[domain] = stats

        json_file = os.path.join(domain_out, f"{domain}_summary.json")
        with open(json_file, "w") as j:
            json.dump(stats, j, indent=2)
        print(f"📊 Summary written to {json_file}")

    # Master summary
    tsv_path = os.path.join(args.outdir, "summary_all_domains.tsv")
    with open(tsv_path, "w") as out:
        out.write("Domain\tTotalHits\tSequences\tMeanLen\tMinLen\tMaxLen\n")
        for d, s in summary_data.items():
            out.write(f"{d}\t{s['total_hits']}\t{s['nseq']}\t{s['mean_length']}\t{s['min_length']}\t{s['max_length']}\n")

    print(f"\n✅ Master summary written to {tsv_path}")
    print("🎯 RNAmodMachinery v3 pipeline completed successfully.")

if __name__ == "__main__":
    main()