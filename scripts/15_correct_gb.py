#!/usr/bin/env python3
"""
link_genes.py

Link gene/CDS identifiers between a reference GenBank file and target GenBank files,
cluster accessory genes via BLASTn, and generate an UpSetR overlap plot.

Features:
- Synteny Constraint: Genes are prevented from cross-mapping to opposite ends.
- Split/Merge support: Multiple hits to a ref gene get numerical suffixes (e.g., LD012.1).
- Novel CDS support: Unmatched genes get alphabetical suffixes (e.g., LD012a).
- UpSetR plotting with large publication-ready fonts and colors.
"""

import argparse
import csv
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

from Bio import SeqIO

# --------------------------------------------------------------------------
# Feature extraction & Positions
# --------------------------------------------------------------------------

def label_for_feature(feat, fallback_index):
    q = feat.qualifiers
    for key in ("locus_tag", "gene", "protein_id"):
        if key in q and q[key]:
            return q[key][0]
    return f"feature_{fallback_index}"

def extract_genes(gb_path, feature_type):
    try:
        records = list(SeqIO.parse(str(gb_path), "genbank"))
    except Exception as e:
        sys.exit(f"FAILED TO PARSE: {gb_path}\nError: {e}")
    if not records:
        raise ValueError(f"No records found in {gb_path}")
    record = records[0]

    genes = []
    idx = 0
    seq_len = len(record.seq)
    
    for feat in record.features:
        if feat.type != feature_type:
            continue
        idx += 1
        seq = feat.extract(record.seq)
        
        start = int(feat.location.start)
        end = int(feat.location.end)
        
        genes.append({
            "label": label_for_feature(feat, idx),
            "feature": feat,
            "seq": str(seq),
            "start": start,
            "end": end,
            # Calculate the relative center (0.0 to 1.0) of the gene
            "rel_pos": ((start + end) / 2.0) / seq_len 
        })
    genes.sort(key=lambda g: g["start"])

    # De-duplicate labels
    seen = {}
    for g in genes:
        base = g["label"]
        seen[base] = seen.get(base, 0) + 1
        if seen[base] > 1:
            g["label"] = f"{base}_{seen[base]}"

    return record, genes


# --------------------------------------------------------------------------
# BLAST & Clustering
# --------------------------------------------------------------------------

def check_blast_available():
    missing = [exe for exe in ("blastn", "makeblastdb") if shutil.which(exe) is None]
    if missing:
        sys.exit(f"ERROR: required BLAST+ executable(s) not found on PATH: {missing}")

def write_fasta(genes, path):
    with open(path, "w") as fh:
        for g in genes:
            fh.write(f">{g['label']}\n{g['seq']}\n")

BLAST_FIELDS = "qseqid sseqid pident length nident mismatch gapopen " \
               "qstart qend sstart send evalue bitscore qlen slen"

def run_blastn(query_fasta, subject_fasta, workdir, threads=1):
    db_prefix = str(workdir / "subject_db")
    subprocess.run(
        ["makeblastdb", "-in", str(subject_fasta), "-dbtype", "nucl", "-out", db_prefix],
        check=True, capture_output=True, text=True
    )
    result = subprocess.run(
        ["blastn", "-task", "blastn", "-query", str(query_fasta), "-db", db_prefix,
         "-outfmt", f"6 {BLAST_FIELDS}", "-num_threads", str(threads)],
        check=True, capture_output=True, text=True
    )
    hits = []
    for line in result.stdout.strip().splitlines():
        if not line: continue
        vals = line.split("\t")
        rec = dict(zip(BLAST_FIELDS.split(), vals))
        for key in ("pident", "length", "nident", "mismatch", "gapopen",
                    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
                    "qlen", "slen"):
            rec[key] = float(rec[key])
        hits.append(rec)
    return hits


def assign_hits_multi(hits, threshold, max_rel_dist, ref_genes, tgt_genes):
    """
    Finds the BEST synteny-valid reference match for each target gene.
    Allows multiple target genes to hit the same reference gene (splits).
    """
    ref_pos = {g["label"]: g["rel_pos"] for g in ref_genes}
    tgt_pos = {g["label"]: g["rel_pos"] for g in tgt_genes}
    
    valid_hits = []
    for h in hits:
        shorter_len = min(h["qlen"], h["slen"])
        if shorter_len <= 0: continue
        
        pct_id = (h["nident"] / shorter_len) * 100.0
        if pct_id < threshold: continue
        
        q, s = h["qseqid"], h["sseqid"] # q is Ref, s is Tgt
        
        # Synteny Constraint: Reject matches that jump across the genome
        if abs(ref_pos[q] - tgt_pos[s]) > max_rel_dist:
            continue
            
        h["pct_id_shorter"] = pct_id
        valid_hits.append(h)
        
    # Sort by Bitscore (desc), then Identity (desc)
    valid_hits.sort(key=lambda x: (x["bitscore"], x["pct_id_shorter"]), reverse=True)
    
    # Assign target gene to its top valid reference hit
    best_for_target = {}
    for h in valid_hits:
        s = h["sseqid"]
        q = h["qseqid"]
        if s not in best_for_target:
            best_for_target[s] = q
            
    return best_for_target

def cluster_unmatched(items, threshold, workdir, threads):
    if not items: return []
    if len(items) == 1: return [[items[0]]]
    
    fasta_path = workdir / "temp_unmatched.fasta"
    with open(fasta_path, "w") as fh:
        for i, (ti, g, col_name) in enumerate(items):
            fh.write(f">seq_{i}\n{g['seq']}\n")
            
    hits = run_blastn(fasta_path, fasta_path, workdir, threads=threads)
    edges = defaultdict(list)
    for h in hits:
        q = int(h["qseqid"].split("_")[1])
        s = int(h["sseqid"].split("_")[1])
        if q == s: continue
        shorter_len = min(h["qlen"], h["slen"])
        if shorter_len > 0 and ((h["nident"] / shorter_len) * 100.0) >= threshold:
            edges[q].append(s)
            edges[s].append(q)
                
    visited = set()
    components = []
    for i in range(len(items)):
        if i not in visited:
            comp = []
            queue = [i]
            visited.add(i)
            while queue:
                curr = queue.pop(0)
                comp.append(curr)
                for neighbor in edges[curr]:
                    if neighbor not in visited:
                        visited.add(neighbor)
                        queue.append(neighbor)
            components.append([items[idx] for idx in sorted(comp)])
    return components

def get_letter_suffix(idx):
    """ Converts 1 -> 'a', 2 -> 'b', 27 -> 'aa' """
    chars = "abcdefghijklmnopqrstuvwxyz"
    res = ""
    idx -= 1
    while idx >= 0:
        res = chars[idx % 26] + res
        idx = (idx // 26) - 1
    return res

# --------------------------------------------------------------------------
# Main pipeline
# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reference", required=True, type=Path)
    ap.add_argument("--targets", required=True, type=Path, nargs="+")
    ap.add_argument("--outdir", required=True, type=Path)
    ap.add_argument("--labels", type=str, nargs="+", 
                    help="Custom labels for the genomes (Ref + Targets)")
    ap.add_argument("--plot-order", type=str, nargs="+",
                    help="Exact order of labels to use in the UpSet plot")
    ap.add_argument("--threshold", type=float, default=69.9)
    ap.add_argument("--max-rel-dist", type=float, default=0.20,
                    help="Max relative genome distance allowed for a match")
    ap.add_argument("--feature-type", default="CDS")
    ap.add_argument("--id-prefix", default="LD")
    ap.add_argument("--threads", type=int, default=1)
    args = ap.parse_args()

    check_blast_available()
    args.outdir.mkdir(parents=True, exist_ok=True)
    
    total_genomes = 1 + len(args.targets)
#    if args.labels and len(args.labels) != total_genomes:
#       sys.exit(f"ERROR: Provided {len(args.labels)} labels, but there are {total_genomes} genomes.")

    print(f"Parsing reference: {args.reference}")
    ref_record, ref_genes = extract_genes(args.reference, args.feature_type)
    width = max(3, len(str(len(ref_genes))))
    for i, g in enumerate(ref_genes, start=1):
        g["id"] = f"{args.id_prefix}{i:0{width}d}"
    ref_id_by_label = {g["label"]: g["id"] for g in ref_genes}

    target_data = []
    for tpath in args.targets:
        print(f"Parsing target: {tpath}")
        rec, genes = extract_genes(tpath, args.feature_type)
        target_data.append((rec, genes, tpath))

    table_rows = {g["id"]: {"ref_label": g["label"]} for g in ref_genes}
    target_assigned_ids = [dict() for _ in target_data]
    unmatched_pool = [] 

    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        ref_fasta = tmp / "reference_genes.fasta"
        write_fasta(ref_genes, ref_fasta)

        for ti, (rec, genes, tpath) in enumerate(target_data):
            col_name = rec.id or Path(tpath).stem
            print(f"BLASTing reference vs {col_name} ...")
            target_fasta = tmp / f"target_{ti}.fasta"
            write_fasta(genes, target_fasta)

            hits = run_blastn(ref_fasta, target_fasta, tmp, threads=args.threads)
            
            # Map target genes to their best valid reference hits
            best_target_to_ref = assign_hits_multi(
                hits, args.threshold, args.max_rel_dist, ref_genes, genes
            )

            # Group target genes by matched Ref ID to handle splits (.1, .2)
            ref_to_tgt = defaultdict(list)
            for g in genes:
                if g["label"] in best_target_to_ref:
                    ref_label = best_target_to_ref[g["label"]]
                    ref_id = ref_id_by_label[ref_label]
                    ref_to_tgt[ref_id].append(g)

            # Assign row IDs (1st match -> LD012, 2nd match -> LD012.1)
            for ref_id, t_genes in ref_to_tgt.items():
                t_genes.sort(key=lambda x: x["start"])
                for i, g in enumerate(t_genes):
                    row_id = ref_id if i == 0 else f"{ref_id}.{i}"
                    if row_id not in table_rows:
                        table_rows[row_id] = {"ref_label": ""}
                    table_rows[row_id][col_name] = g["label"]
                    target_assigned_ids[ti][g["label"]] = row_id

            # Keep track of novel/unmatched genes and their upstream neighbors
            prev_base_id = f"{args.id_prefix}000" 
            for g in genes:
                if g["label"] in best_target_to_ref:
                    # Update nearest anchor
                    full_id = target_assigned_ids[ti][g["label"]]
                    prev_base_id = re.match(r"^([A-Za-z]+\d+)", full_id).group(1)
                else:
                    unmatched_pool.append((ti, g, col_name, prev_base_id))

        if unmatched_pool:
            print("\nCross-comparing unassigned ORFs across all targets...")
            unmatched_by_locus = defaultdict(list)
            for item in unmatched_pool:
                prev_base_id = item[3]
                unmatched_by_locus[prev_base_id].append(item[:3]) 

            # Assign unique suffix letters (a, b) to novel gene clusters
            for prev_id, items in unmatched_by_locus.items():
                clusters = cluster_unmatched(items, args.threshold, tmp, args.threads)
                for c_idx, cluster in enumerate(clusters, start=1):
                    letter_suffix = get_letter_suffix(c_idx)
                    new_id = f"{prev_id}{letter_suffix}"
                    table_rows.setdefault(new_id, {"ref_label": ""})
                    
                    for ti, g, col_name in cluster:
                        target_assigned_ids[ti][g["label"]] = new_id
                        existing = table_rows[new_id].get(col_name, "")
                        # Handle very rare within-target repeats of novel genes
                        if existing:
                            table_rows[new_id][col_name] = f"{existing}; {g['label']}"
                        else:
                            table_rows[new_id][col_name] = g["label"]

    # ---- Writing Outputs ----
    col_names = [rec.id or Path(tpath).stem for rec, _, tpath in target_data]
    plot_labels = args.labels if args.labels else [args.reference.stem] + col_names

    # Advanced sort: LD012 < LD012.1 < LD012a < LD012b
    def sort_key(gid):
        m = re.match(r"([A-Za-z]+)(\d+)(?:\.(\d+))?([a-z]+)?", gid)
        if not m: return (gid, 0, "", 0)
        prefix = m.group(1)
        num = int(m.group(2))
        split_num = int(m.group(3)) if m.group(3) else -1
        letter = m.group(4) if m.group(4) else ""
        return (prefix, num, letter, split_num)

    csv_path = args.outdir / "gene_link_table.csv"
    pa_path = args.outdir / "presence_absence_matrix.csv"
    
    with open(csv_path, "w", newline="") as fh_csv, open(pa_path, "w", newline="") as fh_pa:
        writer_csv = csv.writer(fh_csv)
        writer_pa = csv.writer(fh_pa)
        
        writer_csv.writerow(["gene_id", f"{args.reference.stem}_label", *col_names])
        writer_pa.writerow(["gene_id"] + plot_labels)
        
        for gid in sorted(table_rows.keys(), key=sort_key):
            row = table_rows[gid]
            writer_csv.writerow([gid, row.get("ref_label", ""), *[row.get(c, "") for c in col_names]])
            
            pa_row = [gid]
            pa_row.append(1 if row.get("ref_label") else 0) 
            for c in col_names:
                pa_row.append(1 if row.get(c) else 0)       
            writer_pa.writerow(pa_row)
            
    print(f"\nWrote tables: {csv_path} and {pa_path}")

    # Write GenBank files
    def annotate_and_write(record, genes, id_map, out_path):
        for g in genes:
            gid = id_map.get(g["label"])
            if gid:
                g["feature"].qualifiers["gene"] = [gid]
        SeqIO.write(record, str(out_path), "genbank")

    annotate_and_write(ref_record, ref_genes, ref_id_by_label, args.outdir / f"{args.reference.stem}.linked.gb")
    for (rec, genes, tpath), id_map in zip(target_data, target_assigned_ids):
        annotate_and_write(rec, genes, id_map, args.outdir / f"{Path(tpath).stem}.linked.gb")

    # ---- UpSet Plot ----
    if shutil.which("Rscript"):
        print("\nGenerating UpSet plot via Rscript...")
        
        # Configure plotting order directly from args
        if args.plot_order:
            plot_order_str = ", ".join([f'"{x}"' for x in args.plot_order])
            sets_declaration = f"plot_sets <- c({plot_order_str})"
        else:
            sets_declaration = "plot_sets <- colnames(df)[-1]"
            
        r_script_path = args.outdir / "plot_upset.R"
        png_path = args.outdir / "upset_plot.png"
        
        r_code = f"""
        if (!require("UpSetR", quietly = TRUE)) install.packages("UpSetR", repos="http://cran.us.r-project.org")
        library(UpSetR)

        df <- read.csv("{pa_path.name}", check.names=FALSE)
        {sets_declaration}

        # Publication-ready Okabe-Ito inspired palette
        base_colors <- c("#D55E00", "#0072B2", "#009E73", "#E69F00", "#CC79A7", "#56B4E9", "#F0E442", "#999999")
        set_colors <- base_colors[1:length(plot_sets)]

        # Bumped dimensions (12x8) for large fonts
        png("{png_path.name}", width=12, height=8, units="in", res=300)
        p <- upset(df, sets=plot_sets, keep.order=TRUE, order.by="freq",
                   mainbar.y.label = "Number of Shared CDS",
                   sets.x.label = "Total CDS per Genome",
                   sets.bar.color = set_colors,
                   point.size = 4, line.size = 1.5,
                   text.scale = c(2.0, 1.8, 2.0, 1.5, 2.0, 1.6))
        print(p)
        dev.off()
        """
        with open(r_script_path, "w") as rf:
            rf.write(r_code)
            
        subprocess.run(["Rscript", str(r_script_path.name)], cwd=str(args.outdir))
        print(f"Created UpSet plot: {png_path}")

    print("Done.")

if __name__ == "__main__":
    main()
