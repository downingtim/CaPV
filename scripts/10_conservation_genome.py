#!/usr/bin/env python3
import csv
import re

# 1. Faster CDS Parser
def get_gene_features(gbk_file):
    genes = []
    with open(gbk_file, "r") as f:
        content = f.read()
        for part in re.split(r'\n\s{5}(?=CDS)', content):
            gene_match = re.search(r'/gene="([^"]+)"', part)
            cds_match = re.search(r'CDS\s+(?:complement\()?<?(\d+)\.\.>?(\d+)', part)
            if cds_match:
                start, end = int(cds_match.group(1)), int(cds_match.group(2))
                genes.append({
                    "id": gene_match.group(1) if gene_match else "unknown",
                    "start": start, "end": end
                })
    return genes

print("Loading data...")
genes = get_gene_features("KX894508.gb")
seqs = {}
with open("genome_capv_98.aln", "r") as f:
    curr = None
    for line in f:
        if line.startswith(">"): curr = line.strip()[1:]
        elif curr: seqs.setdefault(curr, []).append(line.strip().upper())
for k in seqs: seqs[k] = "".join(seqs[k])

ref_key = next((k for k in seqs if "KX894508" in k), None)
ref_aln = seqs[ref_key]
seq_list = list(seqs.values())

# Map alignment columns to reference positions
pos_tracker = [None] * len(ref_aln)
count = 0
for i, char in enumerate(ref_aln):
    if char != '-':
        count += 1
        pos_tracker[i] = count

# 2. Analyze Alignment Footprint
# stats track: [conserved_site_count, total_aligned_columns]
gene_stats = {g["id"]: [0, 0] for g in genes}
gene_stats["Intergenic"] = [0, 0]

conserved_sites = []
print("Analyzing columns...")
for i, column in enumerate(zip(*seq_list)):
    ref_pos = pos_tracker[i]
    
    # Determine Gene ID
    if ref_pos:
        gid = next((g["id"] for g in genes if g["start"] <= ref_pos <= g["end"]), "Intergenic")
    else:
        # If position is an insertion relative to reference, use the preceding gene
        prev_pos = next((pos_tracker[j] for j in range(i-1, -1, -1) if pos_tracker[j] is not None), None)
        gid = next((g["id"] for g in genes if g["start"] <= prev_pos <= g["end"]), "Intergenic") if prev_pos else "Intergenic"

    # Track alignment footprint
    gene_stats[gid][1] += 1
    
    # Track conservation
    if len(set(column)) == 1 and column[0] != '-':
        gene_stats[gid][0] += 1
        if ref_pos: conserved_sites.append([ref_pos, gid, column[0]])

# 3. Export
with open("universal_per_gene_summary.csv", "w", newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["Gene_ID", "Aln_Len", "Conserved_Sites", "Conservation_Density"])
    for gid, stats in gene_stats.items():
        cons, aln_len = stats
        density = (cons / aln_len) if aln_len > 0 else 0
        writer.writerow([gid, aln_len, cons, f"{density:.4f}"])

print(f"Analysis complete. Found {len(conserved_sites)} universally conserved sites.")
