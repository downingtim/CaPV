#!/usr/bin/env python3
import csv
import sys
import os

# 1. Parse the local GenBank file manually to build a clean index of gapless reference proteins
gb_file = "KX894508.gb"
if not os.path.exists(gb_file):
    print(f"Error: {gb_file} not found in the current directory.")
    sys.exit(1)

print(f"Parsing actual gene structures from {gb_file}...")
genes_metadata = []
current_gene = None

with open(gb_file, "r") as f:
    in_translation = False
    current_translation = []
    
    for line in f:
        if '/gene="' in line:
            raw_gene = line.split('/gene="')[1].split('"')[0]
            current_gene = raw_gene.replace("LD", "LSDV")
            
        elif '/translation="' in line:
            in_translation = True
            raw_trans = line.split('/translation="')[1].strip().replace('"', '')
            current_translation = [raw_trans]
            if line.count('"') == 2:
                in_translation = False
                if current_gene:
                    full_seq = "".join(current_translation)
                    full_seq = "".join([c for c in full_seq if c.isalpha()])
                    genes_metadata.append({"id": current_gene, "len": len(full_seq)})
                    
        elif in_translation:
            cleaned_line = line.strip().replace('"', '')
            current_translation.append(cleaned_line)
            if '"' in line:
                in_translation = False
                if current_gene:
                    full_seq = "".join(current_translation)
                    full_seq = "".join([c for c in full_seq if c.isalpha()])
                    genes_metadata.append({"id": current_gene, "len": len(full_seq)})

print(f"Successfully mapped {len(genes_metadata)} reference genes.")

# Initialize the per-gene statistics dictionary to track conserved positions
gene_stats = {g["id"]: {"len": g["len"], "conserved_sites": 0} for g in genes_metadata}
# Spacer for potential insertions relative to the reference sequence line
gene_stats["CaPV_Insertion"] = {"len": 0, "conserved_sites": 0}

# 2. Build map lookup translated to true gapless positions
def lookup_gene_and_local_pos(gapless_ref_position):
    """
    Takes a 1-based GAPLESS index matching the reference sequence string 
    and returns its exact Gene ID and inner-gene relative coordinate.
    """
    cumulative_length = 0
    for gene in genes_metadata:
        start_idx = cumulative_length + 1
        end_idx = cumulative_length + gene["len"]
        if start_idx <= gapless_ref_position <= end_idx:
            local_pos = gapless_ref_position - cumulative_length
            return gene["id"], local_pos
        cumulative_length += gene["len"]
    return "Unknown_Gene", "N/A"

# 3. Read alignment records
alignment_file = "CAPV_ONLY/total_proteome_capv.aln"
if not os.path.exists(alignment_file):
    print(f"Error: {alignment_file} missing.")
    sys.exit(1)

print("Reading alignment matrix...")
sequences = {}
current_id = None
with open(alignment_file, "r") as f:
    for line in f:
        if line.startswith(">"):
            current_id = line.strip()[1:]
            sequences[current_id] = []
        elif current_id:
            sequences[current_id].append(line.strip())
for seq_id in sequences:
    sequences[seq_id] = "".join(sequences[seq_id])

# 4. Target the specific Reference string tracker to convert coordinates accurately
ref_sample_key = None
for seq_id in sequences:
    if "KX894508" in seq_id:
        ref_sample_key = seq_id
        break

print(f"Using reference sequence key: {ref_sample_key} for absolute coordinates.")

alignment_length = len(next(iter(sequences.values())))
ref_sequence_string = sequences[ref_sample_key]
conserved_sites = []
total_conserved_count = 0
gapless_ref_counter = 0
print("Analyzing alignment matrix for 100% invariant conserved columns...")
for i in range(alignment_length):
    ref_char = ref_sequence_string[i]
    if ref_char != '-':
        gapless_ref_counter += 1
        
    # Gather residues from ALL sequences at position i
    all_residues = [seq[i] for seq in sequences.values()]
    unique_residues = set(all_residues)
    
    # CONDITION FOR 100% INVARIANT CONSERVATION:
    # All sequences must share the exact same character, and it cannot be a gap character
    if len(unique_residues) == 1:
        conserved_aa = all_residues[0]
        if conserved_aa == '-':
            continue # Skip columns where the entire dataset contains a gap
    else:
        continue # Found a mutation/difference at this site, skip it

    global_alignment_pos = i + 1
    total_conserved_count += 1
    
    # Convert global alignment position to absolute reference coordinates
    if ref_char != '-':
        gene_id, local_gene_pos = lookup_gene_and_local_pos(gapless_ref_counter)
    else:
        gene_id, local_gene_pos = "CaPV_Insertion", "N/A"

    # Log the site
    conserved_sites.append([global_alignment_pos, gene_id, local_gene_pos, conserved_aa])
    gene_stats[gene_id]["conserved_sites"] += 1

# Ensure output sorting follows the physical path of the genome strings
conserved_sites.sort(key=lambda x: x[0])

# 5. Output values to a master structural location database CSV
output_sites_csv = "within_cppv_conservation.csv"
with open(output_sites_csv, mode='w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["Global_Alignment_Pos", "Gene_ID", "Position_In_Gene", "Conserved_AA"])
    for row in conserved_sites:
        writer.writerow(row)

# 6. Output the requested per-gene summary table
output_genes_csv = "within_cppv_per_gene_summary.csv"
with open(output_genes_csv, mode='w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(["Gene_ID", "Gene_Length_AA", "Total_Conserved_Sites", "Conservation_Density_Pct"])
    
    for g_id, stats in gene_stats.items():
        if g_id == "CaPV_Insertion" and stats["conserved_sites"] == 0:
            continue
            
        if stats["len"] > 0:
            density = (stats["conserved_sites"] / stats["len"]) * 100
        else:
            density = 0.0
            
        writer.writerow([g_id, stats["len"], stats["conserved_sites"], f"{density:.2f}%"])

# 7. Print summary metrics dashboard to standard on-screen terminal
pct_total_genome_conserved = (total_conserved_count / gapless_ref_counter * 100) if gapless_ref_counter > 0 else 0

print("\n" + "="*50)
print("          CONSERVATION REPORT DASHBOARD          ")
print("="*50)
print(f"Total Alignment Column Length  : {alignment_length:,} columns")
print(f"Total Reference Proteins Length: {gapless_ref_counter:,} AA")
print("-"*50)
print(f"TOTAL 100% INVARIANT SITES FOUND: {total_conserved_count:,}")
print(f"Total Reference Genome Invariant: {pct_total_genome_conserved:.1f}%")
print("="*50)
print(f"Saved site-by-site database to  : {output_sites_csv}")
print(f"Saved per-gene summary table to : {output_genes_csv}\n")
