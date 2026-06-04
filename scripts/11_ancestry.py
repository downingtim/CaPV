#!/usr/bin/env python3
import os
import csv
import random
import re
from collections import defaultdict

# Core Configuration
GB_FILE = "KX894508.gb"

# Complete Standard Genetic Code Dictionary
GENETIC_CODE = {
    'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
    'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
    'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
    'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
    'TCT':'S', 'TCC':'S', 'TCG':'S', 'TCA':'S',
    'CCT':'P', 'CCC':'P', 'CCG':'P', 'CCA':'P',
    'ACT':'T', 'ACC':'T', 'ACG':'T', 'ACA':'T',
    'GCT':'A', 'GCC':'A', 'GCG':'A', 'GCA':'A',
    'TAT':'Y', 'TAC':'Y', 'TAA':'_', 'TAG':'_',
    'CAT':'H', 'CAC':'H', 'CAG':'Q', 'CAA':'Q',
    'AAT':'N', 'AAC':'N', 'AAG':'K', 'AAA':'K',
    'GAT':'D', 'GAC':'D', 'GAG':'E', 'GAA':'E',
    'TGT':'C', 'TGC':'C', 'TGA':'_', 'TGG':'W',
    'CGT':'R', 'CGC':'R', 'CGG':'R', 'CGA':'R',
    'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R'
}

COMPLEMENT = str.maketrans("ACGT", "TGCA")

# Master Mutation Tracker: {taxon: {genome_pos: {detailed_metadata_dict}}}
mutation_store = {}

def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]

def parse_genbank_rich(gb_path):
    """Parses GenBank records tracking explicit CDS locations and strand directions."""
    genes_metadata = {}
    position_to_gene = {} # Map: genome_index (0-based) -> (gene_id, strand, start_0, end_0)
    origin_seq = []
    in_origin = False
    
    if not os.path.exists(gb_path):
        raise FileNotFoundError(f"Missing required GenBank reference genome file: {gb_path}")
        
    with open(gb_path, 'r') as f:
        current_coords = None
        strand = "+"
        
        for line in f:
            if line.startswith("ORIGIN"):
                in_origin = True
                continue
            if in_origin:
                clean_seq = "".join(c for c in line.strip() if c.isalpha())
                if clean_seq:
                    origin_seq.append(clean_seq.upper())
                continue
                
            line_strip = line.strip()
            if line.startswith("     CDS             "):
                strand = "-" if "complement" in line_strip else "+"
                match = re.findall(r'\d+', line_strip)
                if match and len(match) >= 2:
                    current_coords = (int(match[0]), int(match[1]))
            
            if current_coords and ("/locus_tag=" in line_strip or "/gene=" in line_strip):
                gene_id = line_strip.split("=")[1].replace('"', '')
                start, end = current_coords
                genes_metadata[gene_id] = {'start': start, 'end': end, 'strand': strand}
                
                for pos in range(start - 1, end):
                    position_to_gene[pos] = (gene_id, strand, start - 1, end)
                current_coords = None
                
    full_genome_seq = "".join(origin_seq)
    return full_genome_seq, genes_metadata, position_to_gene

def process_lineage(aln_file, taxon_key, summary_csv_name, genome_seq, position_to_gene, genes_metadata):
    """Processes alignments, matches true KX894508 coordinates using unique sliding anchors."""
    mutation_store[taxon_key] = {}
    
    if not os.path.exists(aln_file):
        print(f"Skipping task: Alignment path not found: {aln_file}")
        return

    seqs = {}
    curr = None
    with open(aln_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                curr = line.strip()[1:]
                seqs[curr] = ""
            else:
                seqs[curr] += line.strip().upper()

    ref_key = next((k for k in seqs.keys() if "LSDV" in k.upper() or "KX894508" in k), None)
    if not ref_key:
        raise KeyError(f"Unable to find valid reference tracking header containing 'LSDV' or 'KX894508' inside {aln_file}")

    ref_aln = seqs[ref_key]
    taxon_aln = seqs[taxon_key]
    capv_aln = seqs["CaPV"]
    
    # FIXED: Find a 100% UNIQUE anchor point to avoid getting pulled into wrong ITR regions
    ref_seq_ungapped = ref_aln.replace('-', '')
    genome_start_pos = -1
    
    if ref_seq_ungapped.upper() in genome_seq.upper():
        if genome_seq.upper().count(ref_seq_ungapped.upper()) == 1:
            genome_start_pos = genome_seq.upper().find(ref_seq_ungapped.upper())
            
    if genome_start_pos == -1:
        seed_size = min(50, len(ref_seq_ungapped))
        for step in range(0, len(ref_seq_ungapped) - seed_size + 1, 5):
            seed = ref_seq_ungapped[step:step+seed_size]
            # Ensure the anchor seed occurs exactly once across the master genome
            if genome_seq.upper().count(seed) == 1:
                hit = genome_seq.upper().find(seed)
                genome_start_pos = hit - step
                if 0 <= genome_start_pos < len(genome_seq):
                    break

    if genome_start_pos == -1:
        print(f"WARNING: Reference alignment mapping failed completely for {aln_file}. Defaulting to index 0.")
        genome_start_pos = 0
    else:
        print(f"--> Map Success: Anchored '{aln_file}' uniquely to genome coordinate position: {genome_start_pos + 1}")

    # FIXED: Map alignment positions with a hard wall ceiling matching the physical genome size
    ref_idx_to_genome_pos = []
    current_g_pos = genome_start_pos
    for char in ref_aln:
        if char != '-':
            if 0 <= current_g_pos < len(genome_seq):
                ref_idx_to_genome_pos.append(current_g_pos)
            else:
                ref_idx_to_genome_pos.append(None) # Instantly discard out of bounds coordinates
            current_g_pos += 1
        else:
            ref_idx_to_genome_pos.append(None)

    genome_pos_to_aln_idx = {g_pos: aln_idx for aln_idx, g_pos in enumerate(ref_idx_to_genome_pos) if g_pos is not None}

    detailed_filename = f"{taxon_key}_detailed_mutations.csv"
    detailed_writer = csv.writer(open(detailed_filename, "w", newline=''))
    detailed_writer.writerow([
        "Gene_ID", "Genome_Position", "Gene_Position_NT", "Position_AA", "Mutation_Type", 
        "CaPV_Base_CDS", f"{taxon_key}_Base_CDS", "CaPV_Codon_CDS", f"{taxon_key}_Codon_CDS", 
        "CaPV_AA", f"{taxon_key}_AA"
    ])

    gene_diffs = defaultdict(int)
    gene_lengths = {gid: (m['end'] - m['start'] + 1) for gid, m in genes_metadata.items()}
    
    class_counts = {
        "Synonymous": 0, "Nonsynonymous": 0, "Intergenic": 0,
        "Stop_Creation": 0, "Stop_Loss": 0, "Incomplete_Codon_Context": 0
    }

    for i in range(len(ref_aln)):
        g_pos = ref_idx_to_genome_pos[i]
        if g_pos is None:
            continue
            
        if taxon_aln[i] != '-' and capv_aln[i] != '-':
            if taxon_aln[i] != capv_aln[i]:
                pos_meta = position_to_gene.get(g_pos)
                
                if pos_meta:
                    gene_id, strand, start_0, end_0 = pos_meta
                    gene_diffs[gene_id] += 1
                    
                    if strand == "+":
                        gene_pos_nt = g_pos - start_0 + 1
                        frame_pos = (gene_pos_nt - 1) % 3
                        pos_aa = (gene_pos_nt - 1) // 3 + 1
                        p0 = g_pos - frame_pos
                        p1 = p0 + 1
                        p2 = p0 + 2
                    else:
                        gene_pos_nt = (end_0 - 1) - g_pos + 1
                        frame_pos = (gene_pos_nt - 1) % 3
                        pos_aa = (gene_pos_nt - 1) // 3 + 1
                        n0 = ((gene_pos_nt - 1) // 3) * 3 + 1
                        p0 = (end_0 - 1) - (n0 - 1)
                        p1 = p0 - 1
                        p2 = p0 - 2

                    def get_cds_base(strand_dir, pos_coord, aln_track):
                        if pos_coord in genome_pos_to_aln_idx:
                            b = aln_track[genome_pos_to_aln_idx[pos_coord]]
                        else:
                            return '-'
                        if b == '-': return '-'
                        return b.translate(COMPLEMENT) if strand_dir == '-' else b

                    capv_codon_cds = get_cds_base(strand, p0, capv_aln) + get_cds_base(strand, p1, capv_aln) + get_cds_base(strand, p2, capv_aln)
                    taxon_codon_cds = get_cds_base(strand, p0, taxon_aln) + get_cds_base(strand, p1, taxon_aln) + get_cds_base(strand, p2, taxon_aln)
                    
                    capv_base_cds = get_cds_base(strand, g_pos, capv_aln)
                    derived_base_cds = get_cds_base(strand, g_pos, taxon_aln)
                    
                    if '-' in capv_codon_cds or '-' in taxon_codon_cds:
                        mut_type = "Incomplete_Codon_Context"
                        capv_aa = "-"
                        derived_aa = "-"
                        class_counts["Incomplete_Codon_Context"] += 1
                    else:
                        capv_aa = GENETIC_CODE.get(capv_codon_cds, 'X')
                        derived_aa = GENETIC_CODE.get(taxon_codon_cds, 'X')
                        
                        if capv_aa == derived_aa:
                            mut_type = "Synonymous"
                            class_counts["Synonymous"] += 1
                        else:
                            mut_type = "Nonsynonymous"
                            class_counts["Nonsynonymous"] += 1
                            if capv_aa != '_' and derived_aa == '_':
                                class_counts["Stop_Creation"] += 1
                            elif capv_aa == '_' and derived_aa != '_':
                                class_counts["Stop_Loss"] += 1
                else:
                    gene_id = "Intergenic"
                    gene_pos_nt = "-"
                    pos_aa = "-"
                    mut_type = "Intergenic"
                    capv_base_cds = capv_aln[i]
                    derived_base_cds = taxon_aln[i]
                    capv_codon_cds, taxon_codon_cds = "-", "-"
                    capv_aa, derived_aa = "-", "-"
                    class_counts["Intergenic"] += 1
                
                detailed_writer.writerow([
                    gene_id, g_pos + 1, gene_pos_nt, pos_aa, mut_type,
                    capv_base_cds, derived_base_cds, capv_codon_cds, taxon_codon_cds,
                    capv_aa, derived_aa
                ])
                
                mutation_store[taxon_key][g_pos] = {
                    'gene_id': gene_id, 'mut_type': mut_type, 'pos_genomic': g_pos + 1,
                    'pos_gene_nt': gene_pos_nt, 'pos_aa': pos_aa,
                    'capv_nt': capv_base_cds, 'derived_nt': derived_base_cds,
                    'capv_codon': capv_codon_cds, 'derived_codon': taxon_codon_cds,
                    'capv_aa': capv_aa, 'derived_aa': derived_aa,
                    'length': gene_lengths.get(gene_id, 0)
                }

    summary_writer = csv.writer(open(summary_csv_name, "w", newline=''))
    summary_writer.writerow(["Gene_ID_or_Class", "NT_Length", f"Diffs_{taxon_key}_to_CaPV", "Divergence_%_or_Rate"])
    
    total_genome_size = len(genome_seq)
    total_genic_size = sum(gene_lengths.values())
    total_intergenic_size = total_genome_size - total_genic_size
    
    summary_writer.writerow(["--- FUNCTIONAL BREAKDOWN MATRIX ---", "", "", ""])
    summary_writer.writerow(["CLASS: Intergenic", total_intergenic_size, class_counts["Intergenic"], f"{class_counts['Intergenic']/total_intergenic_size:.6f}" if total_intergenic_size > 0 else 0])
    summary_writer.writerow(["CLASS: Synonymous", total_genic_size, class_counts["Synonymous"], f"{class_counts['Synonymous']/total_genic_size:.6f}" if total_genic_size > 0 else 0])
    summary_writer.writerow(["CLASS: Nonsynonymous", total_genic_size, class_counts["Nonsynonymous"], f"{class_counts['Nonsynonymous']/total_genic_size:.6f}" if total_genic_size > 0 else 0])
    summary_writer.writerow(["CLASS: Nonsynonymous (Stop Codon Creation)", total_genic_size, class_counts["Stop_Creation"], f"{class_counts['Stop_Creation']/total_genic_size:.6f}" if total_genic_size > 0 else 0])
    summary_writer.writerow(["CLASS: Nonsynonymous (Stop Codon Loss)", total_genic_size, class_counts["Stop_Loss"], f"{class_counts['Stop_Loss']/total_genic_size:.6f}" if total_genic_size > 0 else 0])
    summary_writer.writerow(["CLASS: Incomplete Codon Context (Gaps/No Info)", total_genic_size, class_counts["Incomplete_Codon_Context"], f"{class_counts['Incomplete_Codon_Context']/total_genic_size:.6f}" if total_genic_size > 0 else 0])
    summary_writer.writerow(["-----------------------------------", "", "", ""])
    
    for gid in sorted(gene_lengths.keys()):
        diffs = gene_diffs[gid]
        l = gene_lengths[gid]
        summary_writer.writerow([gid, l, diffs, f"{diffs/l:.4f}" if l > 0 else 0])

def write_shared_mutations_filtered(genome_seq, genes_metadata):
    """Compiles strictly meaningful convergent changes, discarding variations that change the AA destination."""
    pairs = [
        ("LSDV", "SPPV"), ("LSDV", "SPPV_modern"), ("LSDV", "GTPV"), 
        ("GTPV", "SPPV"), ("GTPV", "SPPV_modern")
    ]
    sppv_anc = mutation_store.get("SPPV", {})
    
    total_genome_size = len(genome_seq)
    total_genic_size = sum((m['end'] - m['start'] + 1) for m in genes_metadata.values())
    total_intergenic_size = total_genome_size - total_genic_size

    shared_metrics = {}

    f_all = open("all_shared_mutations.csv", "w", newline='')
    f_diff_base = open("shared_same_aa_diff_base.csv", "w", newline='')
    
    w_all = csv.writer(f_all)
    w_diff = csv.writer(f_diff_base)
    
    headers = [
        "Gene_ID", "Gene_Length_NT", "Comparison_L1", "Comparison_L2", 
        "Genome_Site_Pos", "Gene_Position_NT", "Position_AA", "Mutation_Type", 
        "AA_Evolution_Classification", "Base_Match_Status", "Base_L1_CDS", "Base_L2_CDS", "CaPV_Base_CDS", 
        "Codon_L1_CDS", "Codon_L2_CDS", "CaPV_Codon_CDS", 
        "AA_L1", "AA_L2", "CaPV_AA"
    ]
    w_all.writerow(headers)
    w_diff.writerow(headers)
    
    for t1, t2 in pairs:
        comp_key = f"{t1} vs {t2}"
        shared_metrics[comp_key] = {
            "Intergenic_Same_Base": 0,
            "Synonymous_Same_Base": 0,
            "Nonsynonymous_Same_Base_Same_AA": 0,
            "Nonsynonymous_Diff_Base_Same_AA": 0,
            "Stop_Creation": 0, "Stop_Loss": 0
        }
        
        muts1, muts2 = mutation_store.get(t1, {}), mutation_store.get(t2, {})
        shared_positions = set(muts1.keys()) & set(muts2.keys())
        
        for pos in sorted(shared_positions):
            m1, m2 = muts1[pos], muts2[pos]
            
            if m1['mut_type'] == 'Incomplete_Codon_Context' or m2['mut_type'] == 'Incomplete_Codon_Context':
                continue
            
            if t1 == "SPPV_modern" and pos in sppv_anc and sppv_anc[pos]['derived_nt'] == m1['derived_nt']:
                continue
            if t2 == "SPPV_modern" and pos in sppv_anc and sppv_anc[pos]['derived_nt'] == m2['derived_nt']:
                continue
                
            bases_match = (m1['derived_nt'] == m2['derived_nt'])
            
            if m1['mut_type'] == 'Intergenic':
                if not bases_match: continue
                shared_metrics[comp_key]["Intergenic_Same_Base"] += 1
                aa_class, base_status = "Silent", "Same_Base"
                
            elif m1['mut_type'] == 'Synonymous':
                if not bases_match: continue
                shared_metrics[comp_key]["Synonymous_Same_Base"] += 1
                aa_class, base_status = "Silent", "Same_Base"
                
            else: # Nonsynonymous
                if m1['derived_aa'] == m2['derived_aa']:
                    if bases_match:
                        shared_metrics[comp_key]["Nonsynonymous_Same_Base_Same_AA"] += 1
                        base_status = "Same_Base"
                    else:
                        shared_metrics[comp_key]["Nonsynonymous_Diff_Base_Same_AA"] += 1
                        base_status = "Different_Base"
                    aa_class = "Amino_acid"
                    
                    if m1['capv_aa'] != '_' and m1['derived_aa'] == '_':
                        shared_metrics[comp_key]["Stop_Creation"] += 1
                    elif m1['capv_aa'] == '_' and m1['derived_aa'] != '_':
                        shared_metrics[comp_key]["Stop_Loss"] += 1
                else:
                    continue
                        
            row = [
                m1['gene_id'], m1['length'], t1, t2, 
                m1['pos_genomic'], m1['pos_gene_nt'], m1['pos_aa'], m1['mut_type'],
                aa_class, base_status, m1['derived_nt'], m2['derived_nt'], m1['capv_nt'],
                m1['derived_codon'], m2['derived_codon'], m1['capv_codon'],
                m1['derived_aa'], m2['derived_aa'], m1['capv_aa']
            ]
            
            w_all.writerow(row)
            if aa_class == "Amino_acid" and not bases_match:
                w_diff.writerow(row)
                
    f_all.close()
    f_diff_base.close()

    with open("shared_mutations_summary.csv", "w", newline='') as sf:
        sw = csv.writer(sf)
        sw.writerow(["Comparison_Pair", "Metric_Category", "Shared_Count", "Rate_Per_Available_Sites"])
        
        print(f"\n{'Comparison Pair':<20} | {'Metric Category':<45} | {'Count':<8} | {'Rate':<10}")
        print("-" * 92)
        
        for comp, counts in shared_metrics.items():
            metrics_list = [
                ("Shared Intergenic (Same Base)", counts["Intergenic_Same_Base"], total_intergenic_size),
                ("Shared Synonymous (Same Base)", counts["Synonymous_Same_Base"], total_genic_size),
                ("Shared Nonsynonymous (Same Base, Same AA)", counts["Nonsynonymous_Same_Base_Same_AA"], total_genic_size),
                ("Shared Nonsynonymous (Diff Base, Same AA) *TRUE CONVERGENCE*", counts["Nonsynonymous_Diff_Base_Same_AA"], total_genic_size),
                ("Shared Nonsynonymous Stop Creations", counts["Stop_Creation"], total_genic_size),
                ("Shared Nonsynonymous Stop Losses", counts["Stop_Loss"], total_genic_size),
            ]
            for label, count, space in metrics_list:
                rate = count / space if space > 0 else 0.0
                sw.writerow([comp, label, count, f"{rate:.6f}"])
                print(f"{comp:<20} | {label:<45} | {count:<8} | {rate:.6f}")
            print("-" * 92)

def run_context_aware_simulation(genome_size, iterations=10000):
    """Runs context-aware permutation tests based explicitly on meaningful convergence metrics."""
    pairs = [
        ("LSDV", "SPPV"), ("LSDV", "SPPV_modern"), ("LSDV", "GTPV"), 
        ("GTPV", "SPPV"), ("GTPV", "SPPV_modern")
    ]
    sppv_anc = mutation_store.get("SPPV", {})
    results, raw_p_values = [], []

    for t1, t2 in pairs:
        comp_key = f"{t1} vs {t2}"
        muts1, muts2 = mutation_store.get(t1, {}), mutation_store.get(t2, {})
        
        obs_count = 0
        shared_positions = set(muts1.keys()) & set(muts2.keys())
        for pos in shared_positions:
            m1, m2 = muts1[pos], muts2[pos]
            if m1['mut_type'] == 'Incomplete_Codon_Context' or m2['mut_type'] == 'Incomplete_Codon_Context':
                continue
            if t1 == "SPPV_modern" and pos in sppv_anc and sppv_anc[pos]['derived_nt'] == m1['derived_nt']:
                continue
            if t2 == "SPPV_modern" and pos in sppv_anc and sppv_anc[pos]['derived_nt'] == m2['derived_nt']:
                continue
                
            if m1['mut_type'] in ['Intergenic', 'Synonymous']:
                if m1['derived_nt'] == m2['derived_nt']:
                    obs_count += 1
            else: # Nonsynonymous
                if m1['derived_aa'] == m2['derived_aa']:
                    obs_count += 1

        pool1_genic = [(pos, m['derived_nt']) for pos, m in muts1.items() if m['mut_type'] not in ['Intergenic', 'Incomplete_Codon_Context']]
        pool1_inter = [(pos, m['derived_nt']) for pos, m in muts1.items() if m['mut_type'] == 'Intergenic']
        pool2_genic = [(pos, m['derived_nt']) for pos, m in muts2.items() if m['mut_type'] not in ['Intergenic', 'Incomplete_Codon_Context']]
        pool2_inter = [(pos, m['derived_nt']) for pos, m in muts2.items() if m['mut_type'] == 'Intergenic']

        genic_bound = int(genome_size * 0.90)
        inter_bound = genome_size - genic_bound

        sim_counts = []
        for _ in range(iterations):
            r_pos1_g = random.sample(range(1, genic_bound + 1), len(pool1_genic))
            r_pos2_g = random.sample(range(1, genic_bound + 1), len(pool2_genic))
            sim1_g = dict(zip(r_pos1_g, [nt for _, nt in pool1_genic]))
            sim2_g = dict(zip(r_pos2_g, [nt for _, nt in pool2_genic]))
            
            r_pos1_i = random.sample(range(1, inter_bound + 1), len(pool1_inter))
            r_pos2_i = random.sample(range(1, inter_bound + 1), len(pool2_inter))
            sim1_i = dict(zip(r_pos1_i, [nt for _, nt in pool1_inter]))
            sim2_i = dict(zip(r_pos2_i, [nt for _, nt in pool2_inter]))
            
            overlap_g = sum(1 for p in (set(sim1_g.keys()) & set(sim2_g.keys())) if sim1_g[p] == sim2_g[p])
            overlap_i = sum(1 for p in (set(sim1_i.keys()) & set(sim2_i.keys())) if sim1_i[p] == sim2_i[p])
            sim_counts.append(overlap_g + overlap_i)

        expected = sum(sim_counts) / iterations
        p = sum(1 for c in sim_counts if c >= obs_count) / iterations if (obs_count > 0 or expected > 0) else 1.0
        raw_p_values.append(p)
        results.append({'comparison': comp_key, 'observed': obs_count, 'expected': expected, 'raw_p': p})

    m = len(raw_p_values)
    indexed_p = sorted(list(enumerate(raw_p_values)), key=lambda x: x[1])
    bh_p_values = [0.0] * m
    prev_adj = 1.0
    for rank, (orig_idx, p) in reversed(list(enumerate(indexed_p, start=1))):
        adj_p = min(p * m / rank, prev_adj, 1.0)
        bh_p_values[orig_idx] = adj_p
        prev_adj = adj_p

    print(f"\n{'Comparison':<20} | {'Observed':<10} | {'Expected (Rand)':<16} | {'Raw P-Value':<18} | {'BH Adj P-Value':<18}")
    print("-" * 92)
    for i, res in enumerate(results):
        print(f"{res['comparison']:<20} | {res['observed']:<10} | {res['expected']:<16.4f} | {str(res['raw_p']):<18} | {str(bh_p_values[i]):<18}")

if __name__ == "__main__":
    print("Reading reference tracking structural properties from GenBank...")
    genome_seq, genes_metadata, position_to_gene = parse_genbank_rich(GB_FILE)
    genome_size = len(genome_seq)
    print(f"Loaded reference successfully. Genome Size: {genome_size} bp. CDS count: {len(genes_metadata)}\n")

    tasks = [
        ("lsdv_ancestral3.aln", "LSDV", "lsdv_summary.csv"),
        ("sppv_ancestral3.aln", "SPPV", "sppv_summary.csv"),
        ("sppv_modern_ancestral3.aln", "SPPV_modern", "sppv_modern_summary.csv"),
        ("gtpv_ancestral3.aln", "GTPV", "gtpv_summary.csv")
    ]

    for aln, taxon, summary in tasks:
        process_lineage(aln, taxon, summary, genome_seq, position_to_gene, genes_metadata)
        
    write_shared_mutations_filtered(genome_seq, genes_metadata)
    print("\nRunning background context-preserving permutation engine simulations...")
    run_context_aware_simulation(genome_size)
