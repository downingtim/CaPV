#!/bin/bash

# Total consolidated list of unique NCBI accessions (completely cleaned)
accessions=( "PV877838"
   # "NC_002642" "NC_005179" "NC_035460" "KX894508" 
   # "NC_003389" "NC_006966" "NC_001132" "NC_004002" 
   # "NC_016924" "NC_001266" "NC_004003" "NC_032111"
   # "AY077835" "AY077836" "KC951854" "KX576657" 
   # "MH381810" "MN072620" "MN072621" "MN072622" 
   # "MN072623" "MN072624" "MN072625" "MW020570" 
   # "PV167794" "AY077832" "AY077833" "AY077834" 
   # "MN072626" "MN072627" "MN072628" "MN072629" 
   # "MN072630" "MN072631" "MT137384" "MW020571" 
   # "MW167070" "MW167071" "ON961655" "ON961656" 
   # "ON961657" "OQ434235" "OQ434236" "OQ434237" 
   # "OQ434238" "OQ434239" "OR239060" "PP886236" 
   # "PP886237" "PP886238" "PP886239" "PQ014465" 
   # "PV167793" "PV434148"
#    "AF325528.1" "AF409137.1" 
#    "KX683219.1" "KX894508.1" "KY702007.1" "KY829023.3" 
#    "MH893760.2" "MN072619.1" "MN642592.1" "MN995838.1" 
#    "MT643825.1" "MW631933.1" "MW656253.1" "MW699032.1" 
#    "NC_003027.1" "OK318001.1" "OP297402.1" "OP688128.1" 
#    "OP688129.1" "OQ588787.1" "OR134832.1" "OR134833.1" 
#    "OR134835.1" "OR134836.1" "OR134837.1" "OR134838.1" 
#    "OR134839.1" "OR134840.1" "OR134841.1" "OR134843.1" 
#    "OR134844.1" "OR134845.1" "OR134846.1" "OR134847.1" 
#    "OR134849.1" "OR393169.1" "OR393170.1" "OR393171.1" 
#    "OR393172.1" "OR393173.1" "OR393175.1" "OR393176.1" 
#    "OR393177.1" "OR393178.1" "OR520147.1"
)

echo "Starting genome downloads for ${#accessions[@]} unique records..."

for acc in "${accessions[@]}"; do
    # Name the output file using the clean accession
    output_file="${acc}.fasta"
    
    # Skip download if the file is already sitting in the folder
    if [[ ! -f "$output_file" ]]; then
        echo "Downloading genome for: $acc..."
        
        # Pull down the complete nucleotide sequence 
        curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${acc}&rettype=fasta&retmode=text" > "$output_file"
        
        # Polite 1-second pause so NCBI servers don't throttle you
        sleep 1
    else
        echo "File already exists, skipping: $output_file"
    fi
done

echo "All genome downloads complete!"
