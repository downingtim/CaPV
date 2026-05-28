#!/bin/bash

root_dir="SEQS4"
ref_file="KX894508.gb"
final_output="total_proteome_deerswinepox.fasta"

# Clear any pre-existing run file
> "$final_output"

echo "Step 1: Processing Prokka sample folders..."
echo "------------------------------------------"

for folder in "$root_dir"/*_PROKKA; do
    [ -d "$folder" ] || continue
    
    sample_name=$(basename "$folder" "_PROKKA")
    faa_file="${folder}/${sample_name}.faa"
    
    if [[ -f "$faa_file" ]]; then
        echo "Appending: $sample_name"
        echo ">${sample_name}" >> "$final_output"
        grep -v "^>" "$faa_file" | tr -d '\n\r ' >> "$final_output"
        echo "" >> "$final_output" 
    fi
done

echo "------------------------------------------"
echo "Step 2: Processing reference file cleanly (${ref_file})..."
echo "------------------------------------------"

if [[ -f "$ref_file" ]]; then
    ref_name=$(basename "$ref_file" ".gb")
    echo "Extracting and appending reference via text parsing: $ref_name"
    echo ">${ref_name}_concatenated_proteome" >> "$final_output"
    
    # 1. Use sed to capture only the text trapped between /translation=" and the closing "
    # 2. Use tr to delete all spaces, line breaks, quotes, and metadata tags
    sed -n '/\/translation="/,/"/p' "$ref_file" | \
        sed 's/.*\/translation="//' | \
        tr -d '[:space:]"' >> "$final_output"
        
    echo "" >> "$final_output" # Final trailing newline
    echo "Successfully completed master matrix build."
else
    echo "[ ERROR ] Reference file '$ref_file' not found!"
fi
