#!/bin/bash

# Define the final output filename
output_file="LSDV_combined_genomes.fasta"
temp_file="combined.tmp"

# Clear out any pre-existing temp file
> "$temp_file"

echo "Merging and modifying FASTA files..."

# Loop through all .fasta files in the current directory
for file in *.fasta; do
    # Skip the final output file if it already exists in the folder
    if [[ "$file" == "$output_file" ]]; then
        continue
    fi

    echo "Processing: $file"
    
    # Use sed to replace the start of headers (">") with ">LSDV_"
    sed 's/^>/>LSDV_/' "$file" >> "$temp_file"
done

# Move the temporary file to the final output name
mv "$temp_file" "$output_file"

echo "Done! All sequences combined into: $output_file"
