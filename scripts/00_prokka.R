# conda install -c bioconda blast=2.2.31
# conda create -n prokka_env2 -c conda-forge -c bioconda prokka
# conda activate prokka_env2
# prokka --version  # Confirm it's working

# prokka 1.15.6
# Define the input FASTA files (in this case, just "x.fasta")
files <- list.files(pattern = "\\.fasta$")
# avoid 'combined'

# Loop through each file
for (file in files) {
  # Remove .fa extension for prefix
  prefix <- sub("\\.fasta$", "", file)
  
  # Define output directory name
  outdir <- paste0( prefix, "_PROKKA")
  
  # Create the output log path
  log_file <- file.path("OUT/", paste0(prefix, ".out"))
  
  # Create the Prokka command
  cmd <- sprintf(
    "prokka --kingdom Viruses --gffver 3 --usegenus --outdir %s --genus Capripoxvirus --prefix %s %s --force &> %s",
    outdir, prefix, file, log_file  )
  print(cmd)  
  # Run the command
#   system(cmd)
}
