#!/bin/bash
# This script is used to loop through all GeneBank Files under a folder to perform sequence extraction

# Set up logging
log_dir="./logs/"
mkdir -p "$log_dir"
log_file="$log_dir/extract_antiSMASH_gbk_seq_$(date +%Y%m%d_%H%M%S).log"
echo "Script started at: $(date)" | tee "$log_file"
echo "Log file: $log_file" | tee -a "$log_file"
echo "----------------------------------------" | tee -a "$log_file"

# Function to log messages
log_message() {
    echo "$1" | tee -a "$log_file"
}

## User input:
# Define the strain name and reference GenBank file.
antiSMASH_gbk_dir="/Users/linyusheng/Grazing_Resist_Genes/BGCs/annotate_ref_bgc/antiSMASH_results/BS13-10/"
# Define the name for the output FASTA file
antiSMASH_fasta_output="/Users/linyusheng/Grazing_Resist_Genes/BGCs/annotate_ref_bgc/antiSMASH_results/extracted_fna/BS13-10_antiSMASH_BGC_genes.fasta"
# Path to the gbk seq extraction python script
path_to_python_scripts="/Users/linyusheng/Grazing_Resist_Genes/scripts/"

## Check if there are any gbk files in the folder
if [ ! -d "$antiSMASH_gbk_dir" ]; then
    log_message "Error: gbk directory not found: $antiSMASH_gbk_dir"
    exit 1
fi

gbk_count=$(ls -1 $antiSMASH_gbk_dir/*.gbk 2>/dev/null | wc -l)
if [ $gbk_count -eq 0 ]; then
    log_message "Error: No .gbk files found in gbk directory: $antiSMASH_gbk_dir"
    exit 1
fi

log_message "Found $gbk_count .gbk files in gbk directory."

## Start processing each gbk file

for gbk_file in $antiSMASH_gbk_dir/*.gbk; do
    target_bgc_name=$(basename "$gbk_file" .gbk)
    output_ref_geno_fasta="${antiSMASH_gbk_dir}/${target_bgc_name}.fasta"

    log_message "Extracting gene sequences from GenBank file: $gbk_file"
    # Call the Python script for sequence extraction
    python "$path_to_python_scripts/get_gene_nt_from_gbk.py" \
        -g "$gbk_file" \
        -o "$output_ref_geno_fasta" 2>&1 | tee -a "$log_file"
done

## Concatenate all extracted fasta files into a single file
cat ${antiSMASH_gbk_dir}*.fasta > ${antiSMASH_fasta_output}
log_message "All gene sequences have been extracted and concatenated into: $antiSMASH_fasta_output"