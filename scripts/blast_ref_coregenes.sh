#!/bin/bash
# This is a script for the sequence-trait correlation of the core gene sequences. Similar to the `main_script.sh` pipeline, core gene sequences are stored to a singel fasta file and then blasted against the LE database genomes. The only difference is that the core gene seqeucnces is retrieved by the `pull_seq_from_pangenomes.rmd` script, so this script will skip the genebank retrieval section.

# Set up logging
log_dir="./logs/"
mkdir -p "$log_dir"
log_file="$log_dir/blast_coregenes_to_refgeno_$(date +%Y%m%d_%H%M%S).log"
echo "Script started at: $(date)" | tee "$log_file"
echo "Log file: $log_file" | tee -a "$log_file"
echo "----------------------------------------" | tee -a "$log_file"

# Function to log messages
log_message() {
    echo "$1" | tee -a "$log_file"
}

## User input:
# Define the strain name and core gene fasta file.
strain_name="Microcystis_aeruginosa_BS13-10"
# Path to the reference genotype
path_to_reference="/Users/linyusheng/Grazing_Resist_Genes/BGCs/reference_core_genes/${strain_name}/"
# Path to the database of bacterial genomes, i.e. fasta files for predicted gene nucleotide sequences.
path_to_db_fasta_dir="/Users/linyusheng/Desktop/LE_annotation/bakta"

## Directories and files:
# Path to python scripts
path_to_python_scripts="/Users/linyusheng/Grazing_Resist_Genes/scripts/"


# Check if database directory exists and contains .ffn files
if [ ! -d "$path_to_db_fasta_dir" ]; then
    log_message "Error: Database directory not found: $path_to_db_fasta_dir"
    exit 1
fi

ffn_count=$(ls -1 $path_to_db_fasta_dir/*.ffn 2>/dev/null | wc -l)
if [ $ffn_count -eq 0 ]; then
    log_message "Error: No .ffn files found in database directory: $path_to_db_fasta_dir"
    exit 1
fi

log_message "Found $ffn_count .ffn files in database directory."

# Check if there are files ending in .fasta under `path_to_reference`
if [ ! -d "$path_to_reference" ]; then
    log_message "Error: Reference directory not found: $path_to_reference"
    exit 1
fi

fasta_count=$(ls -1 $path_to_reference/*.fasta 2>/dev/null | wc -l)
if [ $fasta_count -eq 0 ]; then
    log_message "Error: No .fasta files found in reference directory: $path_to_reference"
    exit 1
fi

# Start Script
log_message "========================================="
log_message "Processing strain: $strain_name"
log_message "========================================="

## Set paths for current BGC
reference_coreg_fasta=$(find "$path_to_reference/" -name "*.fasta" -type f | head -1)
output_blast_results_dir="$path_to_reference/blast_results"
output_summary_table="$output_blast_results_dir/${strain_name}_coregenes_targets_summary.txt"

## Create necessary directories if they don't exist
mkdir -p "$output_blast_results_dir"

## BLAST reference core genes to LE database

log_message "Starting BLAST searches for {$strain_name} against LE genome databases..."
for db_fasta in $path_to_db_fasta_dir/*.ffn; do
    if [ -f "$db_fasta" ]; then
        db_name=$(basename "$db_fasta" .ffn)
        log_message "Running BLAST for $target_BGC against database: $db_name"
            
        # Run BLAST search
        python "$path_to_python_scripts/blast_targets.py" \
            --mode blastn \
            --query_fasta "$reference_coreg_fasta" \
            --db_fasta "$db_fasta" \
            --output_dir "$output_blast_results_dir" \
            --min_query_coverage 90 2>&1 | tee -a "$log_file" # higher coverage criteria because its withinspecies comparison
    fi
done

log_message "BLAST searches completed."
    
# Combine the BLAST results into a summary table for this BGC
log_message "Combining BLAST results into a summary table."
if [ -n "$(ls -A "$output_blast_results_dir" 2>/dev/null)" ]; then
    ls -1 "$output_blast_results_dir"/*_blastn_table.txt > "$output_blast_results_dir/blast_results.txt" 2>/dev/null
    if [ -f "$output_blast_results_dir/blast_results.txt" ]; then
        python "$path_to_python_scripts/compile_blast_tables.py" \
            -f "$output_blast_results_dir/blast_results.txt" \
            -o "$output_summary_table" 2>&1 | tee -a "$log_file"
        log_message "Summary table created for $strain_name: $output_summary_table"
    else
        log_message "Warning: No BLAST result files found."
    fi
else
    log_message "Warning: No BLAST results found in $output_blast_results_dir"
fi

log_message "Completed core gene blasting for $strain_name"

# Log completion
log_message "========================================="
log_message "Script completed at: $(date)"
log_message "Log file saved as: $log_file"
    