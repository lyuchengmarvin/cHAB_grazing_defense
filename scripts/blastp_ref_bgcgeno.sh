#!/bin/bash
# This script will perform BLASTP searches against a database of bacterial genomes for multiple BGCs.

# Set up logging
log_dir="./logs/"
mkdir -p "$log_dir"
log_file="$log_dir/blast_BGCs_to_refgeno_$(date +%Y%m%d_%H%M%S).log"
echo "Script started at: $(date)" | tee "$log_file"
echo "Log file: $log_file" | tee -a "$log_file"
echo "----------------------------------------" | tee -a "$log_file"

# Function to log messages
log_message() {
    echo "$1" | tee -a "$log_file"
}

## User input:
# Define the strain name and reference GenBank file.
reference_name="LE19-84"
# Path to the reference genotype
path_to_reference="/Users/linyusheng/Grazing_Resist_Genes/BGCs/curated_bgc_genotype/${reference_name}"
# Path to the database of bacterial genomes, i.e. fasta files for predicted gene nucleotide sequences.
path_to_db_fasta_dir="/Users/linyusheng/Desktop/LE_annotation/bakta"
# Coverage threshold for BLAST search (e.g., 20 means 20% coverage of the query sequence)
coverage_threshold=40

## Directories and files:
# Path to python scripts
path_to_python_scripts="/Users/linyusheng/Grazing_Resist_Genes/scripts/"

# Check if database directory exists and contains .faa files
if [ ! -d "$path_to_db_fasta_dir" ]; then
    log_message "Error: Database directory not found: $path_to_db_fasta_dir"
    exit 1
fi

faa_count=$(ls -1 $path_to_db_fasta_dir/*.faa 2>/dev/null | wc -l)
if [ $faa_count -eq 0 ]; then
    log_message "Error: No .faa files found in database directory: $path_to_db_fasta_dir"
    exit 1
fi

log_message "Found $faa_count .faa files in database directory."

# Find all target BGC directories
bgc_dirs=$(find "$path_to_reference" -maxdepth 1 -type d -not -path "$path_to_reference" 2>/dev/null)
bgc_count=$(echo "$bgc_dirs" | wc -l)

if [ -z "$bgc_dirs" ] || [ $bgc_count -eq 0 ]; then
    log_message "Error: No BGC directories found in: $path_to_reference"
    exit 1
fi

log_message "Found $bgc_count BGC directories to process."

# Process each BGC directory
for bgc_dir in $bgc_dirs; do
    target_BGC=$(basename "$bgc_dir")
    if [ "$target_BGC" == "miscallaneous" ]; then
        log_message "Skipping directory: $target_BGC"
        continue
    fi
    
    log_message "========================================="
    log_message "Processing BGC: $target_BGC"
    log_message "========================================="
    
    # Set paths for current BGC
    output_blast_results_dir="$path_to_reference/$target_BGC/blast_results"
    output_summary_table="$output_blast_results_dir/${reference_name}_${target_BGC}_targets_summary.txt"
    
    # Create necessary directories if they don't exist
    log_message "Creating directories for $target_BGC..."
    mkdir -p "$output_blast_results_dir"

    # Verify that fasta file exists for this BGC
    fasta_count=$(ls -1 $path_to_reference/$target_BGC/*.faa 2>/dev/null | wc -l)
    if [ $fasta_count -eq 0 ]; then
        log_message "Error: No .faa files found in directory: $path_to_reference/$target_BGC"
        exit 1
    fi

    log_message "Found $fasta_count .faa files in directory: $path_to_reference/$target_BGC"
    # Concatenate all fasta files in the BGC directory to create a reference genotype fasta file
    if [ -f "$path_to_reference/$target_BGC/${reference_name}_${target_BGC}.faa" ]; then
        rm -f "$path_to_reference/$target_BGC/${reference_name}_${target_BGC}.faa"
        log_message "Removing already existed file: $path_to_reference/$target_BGC/${reference_name}_${target_BGC}.faa"
        log_message "Creating reference genotype fasta file for $target_BGC..."
    else
        log_message "Creating reference genotype fasta file for $target_BGC..."
    fi
    cat $path_to_reference/$target_BGC/*.faa > "$path_to_reference/$target_BGC/${reference_name}_${target_BGC}.faa"
    output_ref_geno_fasta="${reference_name}_${target_BGC}.faa"

    # Loop through all .faa files in the database directory and run BLAST searches
    log_message "Starting BLAST searches for $target_BGC against multiple databases..."
    for db_fasta in $path_to_db_fasta_dir/*.faa; do
        if [ -f "$db_fasta" ]; then
            db_name=$(basename "$db_fasta" .faa)
            log_message "Running BLAST for $target_BGC against database: $db_name"
            
            # Run BLAST search
            python "$path_to_python_scripts/blast_targets.py" \
                --mode blastp \
                --query_fasta "$path_to_reference/$target_BGC/$output_ref_geno_fasta" \
                --db_fasta "$db_fasta" \
                --output_dir "$output_blast_results_dir" \
                --min_query_coverage $coverage_threshold 2>&1 | tee -a "$log_file"
        fi
    done
    
    log_message "BLAST searches completed for $target_BGC."
    
    # Combine the BLAST results into a summary table for this BGC
    log_message "Combining BLAST results into summary table for $target_BGC..."
    if [ -n "$(ls -A "$output_blast_results_dir" 2>/dev/null)" ]; then
        ls -1 "$output_blast_results_dir"/*_blastp_table.txt > "$output_blast_results_dir/blast_results.txt" 2>/dev/null
        if [ -f "$output_blast_results_dir/blast_results.txt" ]; then
            python "$path_to_python_scripts/compile_blast_tables.py" \
                -f "$output_blast_results_dir/blast_results.txt" \
                -o "$output_summary_table" 2>&1 | tee -a "$log_file"
            log_message "Summary table created for $target_BGC: $output_summary_table"
        else
            log_message "Warning: No BLAST result files found for $target_BGC"
        fi
    else
        log_message "Warning: No BLAST results found in $output_blast_results_dir for $target_BGC"
    fi
    
    log_message "Completed processing BGC: $target_BGC"
done

# Log completion
log_message "========================================="
log_message "All BGCs processed successfully!"
log_message "Script completed at: $(date)"
log_message "Log file saved as: $log_file"
