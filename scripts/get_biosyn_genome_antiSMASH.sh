#!/bin/bash

# This script iterate over the Microcystis strain antiSMASH result folders and concatenate the genebank files into a single genebank file. This file will hold all the genomic regions prediced to code for biosynthetic gene clusters (BGCs) in specific Microcystis strain genome. Then, it runs the python script that converts the genebank file to fasta and gff files.
# Usage: ./get_biosyn_genome_antiSMASH.sh /path/to/antiSMASH/results /path/to/output/

# a list of antiSMASH result folders
ANTISMASH_RESULTS_DIR=$1
OUTPUT_GENEBANK_DIR=$2

# list the antiSMASH result folders
RESULT_FOLDERS=$(find "$ANTISMASH_RESULTS_DIR" -maxdepth 1 -type d -not -path "$ANTISMASH_RESULTS_DIR")

# Iterate over each result folder and concatenate the genebank files
for FOLDER in $RESULT_FOLDERS; do
    GENEBANK_FILES=$(find "$FOLDER" -maxdepth 1 -type f -name "*region001.gbk" -or -name "*region001.gbff")
    if [ -z "$GENEBANK_FILES" ]; then
        echo "No genebank files found in $FOLDER, skipping."
        continue
    else
        # concatenate multiple genebank files into a single genebank file
        OUTPUT_FILE="$OUTPUT_GENEBANK_DIR/$(basename "$FOLDER")_combined.gbk"
        echo "Concatenating genebank files in $FOLDER into $OUTPUT_FILE"
        cat $GENEBANK_FILES > "$OUTPUT_FILE"
        echo "Finished processing $FOLDER."
        # Run the python script to convert genebank to fasta and gff
        python3 scripts/convert_antiSMASH_gbk.py "$OUTPUT_FILE" -o "$OUTPUT_GENEBANK_DIR/$(basename "$FOLDER")"
    fi
done