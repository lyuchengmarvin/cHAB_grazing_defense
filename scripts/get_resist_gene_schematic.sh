#!/bin/bash
# This is the master script for creating gene schematic plot for grazing resistance genes
# Author: Yu-Cheng (Marvin) Lin
# Time: 2025-12-13

# Usage: bash get_resist_gene_schematic.sh
# Input:
# 1. The path to the folder containing genome FASTA and GFF files
# 2. A tab-delimited file with genome name and gene name per line
# Output:
# 1. GenBank files for the specified genes
# This can then be used for gene schematic plotting using clinker or geneviewer.

# Step 1: Define input and output paths
GENOME_DIR="/Users/linyusheng/Desktop/LE_annotation/bakta"  # Directory containing genome FASTA and GFF files
resistant_BGC="Cyanopeptolin" # Example resistant Biosynthetic Gene Cluster
resistant_reference_genoytpe="LE19-388" # Example resistant reference genotype
GENE_LIST_DIR="/Users/linyusheng/Grazing_Resist_Genes/BGCs/curated_bgc_genotype/$resistant_reference_genoytpe/$resistant_BGC/blast_results"  # Tab-delimited file with genome name and gene name per line
OUTPUT_DIR="/Users/linyusheng/Grazing_Resist_Genes/BGCs/resistant_gene_alignment/gene_schematic_gbks/$resistant_reference_genoytpe/$resistant_BGC"  # Directory to save output GenBank files

# Step 2: Create output directory if it doesn't exist
mkdir -p $OUTPUT_DIR

# Step 3: Loop through each genome and extract specified genes
GENE_TABLES=$(find $GENE_LIST_DIR -name "*_blastn_table.txt")
for GENE_TABLE in $GENE_TABLES; do
    GENOME_NAME=$(basename $GENE_TABLE | sed 's/_blastn_table.txt//')
    FASTA_FILE="$GENOME_DIR/${GENOME_NAME}.fna"
    GFF_FILE="$GENOME_DIR/${GENOME_NAME}.gff3"
    OUTPUT_GBK="$OUTPUT_DIR/${GENOME_NAME}_${resistant_BGC}.gbk"

    # First, remove the gene names that does not have a hit in the resistant biosynthetic gene cluster by removing those that has "Unmatched" in the second column

    awk -F'\t' '$2 != "Unmatched"' $GENE_TABLE > ${GENE_TABLE}_filtered.txt
    FILTERED_TABLE="${GENE_TABLE}_filtered.txt"
    
    # Run the Python script to extract genes and create GenBank file
    python3 ~/Grazing_Resist_Genes/scripts/resistant_gene2gbk.py --fasta $FASTA_FILE --gff $GFF_FILE --gene_list $FILTERED_TABLE -o $OUTPUT_GBK

    # Clean up temporary filtered table
    rm $FILTERED_TABLE
done
