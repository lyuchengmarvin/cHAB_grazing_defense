# This script find pattern in headers and extracts that particular sequence to a new fasta file
# Author: Lin, Yu-Cheng (Marvin)
# Date: 2025-08-07
import sys
from Bio import SeqIO
import re
import argparse

# Function to put target gene IDs into a list
def put_target_gene_ids(gene_id_txt):
    """
    Ideal input format is a text file with gene IDs line by line or comma-separated.
    This function reads the file and returns a list of gene IDs.
    """
    with open(gene_id_txt) as handle:
        read = handle.readlines()
        handle.close()
        # Check if the file is empty
        if not read:
            print('The input file is empty. Please provide a valid file with gene IDs.')
            sys.exit(1)
        # Remove whitespace characters
        read = [line.strip() for line in read if line.strip()]
        # If the file has one line, split it into a list
        if len(read) == 1:
            if ',' in read[0]:
                query_gene_id = read[0].split(',')
            elif '\t' in read[0]:
                query_gene_id = read[0].split('\t')
        else:
            # If each line is a gene ID, return it as a list
            query_gene_id = [line for line in read if line]
            # If the file has multiple lines, yet still each line has separators, its format is not as expected
            if query_gene_id and (',' in query_gene_id[0] or '\t' in query_gene_id[0]):
                print('The input file format is not as expected. Please provide a valid file with gene IDs.')
                sys.exit(1)
    return query_gene_id

# Function to extract sequences from fasta based on gene IDs
def extract_sequences_from_fasta(gene_ids, fasta_file, output_file):
    """
    This function finds sequences in a fasta file based on the provided gene IDs. 
    It should treat the gene IDs as patterns to match against the fasta headers.
    This will return a new fasta file with the sequences that match the gene IDs.
    """
    genes_from_fasta = SeqIO.index(fasta_file, "fasta")
    with open(output_file, 'w') as output_handle:
        gn = 0
        for gene in gene_ids:
            found = False
            pattern = re.compile(gene)
            for record_id in genes_from_fasta:
                if pattern.search(record_id):
                    SeqIO.write(genes_from_fasta[record_id], output_handle, "fasta")
                    found = True
                    gn += 1
            if not found:
                print(f"Gene ID pattern '{gene}' not found in fasta file.")
    output_handle.close()
    print(f"A total of {gn} sequences were extracted from the input fasta with {len(genes_from_fasta)} sequences.")

# Main function
def main():
    """
    Main function to handle command line arguments and execute the sequence extraction.
    Args:
        gene_id_txt (str): Path to the text file containing gene IDs. It should be a text file with gene IDs line by line or comma-separated.
        fasta (str): Path to the input fasta file.
        output_fasta (str): Path to the output fasta file where extracted sequences will be saved.
    """
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Extract sequences from a fasta file with partial/complete match to the provided gene IDs.")
    parser.add_argument("-g", "--gene_id_txt", type=str, help="Path to the text file containing gene IDs.", required=True)
    parser.add_argument("-f", "--input_fasta", type=str, help="Path to the input fasta file.", required=True)
    parser.add_argument("-o", "--output_fasta", type=str, help="Path to the output fasta file.", required=True)
    args = parser.parse_args()
    
    gene_id_txt = args.gene_id_txt
    fasta = args.input_fasta
    output_fasta = args.output_fasta

    query_gene_id = put_target_gene_ids(gene_id_txt)

    extract_sequences_from_fasta(query_gene_id, fasta, output_fasta)
    print(f"Sequences were output to {output_fasta}.")
if __name__ == '__main__':
   main()