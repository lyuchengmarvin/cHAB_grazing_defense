# This is a script to extract nucleotide sequences of genes from a genebank (.gbk) file
# Author: Yu-Cheng (Marvin) Lin
# Date: 2024-10-10

# Load pakages
from Bio import SeqIO
import glob
import argparse

# Define the function to put all gene names in a genebank file into a list
def get_locus_tag(gbk_file):
    """Extracts locus tags from a GeneBank file.

    Args:
        file_path (str): Path to the GeneBank file.
        
    Returns:
        list: A list of locus tags.
    """
    locus_tags = []
    
    # parse the GeneBank file
    for record in SeqIO.parse(gbk_file, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                if "locus_tag" in feature.qualifiers:
                    locus_tags.append(feature.qualifiers["locus_tag"][0])
                    
    return locus_tags

# Define the function to extract gene sequences
def extract_cds_sequence(file_path, locus_tag):
    """Extracts the nucleotide sequence of a CDS with a specific locus tag from a GeneBank file.
    
    Args:
        file_path (str): Path to the GeneBank file.
        locus_tag (str): The locus tag of the desired CDS feature.
        
    Returns:
        str: The nucleotide sequence of the specified CDS, or None if not found.
    """
    # Parse the GeneBank file
    for record in SeqIO.parse(file_path, "genbank"):
        for feature in record.features:
            if feature.type == "CDS":
                if "locus_tag" in feature.qualifiers and locus_tag in feature.qualifiers["locus_tag"]:
                    # Extract the location of the CDS feature
                    cds_location = feature.location
                    # Extract the nucleotide sequence at the specified location
                    cds_sequence = cds_location.extract(record).seq
                    return str(cds_sequence)
    
    # Return None if the locus tag or CDS is not found
    return None
# Define function that writes gene sequences to a fasta file
def write_fasta(locus_tag, sequence, output_file):
    """Writes a nucleotide sequence to a fasta file.
    
    Args:
        locus_tag (str): The locus tag of the gene.
        sequence (str): The nucleotide sequence of the gene.
        output_file (str): The output file path.
    """
    with open(output_file, "a") as file:
        file.write(f">{locus_tag}\n")
        file.write(f"{sequence}\n")
# Main function to extract gene sequences and write to fasta files
def main():
    """ The function takes in gene names (.csv) and a genebank file, and write the gene sequences to fasta files.
    If gene names are not specified, the function will extract all gene sequences in the genebank file.

    Args:
        gbk_file (str): Path to the genebank file.
        output_dir (str): Path to the output directory.
        gene_list (str): A list of gene names separated by commas
    """
    # Set up argument parser
    parser = argparse.ArgumentParser(description="Extract gene sequences from a GeneBank file and write to fasta files.")
    parser.add_argument("-g", "--gbk_file", type=str, help="Path to the GeneBank file.", required=True)
    parser.add_argument("-o","--output_dir", type=str, help="Path to the output directory.")
    parser.add_argument("-l","--locus_tags", type=str, help="Comma-separated list of locus tags.", default=None)

    # Parse arguments
    args = parser.parse_args()
    gbk_file = args.gbk_file
    output_dir = args.output_dir
    locus_tags = args.locus_tags

    # Check if the genebank file exists
    if not glob.glob(gbk_file):
        raise FileNotFoundError(f"The genebank file {gbk_file} does not exist.")
    # If the output directory exists, delete it
    if glob.glob(output_dir):
        raise FileExistsError(f"The output directory {output_dir} already exists.")
    # Check if the gene_list is specified
    if locus_tags:
        with open(locus_tags, "r") as file:
            locus_tags = file.read().strip().split(",")
    else:
        # Get all gene names in the genebank file
        locus_tags = get_locus_tag(gbk_file)

    # Extract gene sequences for each gene
    for gene_name in locus_tags:
        sequence = extract_cds_sequence(gbk_file, gene_name)
        if sequence:
            write_fasta(gene_name, sequence, output_dir)
        else:
            print(f"Gene sequence for {gene_name} not found.")
            # Main function to extract gene sequences

if __name__ == '__main__':
   main()
