# This python script filter the presence absence table based on user-defined minor allele frequency
# Author: Yu-Cheng (Marvin) Lin
# Date: 2025-12-10

import os
import argparse
import pandas as pd


def cal_af (presence_absence_list):
    """ Calculate the allele frequency from presence absence list
    Args:
        presence_absence_list (list): A list of 0s and 1s representing absence and presence of a gene.
    Returns:
        float: The allele frequency of the gene.
    """
    total = len(presence_absence_list)
    count_ones = sum(presence_absence_list)
    af = count_ones / total
    return af

def filter_by_maf(pres_abs_df, maf_threshold):
    """
    Filter the presence-absence DataFrame based on minor allele frequency threshold.
    Args:
        pres_abs_df (pd.DataFrame): DataFrame containing presence-absence data.
        maf_threshold (float): Minor allele frequency threshold for filtering.
    Returns:
        pd.DataFrame: Filtered DataFrame with genes meeting the MAF criteria.
    """
    # Calculate allele frequencies for each gene (row)
    af_series = pres_abs_df.apply(cal_af, axis=1)
    
    # Determine which genes meet the MAF criteria
    genes_to_keep = af_series[(af_series >= maf_threshold)].index
    
    # Filter the DataFrame to keep only the selected genes
    filtered_df = pres_abs_df.loc[genes_to_keep]

    # Count and report the number of genes before and after filtering
    original_gene_count = pres_abs_df.shape[0]
    filtered_gene_count = filtered_df.shape[0]
    print(f"Number of genes before filtering: {original_gene_count}")
    print(f"Number of genes after filtering: {filtered_gene_count}")
    
    return filtered_df

def main():
    """
    Main function to filter presence-absence table based on minor allele frequency.
    Args: 
        pres_abs.tsv: Command line input presence-absence table file path
        maf: Command line input minor allele frequency threshold
    returns:
        Filtered presence-absence table saved as 'filtered_pres_abs.tsv'.
    """
    parser = argparse.ArgumentParser(
        description="Filter pyseer-outputed presence-absence table based on minor allele frequency."
    )
    parser.add_argument(
        "pres_abs",
        type=str,
        help="Tab-delimited presence-absence table in numeric format(0s and 1s)."
    )
    parser.add_argument(
        "-maf", "--minor_allele_frequency",
        type=float,
        default=0.01,
        help="Minor allele frequency threshold for filtering."
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default="filtered_pres_abs.tsv",
        help="Output file path for the filtered presence-absence table."
    )
    args = parser.parse_args()

    # 1. Read presence-absence table
    try:
        pres_abs_df = pd.read_csv(args.pres_abs, sep = "\t", index_col=0)
    except FileNotFoundError:
        print(f"Error: Presence-absence file not found at '{args.pres_abs}'")
        return
    except Exception as e:
        print(f"Error reading presence-absence file: {e}")
        return
    
    # 2. Filter by minor allele frequency
    filtered_df = filter_by_maf(pres_abs_df, args.minor_allele_frequency)

    # 3. Save filtered DataFrame to output file
    try:
        filtered_df.to_csv(args.output, sep="\t")
        print(f"Filtered presence-absence table saved to '{args.output}'")
    except Exception as e:
        print(f"Error saving filtered presence-absence table: {e}")
        return
    
if __name__ == "__main__":
    main()