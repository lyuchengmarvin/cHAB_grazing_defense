# This script will output a knship matrix for population corrrection in pyseer
# Author: Yucheng Marvin Lin
# Time: 2025-12-03

## This script run pyseer externally
import argparse
import os
import dendropy
import pandas as pd

## User input:
# phylogeny_dir = "path_to/phylogeny_file.nwk"
# output_dir = "path_to/output_directory/phylogeny_distances.tsv"

## Phylogeny-based method: Extract distance matrix from a phylogeny
def phylogeny2distmatrix(phylogeny_dir, output_dir, epsilon=1e-8):
    """
    Extract a distance matrix from a phylogeny. Since MugGWAS will implement the FaST-LMM model in pyseer, the similarities based on the shared branch length between each pair's MRCA and the root will be calculated.
    Args:
        phylogeny_dir (str): Path to the phylogeny file.
        output_dir (str): Path to the output directory.
        epsilon (float): Small value to add to every matrix cell to avoid zeros.
    """
    # Check if the phylogeny file exists
    if not os.path.exists(phylogeny_dir):
        raise FileNotFoundError(f"Phylogeny file {phylogeny_dir} does not exist.")
    
    # Load the phylogeny
    tree = dendropy.Tree.get(path=phylogeny_dir, schema="newick", preserve_underscores=True)
    
    dist_matrix = {}

    # Extract distance matrix
    pdm = tree.phylogenetic_distance_matrix()

    for idx1, taxon1 in enumerate(tree.taxon_namespace):
        dist_matrix[taxon1.label] = {}
        for taxon2 in tree.taxon_namespace:
            if taxon2.label not in dist_matrix[taxon1.label].keys():
                mrca = pdm.mrca(taxon1, taxon2)
                dist_matrix[taxon1.label][taxon2.label] = mrca.distance_from_root()

    # Convert to DataFrame, reindex, and save to tsv
    dist_df = pd.DataFrame(dist_matrix)
    dist_df = dist_df.reindex(dist_df.columns)

    # ensure numeric and add small epsilon to avoid exact zeros
    dist_df = dist_df.astype(float) + float(epsilon)

    dist_df.to_csv(output_dir, sep="\t", index=True)

def main():
    """
    Main function to run the phylogeny to distance matrix conversion.
    """
    parser = argparse.ArgumentParser(
        description="Extract a Kinship matrix from a phylogeny for population structure estimation."
    )
    parser.add_argument(
        "-nwk", "--phylogeny-file",
        help="The input phylogeny file in Newick format."
    )
    parser.add_argument(
        "-o","--output-kinship-tsv",
        help="Path to the output TSV file for the distance matrix."
    )
    parser.add_argument(
        "--epsilon", type=float, default=1e-8,
        help="Small value to add to every matrix cell to avoid zeros (default: 1e-8)."
    )
    args = parser.parse_args()
    # Example usage
    phylogeny_dir = args.phylogeny_file
    output_dir = args.output_kinship_tsv
    phylogeny2distmatrix(phylogeny_dir, output_dir, args.epsilon) 

if __name__ == "__main__":
    main()