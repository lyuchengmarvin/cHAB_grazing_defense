# This script extracts specified genes from a FASTA and GFF file
# and creates a GenBank file containing only those genes and their associated features.
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import argparse
import gffutils
import os

def read_gene_list(gene_file):
    """
    Reads a list of gene names from a tab-delimitedfile, one per line.
    Args:
        gene_file (str): A tab-delimited file with genome name, gene name per line.
    """
    gene_names = set()
    with open(gene_file, 'r') as f:
        # skip the header if present
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                gene_names.add(parts[1])  # assuming gene name is in the second column
    return list(gene_names)

def filter_and_create_genbank(fasta_file, gff_file, gene_list, output_gbk):
    """
    Filters FASTA/GFF by a list of gene names and creates a GenBank file.
    """
    
    print(f"Starting process to extract genes: {gene_list}")
    
    # 1. Load the GFF features and filter by gene name
    # We use gffutils to map gene names (often in the 'ID' or 'Name' attribute)
    try:
        # Create an in-memory database of the GFF file
        db = gffutils.create_db(
            gff_file, 
            dbfn=':memory:', 
            force=True, 
            keep_order=True
            # sort_by='start'  <-- Removed this line
        )
    except Exception as e:
        print(f"Error parsing GFF file: {e}")
        return

    filtered_features = []
    relevant_seq_ids = set()
    
    # Identify and collect all features (gene, CDS, mRNA) associated with the target genes
    for gene_name in gene_list:
        # Check if the gene itself exists (searching by ID or Name attribute)
        try:
            gene_feature = db[gene_name]
            relevant_seq_ids.add(gene_feature.seqid)
            
            # Find all children features (CDS, mRNA) related to this gene
            children = list(db.children(gene_feature, featuretype=None, order_by='start'))
            
            # Add the gene feature itself and its children
            filtered_features.extend([gene_feature] + children)
        
        except gffutils.FeatureNotFoundError:
            print(f"Warning: Gene '{gene_name}' not found in GFF file. Skipping.")
        except Exception as e:
            print(f"Error processing gene '{gene_name}': {e}")
            
    if not filtered_features:
        print("No features found for the specified gene list. Exiting.")
        return

    print(f"Found {len(filtered_features)} related features on sequence IDs: {relevant_seq_ids}")
    
    # 2. Extract sequences for the relevant contigs/scaffolds
    relevant_sequences = {}
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in relevant_seq_ids:
                relevant_sequences[record.id] = record
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        return

    # 3. Combine filtered features and sequence into a new GenBank file
    
    # Process one sequence ID at a time (e.g., one contig)
    for seq_id in relevant_seq_ids:
        
        sequence_record = relevant_sequences.get(seq_id)
        if not sequence_record:
            print(f"Sequence ID {seq_id} found in GFF but not in FASTA. Skipping.")
            continue
            
        # Create a new SeqRecord for the GenBank output
        genbank_record = SeqRecord(
            sequence_record.seq,
            id=sequence_record.id,
            name=sequence_record.id,
            description=f"Extracted features for: {', '.join(gene_list)}",
            annotations={"molecule_type": "DNA"}
        )
        
        # Add the filtered features
        for feature in [f for f in filtered_features if f.seqid == seq_id]:
            
            start = feature.start - 1  # GFF 1-based to Biopython 0-based
            end = feature.end
            strand = 1 if feature.strand == '+' else -1 if feature.strand == '-' else 0

            location = FeatureLocation(start, end, strand=strand)
            
            # Use 'attributes' for qualifiers
            qualifiers = dict(feature.attributes)
            
            seq_feature = SeqFeature(
                location=location,
                type=feature.featuretype, 
                qualifiers=qualifiers
            )
            genbank_record.features.append(seq_feature)

        # 4. Write the final SeqRecord object to GenBank format
        # If there are multiple contigs, we append them to the same output file
        mode = "a" if os.path.exists(output_gbk) else "w"
        try:
            with open(output_gbk, mode) as handle:
                # Use a list to write multiple records if needed
                SeqIO.write([genbank_record], handle, "genbank")
            print(f"Successfully added features from {seq_id} to {output_gbk}")
        except Exception as e:
            print(f"Error writing GenBank file: {e}")


def main():
    """
    Main function to parse command-line arguments and run the conversion.
    """
    parser = argparse.ArgumentParser(
        description="Extract specified genes from FASTA/GFF and create a GenBank file."
    )
    parser.add_argument("--fasta", help="Input FASTA file")
    parser.add_argument("--gff", help="Input GFF file")
    parser.add_argument(
        "--gene_list", 
        help="A tab-delimited file with genome name and gene name per line."
    )
    parser.add_argument(
        "-o", 
        "--output", 
        default="filtered_genes.gbk", 
        help="Output GenBank file name (default: filtered_genes.gbk)"
    )
    args = parser.parse_args()
    
    # Read gene names from the provided file
    gene_names = read_gene_list(args.gene_list)
    # Filter and create GenBank file
    filter_and_create_genbank(args.fasta, args.gff, gene_names, args.output)

if __name__ == "__main__":
    main()