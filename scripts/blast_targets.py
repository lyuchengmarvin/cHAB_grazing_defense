import subprocess
from collections import defaultdict
from Bio.Blast import NCBIXML
import os
import argparse

# Obtain index file for the nucleotide database
## Run makeblastdb outside of Python
def run_makeblastdb(dbtype, input_fasta, output_prefix):
    """Run makeblastdb to create a nucleotide BLAST database."""
    command = [
        "makeblastdb",
        "-in", input_fasta,
        "-dbtype", dbtype,
        "-out", output_prefix,
        "-title", output_prefix
    ]
    print("Running makeblastdb with command:", " ".join(command))
    subprocess.run(command, check=True)
    print(f"BLAST database has been created.")

## Run blastp outside of Python
def run_blastp(query, db, evalue=0.01):
    """Run BLASTp and save results to an XML file."""
    output = f"{db}.xml"
    command = [
        "blastp",
        "-query", query,
        "-db", db,
        "-evalue", str(evalue),
        "-outfmt", "5",  # Output format 5 corresponds to XML
        "-out", output
    ]
    print("Running blastp with command:", " ".join(command))
    subprocess.run(command, check=True)
    print(f"BLASTp completed.")

## Run blastn outside of Python
def run_blastn(query, db, evalue=0.01):
    """Run BLASTn and save results to an XML file."""
    output = f"{db}.xml"
    command = [
        "blastn",
        "-query", query,
        "-db", db,
        "-evalue", str(evalue),
        "-outfmt", "5",  # Output format 5 corresponds to XML
        "-out", output
    ]
    print("Running blastn with command:", " ".join(command))
    subprocess.run(command, check=True)
    print(f"BLASTn completed.")

# Parse the BLAST XML output
def parse_blast_xml(xml_file, min_query_coverage=90):
    """Parse BLAST XML output and extract relevant information.
        min_query_coverage: minimum query coverage to consider a hit valid (default: 90).
        This function iterates through each BLAST record and finds the best hit based on identity score and query coverage.
    """
    results = {}
    query_lengths = {}
    with open(xml_file) as f:
        blast_records = NCBIXML.parse(f)
        for blast_record in blast_records:
            query_name = blast_record.query
            query_length = blast_record.query_length
            query_lengths[query_name] = query_length  # Store query length
            best_hit = None
            best_identity_score = 0
            best_query_coverage = 0

            for alignment in blast_record.alignments:
                ref_name = alignment.hit_def.split()[0]
                for hsp in alignment.hsps:
                    identity_score = (hsp.identities / hsp.align_length) * 100
                    query_coverage = (hsp.align_length / query_length) * 100

                    if query_coverage >= min_query_coverage and identity_score > best_identity_score:
                        best_hit = ref_name
                        best_identity_score = identity_score
                        best_query_coverage = query_coverage

            if best_hit:
                results[query_name] = {
                    best_hit: {
                        'identity_score': best_identity_score,
                        'query_coverage': best_query_coverage,
                        'query_length': query_length  # Add query length here
                    }
                }
    return results
# Check unmatched results and add to the table
def check_unmatched_results(results, query_fasta, output_file):
    """Check for unmatched queries and add them to the table."""
    query_names = []
    with open(query_fasta) as qf:
        for line in qf:
            if line.startswith(">"):
                query_names.append(line.strip().lstrip(">"))

    if not results:
        with open(output_file, 'a') as f:
            for query_name in query_names:
                f.write(f"{query_name}\tUnmatched\t0.0\t0.0\t0\n")
        print("Unmatched queries have been added to the table.")
    else:
        with open(output_file, 'a') as f:
            for query_name in query_names:
                if query_name not in results:
                    f.write(f"{query_name}\tUnmatched\t0.0\t0.0\t0\n")
        print("Unmatched queries have been added to the table.")

# Write the results to a table
def write_results_to_table(results, output_file):
    """Write the BLAST results to a table format."""
    with open(output_file, 'w') as f:
        f.write("Query_Name\tReference_Name\tIdentity_Score\tQuery_Coverage\tQuery_Length\n")
        for query_name, ref_hits in results.items():
            for ref_name, data in ref_hits.items():
                f.write(f"{query_name}\t{ref_name}\t{data['identity_score']:.2f}\t{data['query_coverage']:.2f}\t{data['query_length']}\n")
    print(f"Results have been written to {output_file}")
def cleanup(db_prefix, output_xml):
    """Remove temporary files created by BLAST."""
    tempfiles = [f"{db_prefix}.ndb",f"{db_prefix}.nhr", f"{db_prefix}.nin", f"{db_prefix}.njs", f"{db_prefix}.not", f"{db_prefix}.nsq", f"{db_prefix}.ntf", f"{db_prefix}.nto", f"{db_prefix}.pdb",f"{db_prefix}.phr", f"{db_prefix}.pin", f"{db_prefix}.pjs", f"{db_prefix}.pot", f"{db_prefix}.psq", f"{db_prefix}.ptf", f"{db_prefix}.pto", output_xml]
    for file in tempfiles:
        if os.path.exists(file):
            os.remove(file)
def main():
    # Specifty variables from command line arguments
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Run BLASTn or BLASTp and process results.")
    parser.add_argument("-m","--mode", type=str, help="Mode of BLAST to run (blastn or blastp).", required=True)
    parser.add_argument("-q","--query_fasta", help="Path to the query FASTA file.")
    parser.add_argument("-d","--db_fasta", help="Path to the database FASTA file.")
    parser.add_argument("-c","--min_query_coverage", type=float, default=90, help="Minimum query coverage to consider (default: 90).")
    parser.add_argument("-e","--evalue", type=float, default=0.01, help="E-value threshold for BLAST (default: 0.01).")
    parser.add_argument("-p", "--output_prefix", type=str, help="Output prefix for BLAST database and results.")
    parser.add_argument("-o","--output_dir", help="Output directory for results (default: ./).", default="./")

    args = parser.parse_args()

    # Assign arguments to variables
    mode = args.mode
    query_fasta = args.query_fasta
    db_fasta = args.db_fasta
    min_query_coverage = args.min_query_coverage
    evalue = args.evalue
    output_dir = args.output_dir

    # Run the main function for each database
    if not args.output_prefix:
        args.output_prefix = db_fasta.split('/')[-1].split('.')[0]
    if mode == "blastp":
        output_table = f"{output_dir}/{args.output_prefix}_blastp_table.txt"
        # Build the BLAST database
        run_makeblastdb('prot', db_fasta, args.output_prefix)
        # Run BLASTp
        run_blastp(query=query_fasta, db=args.output_prefix, evalue=evalue)
    elif mode == "blastn":
        output_table = f"{output_dir}/{args.output_prefix}_blastn_table.txt"
        # Build the BLAST database
        run_makeblastdb('nucl', db_fasta, args.output_prefix)
        # Run BLASTn
        run_blastn(query=query_fasta, db=args.output_prefix, evalue=evalue)
    else:
        print("Invalid mode. Please specify blastn or blastp.")
    
    # Parse BLAST results
    out_xml = f"{args.output_prefix}.xml"
    results = parse_blast_xml(out_xml, min_query_coverage)
    # Write results to a table
    write_results_to_table(results, output_table)
    # Check for unmatched queries
    check_unmatched_results(results, query_fasta, output_table)
    # Clean up temporary files from makeblastdb and BLAST
    cleanup(args.output_prefix, out_xml)
if __name__ == "__main__":
    main()