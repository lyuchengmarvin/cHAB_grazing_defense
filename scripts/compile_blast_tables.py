import os
import argparse

def query_dict(strain_names):
    """
    Create a dictionary that store the matches for each query name.
    
    Args:
    strain_names (list): a list of strain names.
    """
    dictionary = {}
    for strain in strain_names:
        dictionary[strain] = {}

    return dictionary

def add_match(dictionary, strain_name, match_dict):
    """
    Add a match to the query dictionary strain-by-strain.
    
    Args:
    dictionary (dict): the dictionary to update.
    strain_name (str): to call the strain dictionary that has all query names.
    match_dict (dict): the match dictionary containing reference_name identity score, and query coverage to add.
    """
    for query_name, match in match_dict.items():
        dictionary[strain_name][query_name] = match


def process_file(filename, dictionary):
    """ 
    Process a single tab delimited file and update the dictionary.

    Args:
    filename (str): the name of the file to process.
    dictionary (dict): the dictionary to update.
    """
    strain_name = os.path.splitext(os.path.basename(filename))[0]
    strain_name = strain_name.replace('_blastn_table', '')

    with open(filename, 'r') as file:
        headers = file.readline().strip().split('\t')

        for line_num, line in enumerate(file, start=2):  # Start at 2 since we skip header
            fields = line.strip().split('\t')
            
            # Skip empty lines
            if not line.strip():
                continue
                
            # Check if line has enough fields
            if len(fields) < 5:
                print(f"Warning: Line {line_num} in {filename} has only {len(fields)} fields, skipping: {line.strip()}")
                continue
                
            try:
                query_name = fields[0]
                reference_name = fields[1]
                identity_score = float(fields[2])
                query_coverage = float(fields[3])
                query_length = int(fields[4])
                match_dict = {query_name: {
                    'Reference_name': reference_name,
                    'Identity_score': identity_score,
                    'Query_coverage': query_coverage,
                    'Query_length': query_length
                }}
                add_match(dictionary, strain_name, match_dict)
            except (ValueError, IndexError) as e:
                print(f"Warning: Error processing line {line_num} in {filename}: {e}")
                print(f"Line content: {line.strip()}")
                continue

def output_combined_table(dictionary, output_filename):
    """
    Output the combined table from the dictionary. Each row should contain query name, reference name, identity score, query coverage, and query length for each strain.

    Args:
    dictionary (dict): the dictionary to output.
    filenames (list): the list of filenames to use as column headers.
    """
    with open(output_filename, 'w') as file:
        # Print header
        file.write("Strain_Name\tQuery_Name\tReference_Name\tIdentity_Score\tQuery_Coverage\tQuery_Length\n")
        # Print rows
        for strain_name, match_dict in dictionary.items():
            for query_name, match in match_dict.items():
                file.write(f"{strain_name}\t{query_name}\t{match['Reference_name']}\t{match['Identity_score']}\t{match['Query_coverage']}\t{match['Query_length']}\n")

def main():
    # Define command line options
    parser = argparse.ArgumentParser(description='Combine BLAST results for multiple strains.')
    parser.add_argument('-f','--blast_result_filenames', help='A list of BLAST result files to process.')
    parser.add_argument('-o','--output_filename', help='Output file name.')
    args = parser.parse_args()
    blast_files = args.blast_result_filenames
    output_file = args.output_filename

    # Open the blast result files
    with open(blast_files,'r') as file:
        filenames = file.readlines()
        filenames = [filename.strip() for filename in filenames]

    # strip off the file extension '_blastn_table.txt' and path to get the strain names
    strain_names = [os.path.splitext(os.path.basename(filename))[0] for filename in filenames]
    strain_clean_names = [strain_name.replace('_blastn_table', '') for strain_name in strain_names]
    print(f"Strain names: {strain_clean_names}")

    # Loop over files and process each one
   
    dictionary = query_dict(strain_clean_names)
    for filename in filenames:
        process_file(filename, dictionary)
    output_combined_table(dictionary, output_file)

if __name__ == "__main__":
    main()