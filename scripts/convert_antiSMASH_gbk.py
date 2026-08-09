import re
import os
import argparse
from io import StringIO
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, CompoundLocation
import collections

def genbank_to_fasta(genbank_data):
    """
    Converts GenBank format data into FASTA format.

    Uses Biopython's SeqIO.write to correctly format the sequence record(s)
    into a FASTA string, including the definition line.

    Args:
        genbank_data (str): A string containing the content of a GenBank file.

    Returns:
        str or None: The sequence(s) formatted in FASTA, or None if no record is found.
    """
    handle_in = StringIO(genbank_data)
    handle_out = StringIO()
    
    # Parse the GenBank records
    records = list(SeqIO.parse(handle_in, "genbank"))
    
    if not records:
        return None

    # Use SeqIO.write to format the records into FASTA
    SeqIO.write(records, handle_out, "fasta")
    return handle_out.getvalue()


def genbank_to_gff3(genbank_data):
    """
    Converts GenBank format data features into GFF3 format annotations,
    handling multiple loci (records) within the input data.

    This function iterates through all records, extracts key features, and converts
    them into the tab-separated GFF3 format, correctly handling the 0-to-1 indexing.
    
    Args:
        genbank_data (str): A string containing the content of a GenBank file.

    Returns:
        str: A string containing the features formatted in GFF3.
    """
    gff3_output = ["##gff-version 3"]
    
    # Use StringIO to treat the string as a file handle for SeqIO.parse
    handle = StringIO(genbank_data)
    
    # Parse ALL GenBank records
    records = list(SeqIO.parse(handle, "genbank"))

    if not records:
        return "##gff-version 3\n# Error: No GenBank records found."
        
    for record in records:
        seqid = record.id
        
        # Add the sequence-region directive for the current record
        gff3_output.append(f"##sequence-region {seqid} 1 {len(record.seq)}")
        
        # Iterate over all features in the record
        feature_count = 0
        for feature in record.features:
            # Skip features without a defined location
            if not feature.location:
                continue

            # Standard GFF3 fields
            source = "Biopython"
            ftype = feature.type
            score = "."
            phase = "."

            # Strand conversion: 1 -> '+', -1 -> '-', 0 -> '.'
            strand = {1: '+', -1: '-', 0: '.'}.get(feature.location.strand, '.')

            # Attributes (Column 9)
            attributes = collections.OrderedDict()
            
            # Use 'locus_tag' or 'gene' as the primary Name/ID
            name = feature.qualifiers.get('locus_tag', [None])[0] or feature.qualifiers.get('gene', [None])[0]
            if name:
                # Create a unique ID for the feature (unique across the whole file)
                feature_count += 1
                feature_id = f"{seqid}_{ftype}_{feature_count}"
                attributes['ID'] = feature_id
                attributes['Name'] = name
            
            # Add product and other standard qualifiers
            product = feature.qualifiers.get('product', [None])[0]
            if product:
                # GFF3 attributes cannot contain unescaped semicolons, so we replace them
                attributes['product'] = product.replace(';', ',') 

            # Convert attributes dictionary to the GFF3 attribute string
            attr_list = [f"{k}={v}" for k, v in attributes.items()]
            attr_string = ";".join(attr_list)

            # Handle GenBank location complexity (e.g., 'join' locations)
            # A GenBank feature can correspond to multiple GFF3 lines (e.g., a CDS 
            # made of multiple exons/parts).
            parts = []
            if isinstance(feature.location, FeatureLocation):
                parts.append(feature.location)
            elif isinstance(feature.location, CompoundLocation):
                parts.extend(feature.location.parts)
            
            for part_location in parts:
                # Biopython locations are 0-indexed (start is excluded, end is included)
                # GFF3 is 1-indexed (start and end are included)
                start_1_indexed = int(part_location.start) + 1
                end_1_indexed = int(part_location.end)

                # Ensure we skip if start/end is invalid after conversion
                if start_1_indexed > end_1_indexed:
                    continue

                gff_line = [
                    seqid,
                    source,
                    ftype,
                    str(start_1_indexed),
                    str(end_1_indexed),
                    score,
                    strand,
                    phase,
                    attr_string
                ]
                gff3_output.append("\t".join(gff_line))

    return "\n".join(gff3_output)


def main():
    """
    Main function to parse command-line arguments and run the conversion.
    """
    parser = argparse.ArgumentParser(
        description="Convert a GenBank file to FASTA sequence and GFF3 annotation files."
    )
    parser.add_argument(
        "input_file",
        help="Path to the input GenBank (.gb or .gbk) file."
    )
    parser.add_argument(
        "-o", "--output-prefix",
        help="Optional prefix for output files. If not provided, the input filename (without extension) will be used.",
        default=None
    )
    args = parser.parse_args()

    # 1. Read input file
    try:
        with open(args.input_file, "r") as f:
            genbank_content = f.read()
    except FileNotFoundError:
        print(f"Error: Input file not found at '{args.input_file}'")
        return
    except Exception as e:
        print(f"Error reading input file: {e}")
        return

    # 2. Determine output prefix and paths
    if args.output_prefix:
        output_prefix = args.output_prefix
    else:
        # Use filename without extension as default prefix
        base_name = os.path.basename(args.input_file)
        output_prefix = os.path.splitext(base_name)[0]

    fasta_path = f"{output_prefix}.fa"
    gff3_path = f"{output_prefix}.gff3"

    print(f"Processing '{args.input_file}'...")
    
    # 3. Perform conversions
    fasta_result = None
    gff3_result = None
    try:
        fasta_result = genbank_to_fasta(genbank_content)
        gff3_result = genbank_to_gff3(genbank_content)
    except Exception as e:
        print(f"Error during Biopython parsing or conversion: {e}")
        return

    # 4. Write outputs
    if fasta_result:
        try:
            with open(fasta_path, "w") as f:
                f.write(fasta_result)
            print(f"✅ Successfully wrote FASTA sequence to '{fasta_path}'")
        except Exception as e:
            print(f"Error writing FASTA file: {e}")
    else:
        print("⚠️ Warning: Could not generate FASTA output (no records found).")

    if gff3_result and "Error" not in gff3_result:
        try:
            with open(gff3_path, "w") as f:
                f.write(gff3_result)
            print(f"✅ Successfully wrote GFF3 annotations to '{gff3_path}'")
        except Exception as e:
            print(f"Error writing GFF3 file: {e}")
    elif gff3_result:
        # This handles the case where genbank_to_gff3 returns an error string
        print(f"⚠️ Warning: Could not generate GFF3 output: {gff3_result.splitlines()[-1]}")
    else:
        print("⚠️ Warning: Could not generate GFF3 output.")


if __name__ == "__main__":
    main()