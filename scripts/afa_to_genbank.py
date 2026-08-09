# This script reads an aligned fasta file and output each sequence to separate genbank files
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import os
import re
import sys

def sanitize_name(name):
    return re.sub(r'[^\w\.-]', '_', name)

def main():
    parser = argparse.ArgumentParser(description="Convert aligned FASTA to one GenBank file per sequence.")
    parser.add_argument("-i", "--input", required=True, help="Input aligned FASTA file")
    parser.add_argument("-o", "--outdir", default=".", help="Output directory for GenBank files")
    parser.add_argument("-p", "--prefix", default="", help="Optional prefix for output filenames (will be prepended to sequence id)")
    parser.add_argument("--organism", default="", help="Organism string to place in the GenBank record (optional)")
    parser.add_argument("--remove-gaps", action="store_true", help="Remove gap characters ('-' and '.') from sequences before writing")
    args = parser.parse_args()

    input_fasta = args.input
    outdir = args.outdir
    prefix = args.prefix
    organism = args.organism
    remove_gaps = args.remove_gaps

    if not os.path.isfile(input_fasta):
        print(f"Error: input file not found: {input_fasta}", file=sys.stderr)
        sys.exit(1)
    os.makedirs(outdir, exist_ok=True)

    for rec in SeqIO.parse(input_fasta, "fasta"):
        try:
            seq_str = str(rec.seq)
            if remove_gaps:
                seq_str = seq_str.replace("-", "").replace(".", "")
            seq = Seq(seq_str)

            seq_id = sanitize_name(rec.id)
            filename_id = f"{prefix}{seq_id}" if prefix else seq_id
            out_path = os.path.join(outdir, f"{filename_id}.gb")

            gb_rec = SeqRecord(seq,
                               id=seq_id,
                               name=rec.name if rec.name else seq_id,
                               description=rec.description if rec.description else "")
            # Basic annotations
            gb_rec.annotations["molecule_type"] = "DNA"
            if organism:
                gb_rec.annotations["organism"] = organism

            # minimal required fields for GenBank writer
            gb_rec.features = []

            SeqIO.write(gb_rec, out_path, "genbank")
            print(f"Wrote {out_path}")
        except Exception as e:
            print(f"Warning: failed to write record {rec.id}: {e}", file=sys.stderr)
            continue

if __name__ == "__main__":
    main()