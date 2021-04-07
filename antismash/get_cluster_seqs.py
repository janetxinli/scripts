#!/usr/bin/env python3
"""
Pull MIBiG cluster nucleotide sequences from MIBiG GenBank file
@author: janetxinli
"""

import sys
import argparse
from Bio import SeqIO

def get_cluster_record(gb, id):
    """
    Get specified cluster record from GenBank file from its ID.
    Returns None if no record is present with the name id in the
    GenBank file.
    """
    full = SeqIO.parse(gb, "gb")
    for record in full:
        if record.name == id:
            return record

    return None  # Cluster not found

def print_seqs(record, outfile):
    """Extract and print feature sequences from GenBank record."""
    if outfile is None:  # Print to stdout if no output file is provided
        outfh = sys.stdout
    else:
        outfh = open(outfile, "w")
    
    seq = record.seq
    for feature in record.features:
        if feature.type == "CDS":
            locus_tag = feature.qualifiers["locus_tag"][0]
            product =  feature.qualifiers["product"][0].replace(" ", "_")
            prot_id = feature.qualifiers["protein_id"][0]
            print(f">{record.name}|{feature.location.start+1}-{feature.location.end+1}|{locus_tag}|{product}|{prot_id}", file=outfh)
            print(feature.location.extract(seq), file=outfh)
    
    outfh.close()

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Pull cluster sequences from MIBiG GenBank file")
    parser.add_argument("gb",
                        type=str,
                        help="MIBiG GenBank file (concatenated or individual)")
    parser.add_argument("id",
                        type=str,
                        help="Cluster ID for sequences to be pulled for")
    parser.add_argument("-o", "--output",
                        type=str,
                        default=None,
                        help="Prefix for output fasta file [stdout]")
    return parser.parse_args()

def main():
    args = parse_args()
    record = get_cluster_record(args.gb, args.id)
    if record is None:
        print("error: record {} not found in MIBiG GenBank file".format(args.id))
        sys.exit(1)
    
    print_seqs(record, args.output)

if __name__ == "__main__":
    main()