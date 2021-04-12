#!/usr/bin/env python3
"""
Pull MIBiG cluster nucleotide sequences from MIBiG GenBank file
@author: janetxinli
"""

import sys
import argparse
from Bio import SeqIO
from antismash import get_cluster_record

def print_seqs(record):
    """Extract and print feature sequences from GenBank record."""
    seq = record.seq
    i = 0
    for feature in record.features:
        if feature.type == "CDS":
            locus_tag = feature.qualifiers["locus_tag"][0] if "locus_tag" in feature.qualifiers else "NA"
            product =  feature.qualifiers["product"][0].replace(" ", "_") if "product" in feature.qualifiers else "NA"
            prot_id = feature.qualifiers["protein_id"][0] if "protein_id" in feature.qualifiers else "NA"
            outfile = record.name + ".seq" + str(i) + ".fa"
            outfh = open(outfile, "w")
            print(f">{record.name}|{feature.location.start+1}-{feature.location.end+1}|{locus_tag}|{product}|{prot_id}", file=outfh)
            print(feature.location.extract(seq), file=outfh)
            outfh.close()
            i += 1

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Pull cluster sequences from MIBiG GenBank file")
    parser.add_argument("gb",
                        type=str,
                        help="MIBiG GenBank file (concatenated or individual)")
    parser.add_argument("id",
                        type=str,
                        help="Cluster ID for sequences to be pulled for")
    return parser.parse_args()

def main():
    args = parse_args()
    record = get_cluster_record(args.gb, args.id)
    if record is None:
        print("error: record {} not found in MIBiG GenBank file".format(args.id))
        sys.exit(1)
    
    print_seqs(record)

if __name__ == "__main__":
    main()