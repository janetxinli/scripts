#!/usr/bin/env python3
"""
Extract gene information from a MIBiG cluster Genbank file and print
in .tsv format. Expects the .gbk file to contain a single record.
"""

import argparse
from Bio import SeqIO

def parse_cluster(filename, field):
    cluster = {}  # identifier -> product
    record = SeqIO.read(filename, "genbank")
    for feature in record.features:
        if feature.type == "CDS":
            identifier = feature.qualifiers[field][0]
            cluster[identifier] = feature.qualifiers["product"][0]
    
    return cluster


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("gbk",
                        type=str,
                        help="Cluster Genbank file")
    parser.add_argument("-f", "--field",
                        type=str,
                        choices=["locus_tag", "protein_id"],
                        default="locus_tag",
                        help="Identifier field")
    return parser.parse_args()

def main():
    args = parse_args()
    cluster = parse_cluster(args.gbk, args.field)
    print(args.field, "product", sep="\t")
    for k in cluster:
        print(k, cluster[k], sep="\t")

if __name__ == "__main__":
    main()

