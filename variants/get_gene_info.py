#!/usr/bin/env python3
"""
Get gene names and AED/eAED score for snpEff gene output.
"""

import argparse
import os
import re
import sys
from collections import namedtuple

Gene = namedtuple("Gene", ["name", "aed", "eaed"])

def load_genes(fasta):
    """
    Reads gene names and AED/eAED score from a MAKER fasta 
    file. Returns a dictionary of gene_id -> Gene (name, aed, eaed)
    """
    genes = dict()
    
    with open(fasta, "r") as fh:
        for line in fh:
            if line.startswith(">"):  # fasta header
                line = line.strip()
                gene_id = re.match(">([^\s]+)", line)[1]
                name = re.search("Name:\"(.+)\"", line)[1]
                aed = re.search("AED:(\d*\.\d*)", line)[1]
                eaed = re.search("eAED:(\d*\.\d*)", line)[1]
                genes[gene_id] = Gene(name, aed, eaed)
    
    return genes


def add_gene_info(vargenes, gene_info, outfile):
    if not outfile:
        outfh = sys.stdout
    else:
        outfh = open(outfile, "w+")
    
    with open(vargenes, "r") as fh:
        outfh.write(fh.readline())  # keep first line as-is
        cols = fh.readline().strip().split("\t")
        cols.extend(["maker_gene_name", "AED", "eAED\n"])
        outfh.write("\t".join(cols))
        
        for line in fh:
            line = line.strip().split("\t")
            gene_id = line[2]
            try:
                gene_name, aed, eaed = gene_info[gene_id]
            except KeyError:
                print(f"Error: Gene {gene_id} not found in reference fasta file. Check your inputs.")
                if not outfile is None:
                    os.remove(outfile)
                sys.exit(1)
            line.extend([gene_name, aed, f"{eaed}\n"])
            outfh.write("\t".join(line))


def parse_args():
    parser = argparse.ArgumentParser(description="Get gene names and annotation scores for snpEff gene output")
    
    parser.add_argument("vargenes",
                        type=str,
                        help="snpEff gene summary output txt file")
    parser.add_argument("fasta",
                        type=str,
                        help="MAKER Fasta file for variant calling reference "
                             "containing gene names and scores in headers")
    parser.add_argument("-o", "--outfile",
                        type=str,
                        default=None,
                        help="Output file name [stdout]")
    
    return parser.parse_args()


def main():
    args = parse_args()

    if args.outfile is None:
        if os.path.exists(args.outfile):
            print(f"Error: file {args.outfile} already exists. Exiting.")
    
    gene_info = load_genes(args.fasta)
    add_gene_info(args.vargenes, gene_info, args.outfile)


if __name__ == "__main__":
    main()
