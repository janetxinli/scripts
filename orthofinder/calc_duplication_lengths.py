#!/usr/bin/env python3
"""
Calculate total length of duplicated orthogroups on Terminal nodes (events that
occurred within a species).
"""

import sys
import argparse
import re
from gff import ID_RE

def read_dups(infile="-"):
    """Read duplicated genes names."""
    genes = []
    if infile is None:
        if sys.stdin.isatty():
            print("calc_duplication_lengths.py: error: gene names must "
                  "be piped from stdin or passed as an argument")
            sys.exit(1)
        fh = sys.stdin
    else:
        fh = open(infile, "r")

    for line in fh:
        genes.append(line.strip())
    
    fh.close()
    return genes

def get_gene_lengths(gff):
    """Load all gene lengths from GFF file."""
    gene_lengths = {}  # maker id -> gene length
    with open(gff, "r") as fh:
        for line in fh:
            if not line.startswith("#"):
                line = line.strip().split("\t")
                if line[2] == "gene":
                    gene = re.search(ID_RE, line[8])[1]
                    length = int(line[4]) - int(line[3])
                    gene_lengths[gene] = length
    
    return gene_lengths

def parse_args():
    """Parse the command line arguments."""
    parser = argparse.ArgumentParser(description="Calculation duplication length on "
                                     "a terminal OrthoFinder species tree node")
    parser.add_argument("node",
                        type=str,
                        help="Species tree node")
    parser.add_argument("gff",
                        type=str,
                        help="Annotation GFF file corresponding to species "
                             "tree node")
    parser.add_argument("dups",
                        type=str,
                        nargs="?",
                        default=None,
                        help="Duplicated genes (line-separated) [stdin]")                             
    
    return parser.parse_args()

def main():
    args = parse_args()
    duplicated_genes = read_dups(args.dups)
    lengths = get_gene_lengths(args.gff)
    
    # Calculate total dup length
    total_length = 0
    for gene in duplicated_genes:
        total_length += lengths[gene]
    
    print(total_length)

if __name__ == "__main__":
    main()
