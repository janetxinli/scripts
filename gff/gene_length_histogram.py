#!/usr/bin/env python3
# Print a histogram of gene lengths to stdout
# Expected format:
# scaf926 .       contig  1       6497    .       .       .       ID=scaf926;Name=scaf926

import sys
import argparse
import numpy as np
from utils import print_histogram

def get_gene_lengths(gff):
    """Returns a list of gene lengths from gff file."""
    gene_lengths = []
    with open(gff, "r") as fh:
        for line in fh:
            if line[0] == ">":  # Skip fasta sequences at end of file
                break
            if line[0] != "#":
                #print(line)
                line = line.split("\t")
                if line[2] == "gene":
                    start = int(line[3])
                    end = int(line[4])
                    gene_lengths.append(end-start)
    return gene_lengths


def print_named_histogram(gene_lengths, binwidth, name):
    """Prints histogram of gene lengths with a given name."""
    max_length = max(gene_lengths)
    hist_max = (100 - (max_length % 100)) + max_length
    counts, bins = np.histogram(gene_lengths, bins=hist_max//binwidth, range=(0, hist_max))
    for i, b in enumerate(bins[1:]):
        print(int(b), counts[i], name, sep="\t", flush=True)


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Extract and print gene lengths as a histogram")
    parser.add_argument("-g", "--gff",
                        type=str,
                        default="/dev/stdin",
                        help="Input GFF file")
    parser.add_argument("-n", "--name",
                        type=str,
                        default=None,
                        help="Name of gff file to print to histogram")
    parser.add_argument("-b", "--binwidth",
                        type=int,
                        default=100,
                        help="Bin width for histogram")
    return parser.parse_args()


def main():
    args = parse_args()
    gene_lengths = get_gene_lengths(args.gff)
    if args.name is None:
        print_histogram(gene_lengths, args.binwidth)
    else:
        print_named_histogram(gene_lengths, args.binwidth, args.name)


if __name__ == "__main__":
    main()
