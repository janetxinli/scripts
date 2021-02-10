#!/usr/bin/env python3
# Print a histogram of gene lengths
# Expected format:
# scaf926 .       contig  1       6497    .       .       .       ID=scaf926;Name=scaf926

import sys
import numpy as np

def get_gene_lengths(gff):
    """Returns a list of gene lengths from gff file."""
    gene_lengths = []
    with open(gff, "r") as fh:
        for line in fh:
            if line[0] == ">":
                break
            if line[0] != "#":
                #print(line)
                line = line.split("\t")
                if line[2] == "gene":
                    start = int(line[3])
                    end = int(line[4])
                    gene_lengths.append(end-start)
    return gene_lengths


def print_histogram(gene_lengths, name):
    """Prints histogram of gene lengths with bin size=100."""
    max_length = max(gene_lengths)
    hist_max = (100 - (max_length % 100)) + max_length
    counts, bins = np.histogram(gene_lengths, bins=hist_max//100, range=(0, hist_max))
    for i, b in enumerate(bins[1:]):
        print(int(b), counts[i], name, sep="\t", flush=True)


def main():
    if len(sys.argv) != 3:
        print("Usage: {} <gff file> <name>".format(sys.argv[0]), file=sys.stderr)
        sys.exit(1)
    gff = sys.argv[1]
    name = sys.argv[2]
    gene_lengths = get_gene_lengths(gff)
    print_histogram(gene_lengths, name)


if __name__ == "__main__":
    main()
