#!/usr/bin/env python3
"""Calculate barcode multiplicity for an interleaved linked reads fastq file."""

import argparse
import sys
import numpy as np
from fasta import read_fasta

def bx_multiplicity(rfile):
    """
    Calculate barcode multiplicity given a read file name. Returns a dictionary
    of barcode -> num reads with that barcode.
    """
    bxs = {}
    if rfile == "-":
        rfile = "/dev/stdin"
    with open(rfile) as reads:
        if not reads.isatty():
            for head, _, bx, _ in read_fasta(reads):
                if bx != None:
                    bxs.setdefault(bx, 0)
                    bxs[bx] += 1
        else:
            raise RuntimeError("Reads must be piped from stdin if file name is not provided")
    return bxs


def print_hist(bx_mult, binwidth, outfile):
    """
    Given a dictionary of barcode multiplicity where barcode -> num reads,
    print a histogram as a tsv file to stdout.
    """
    multiplicity = [i for i in bx_mult.values()]
    max_mult = np.max(multiplicity)
    # Calculate bin width and max
    upper_lim = max_mult + (binwidth - (max_mult % binwidth))
    bins = upper_lim // binwidth
    counts, bins = np.histogram(multiplicity, bins, range=(0, upper_lim))
    if outfile == "-":
        of = sys.stdout
    else:
        of = open(outfile, "w+")
    for i, bin in enumerate(bins[1:]):
        print(int(bin), counts[i], sep="\t", file=of)
    of.close()


def get_args():
    """Parse the command line arguments."""
    parser = argparse.ArgumentParser(description="Print a histogram of barcode multiplicity for an interleaved fastq file")
    parser.add_argument("-r", "--reads",
                        type=str,
                        default="-",
                        help="Read file to be parsed [stdin]")
    parser.add_argument("-o", "--output_file",
                        type=str,
                        default="-",
                        help="Output file for histogram to be printed to [stdout]")
    parser.add_argument("-w", "--bin_width",
                        type=int,
                        default=1000,
                        help="Bin width for histogram")
    return parser.parse_args()


def main():
    args = get_args()
    print_hist(bx_multiplicity(args.reads), args.bin_width, args.output_file)

if __name__ == "__main__":
    main()
