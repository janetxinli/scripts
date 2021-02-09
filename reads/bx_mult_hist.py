#!/usr/bin/env python3
"""Calculate barcode multiplicity for an interleaved linked reads fastq file."""

import argparse
import sys
import numpy as np
from utils import print_histogram, read_fasta

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
            for _, _, bx, _ in read_fasta(reads):
                if bx != None:
                    bxs.setdefault(bx, 0)
                    bxs[bx] += 1
        else:
            raise RuntimeError("Reads must be piped from stdin if file name is not provided")
    return bxs


def get_args():
    """Parse the command line arguments."""
    parser = argparse.ArgumentParser(description="Print a histogram of barcode multiplicity for an interleaved fastq file")
    parser.add_argument("-r", "--reads",
                        type=str,
                        default="-",
                        help="Read file to be parsed [stdin]")
    parser.add_argument("-o", "--output_file",
                        type=str,
                        default=None,
                        help="Output file for histogram to be printed to [stdout]")
    parser.add_argument("-w", "--bin_width",
                        type=int,
                        default=1000,
                        help="Bin width for histogram")
    return parser.parse_args()


def main():
    args = get_args()
    bx_mult = bx_multiplicity(args.reads)
    multiplicities = [i for i in bx_mult.values()]
    print_histogram(multiplicities, args.bin_width, args.output_file)

if __name__ == "__main__":
    main()
