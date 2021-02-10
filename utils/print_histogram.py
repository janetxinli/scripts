#!/usr/bin/env python3
"""Print a histogram with a specified bin width given a list of values."""

import sys
import argparse
import numpy as np

def print_histogram(values, binwidth, outfile=None):
    """
    Given a list of values and desired bin width, prints a histogram of values
    in tsv format. By default prints histogram to stdout.
    """
    max_val = np.max(values)
    # Calculate bin width and max
    upper_lim = max_val + (binwidth - (max_val % binwidth))
    bins = upper_lim // binwidth
    counts, bins = np.histogram(values, int(bins), range=(0, upper_lim))
    if outfile is None:
        of = sys.stdout
    else:
        of = open(outfile, "w+")
    for i, bin in enumerate(bins[1:]):
        print(int(bin), counts[i], sep="\t", file=of)
    of.close()


def parse_args():
    """Parse the command line arguments."""
    parser = argparse.ArgumentParser(description="Print a histogram given a list of values")
    parser.add_argument("-i", "--input",
                        default="/dev/stdin",
                        type=str,
                        help="List of values (line-separated) to print a histogram of [stdin]")
    parser.add_argument("-o", "--outfile",
                        default=None,
                        type=str,
                        help="File to print histogram to [stdout]")
    parser.add_argument("-b", "--binwidth",
                        type=int,
                        required=True,
                        help="Bin width for histogram")
    return parser.parse_args()


def main():
    args = parse_args()
    hist_vals = []
    with open(args.input, "r") as infile:
        for line in infile:
            try:
                hist_vals.append(float(line.strip()))
            except ValueError:
                print("error: values must be numeric")
                sys.exit(1)
    
    print_histogram(hist_vals, args.binwidth, args.outfile)


if __name__ == "__main__":
    main()
