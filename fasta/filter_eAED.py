#!/usr/bin/env python3
"""
Filter sequences in a MAKER fasta file by eAED score.
Expected fasta header format:
>augustus_masked-scaf926-processed-gene-0.3-mRNA-1 protein AED:1.00 eAED:1.00 QI:0|0|0|0|1|1|5|0|1308
"""

import sys
import argparse

def filter_seqs(fasta, filter):
    """Filter sequences with eAED equal or greater than a given value."""
    with open(fasta, "r") as fh:
        for line in fh:
            if line[0] == ">":
                score = float(line.strip().split(" ")[3].split(":")[1])
                print_seq = score < filter
                if print_seq:
                    print(line.strip())
            else:
                if print_seq:
                    print(line.strip())


def get_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Filter out sequences by eAED")
    parser.add_argument("input",
                        nargs="?",
                        type=str,
                        default="/dev/stdin",
                        help="Input fasta file [stdin]")
    parser.add_argument("-f", "--filter",
                        type=float,
                        default=1.0,
                        help="Filter sequences with eAED equal to or greater than this value [1]")
    return parser.parse_args()

def main():
    args = get_args()
    filter_seqs(args.input, args.filter)

if __name__ == "__main__":
    main()
