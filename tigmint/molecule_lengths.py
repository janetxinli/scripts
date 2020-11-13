#!/usr/bin/env python3

import sys
import os
import argparse
from reads.read_lengths import print_histogram

def get_molecules(bedfile):
    """Read molecule lengths from a bed file."""
    mol_len = []
    with open(bedfile) as bed:
        for line in bed:
            line_content = line.split("\t")
            mol_len.append(abs(int(line_content[2]) - int(line_content[1])))
    return mol_len


def get_args():
    """Parse the command line arguments."""
    parser = argparse.ArgumentParser(description="Get information from a tigmint-molecule bed file")
    parser.add_argument("data",
                        type=str,
                        choices=["histogram", "lengths"],
                        help="Molecule information (histogram or lengths)")
    parser.add_argument("-b", "--bed",
                        type=str,
                        default=None,
                        help="Bedfile containing molecule extents [stdin]")
    parser.add_argument("-w", "--bin_width",
                        type=int,
                        default=1000,
                        help="Desired bin width for histogram [1000]")
    parser.add_argument("-o", "--outfile",
                        type=str,
                        default=sys.stdout,
                        help="Output file [stdout]")
    return parser.parse_args()


def main():
    args = get_args()
    if args.bed == None:
        args.bed = "/dev/stdin"
    molecule_lengths = get_molecules(args.bed)
    if args.data == "lengths":
        print("\n".join(molecule_lengths), file=args.outfile)
    elif args.data == "histogram":
        print_histogram(molecule_lengths, args.bin_width, args.outfile)


if __name__ == "__main__":
    main()
