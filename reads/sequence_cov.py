#!/usr/bin/env python3
"""
Author: Janet Li (@janetxinli)
Compute sequence coverage for a genome sequencing read file.
"""

import sys
import argparse

def calculate_cov(fh, genome_size, outfile):
    """Calculate the sequence coverage of a given file."""
    line = fh.readline()
    sum_bases = 0

    if line[0] == "@":  # Fastq
        for i, line in enumerate(fh):
            if i % 4 == 0:
                sum_bases += len(line)
    elif line[0] == ">":  #Fasta
        for i, line in enumerate(fh):
            if i % 2 == 0:
                sum_bases += len(line)
    
    outfile.write(str(sum_bases / genome_size) + "\n")


def get_gsize(size_string):
    """Get genome size from a string"""
    try:
        return float(size_string)
    except ValueError:
        print("error: genome size must be an integer or in scientific notation (e.g. 3e9)", file=sys.stderr, flush=True)


def get_args():
    """Parse the command-line arguments."""
    parser = argparse.ArgumentParser(description="Calculate the sequence coverage for a set of reads")
    parser.add_argument("-r", "--reads",
                        type=argparse.FileType("r"),
                        default=sys.stdin,
                        help="Read file [stdin]")
    parser.add_argument("-g", "--gsize",
                        type=str,
                        required=True,
                        help="Genome size (in base pairs)")
    parser.add_argument("-o", "--outfile",
                        type=argparse.FileType("w"),
                        default=sys.stdout,
                        help="File to print sequence coverage to [stdout]")
    return parser.parse_args()


def main():
    args = get_args()
    calculate_cov(args.reads, get_gsize(args.gsize), args.outfile)

if __name__ == "__main__":
    main()
