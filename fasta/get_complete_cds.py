#!/usr/bin/env python3
"""
Finds complete protein sequences (starts with M, ends with *) from GAG output genome.proteins.fasta
Expects sequences are in single line fasta format (one line per sequence).
"""

import sys
import argparse

def find_complete_sequences(fasta, outfile):
    """Prints complete protein sequences (starts with M, ends with *)."""
    if outfile == "-":
        outfh = sys.stdout
    else:
        outfh = open(outfile, "w")
    with open(fasta, "r") as infh:
        for i, line in enumerate(infh):
            if i % 2 == 0:
                header = line.strip()
            else:
                seq = line.strip()
                if seq.startswith("M") and seq.endswith("*") and "*" not in seq[1:-1]:
                    print(header, seq, sep="\n", file=outfh)
    outfh.close()

def parse_args():
    """Parse the command line arguments."""
    parser = argparse.ArgumentParser(description="Filter complete protein sequences in a fasta file")
    parser.add_argument("fasta",
                        type=str,
                        help="Single-line fasta file containing protein sequences")
    parser.add_argument("-o", "--outfile",
                        type=str,
                        default="-",
                        help="Output file for filtered protein sequences [stdout]")
    return parser.parse_args()

def main():
    args = parse_args()
    find_complete_sequences(args.fasta, args.outfile)

if __name__ == "__main__":
    main()

