#!/usr/bin/env python3
"""
Finds complete protein sequences (starts with M, ends with *) or complete gene sequences
(starts with ATG, ends with TAG, TAA, TGA).
"""

import sys
import argparse
from utils import read_fasta

def find_complete_sequences(fasta, outfile, mode):
    """Prints complete protein sequences (starts with M, ends with *)."""
    if mode == "p":  # Protein start and stop codons
        start_codon = "M"
        end_codons = {"*"}
        clen = 1
    else:  # mRNA start and stop codons
        start_codon = "ATG"
        end_codons = {"TAG", "TAA", "TGA"}
        clen = 3

    if outfile == "-":
        outfh = sys.stdout
    else:
        outfh = open(outfile, "w")
    
    with open(fasta, "r") as infile:
        for header, seq, _, _ in read_fasta(infile):
            cur_start = seq[:clen].upper()
            cur_end = seq[-clen:].upper()
            if cur_start == start_codon and cur_end in end_codons:
                print(">" + header, file=outfh)
                print(seq, file=outfh)
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
    parser.add_argument("-m", "--mode",
                        type=str,
                        choices=["p", "m"],
                        default="p",
                        help="Mode for finding complete sequences ('p' for proteins, "
                             "'m' for mrna sequences) [p]")
    return parser.parse_args()

def main():
    args = parse_args()
    find_complete_sequences(args.fasta, args.outfile, args.mode)

if __name__ == "__main__":
    main()

