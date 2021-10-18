#!/usr/bin/env python3
# Get nucleotide sequences for orthogroups

import argparse
import sys
from utils import read_fasta

def load_fasta_sequences(fasta):
    """Returns a dictionary of header -> sequence."""
    seqs = {}
    with open(fasta, "r") as fh:
        for header, seq, _, _ in read_fasta(fh):
            seqs[header] = seq
    
    return seqs


def print_orthogroup_genes(tsv, species_seqs, fastas):
    """Prints nucleotide sequences of each orthogroup to an individual fasta."""
    with open(tsv, "r") as fh:
        fh.readline()  # Skip header
        for line in fh:
            line = line.strip().split("\t")
            hog = line[0]
            outfile = open(f"{hog}.cds.fa", "w")
            for i, gene in enumerate(line[2:]):
                try:
                    seq = species_seqs[i][gene]
                except KeyError:
                    print(f"get_orthogroup_cds.py: error: gene {gene} not found in {fastas[i]}")
                    sys.exit(1)
                print(f">{gene}", file=outfile, flush=True)
                print(seq, file=outfile, flush=True)


def parse_args():
    parser = argparse.ArgumentParser(description="Print orthogroup nucleotide sequences. "
                                     "Will print each OG to a single fasta file with the suffix "
                                     ".cds.fa")
    parser.add_argument("tsv",
                        type=str,
                        help="Orthogroup table (output of get_orthogroups.py)")
    parser.add_argument("fastas",
                        nargs="+",
                        help="mRNA sequences for species in orthogroup table. "
                             "Must be provided in the same order")
    return parser.parse_args()


def main():
    args = parse_args()

    # Check that number of species in tsv match number of fasta files
    with open(args.tsv, "r") as fh:
        header = fh.readline().strip().split("\t")
        if (len(header) - 2) != len(args.fastas):
            print("get_orthogroup_cds.py: error: number of fasta files does not "
                  "match number of species in orthogroup table")
            sys.exit(1)
    
    seqs = []
    for f in args.fastas:
        seqs.append(load_fasta_sequences(f))

    print_orthogroup_genes(args.tsv, seqs, args.fastas)


if __name__ == "__main__":
    main()
