#!/usr/bin/env python3
"""Calculate k-mer coverage and GC content for a FASTA/FASTQ file."""

import argparse
import sys
from utils import calc_gc, read_fasta

def get_kmers(sequence, k):
    """Returns a generator for iterating over k-mers in a sequence."""
    if type(sequence) != str:
        raise ValueError("Input sequence must be a string")
    
    if type(k) != int:
        raise ValueError("k must be a length")
    
    i = 0
    while i <= (len(sequence) - k):
        i += 1
        yield sequence[i:i+k]


def count_kmer_coverage(fh, k):
    """
    Counts the occurrences of kmers in a given FASTA or FASTQ file handle.
    Returns a dictionary of kmer -> coverage.
    """
    kmer_cov = dict()
    for _, seq, _, _ in read_fasta(fh):
        for kmer in get_kmers(seq, k):
            if kmer not in kmer_cov:
                kmer_cov[kmer] = 1
            else:
                kmer_cov[kmer] += 1
    
    return kmer_cov


def print_output(kmer_cov, outfile):
    """
    Given a dictionary of k-mers -> coverage, calculate GC content of k-mers
    and print in tab-separated format to outfile.
    """
    if outfile != sys.stdout:
        outfile = open(outfile)
    
    for kmer in kmer_cov:
        gc = calc_gc(kmer)
        print(kmer_cov[kmer], gc, sep="\t", file=outfile, flush=True)
    
    if outfile != sys.stdout:
        outfile.close()


def parse_args():
    """Parse the command line arguments."""
    parser = argparse.ArgumentParser(description="Calculate the coverage and GC content of "
                                     "k-mers in FASTQ/FASTA sequences")
    parser.add_argument("input",
                        type=str,
                        help="Input FASTA/FASTQ file (use '-' to read from stdin)")
    parser.add_argument("-k", "--kmer_length",
                        type=int,
                        default=80,
                        help="K-mer length [80]")
    parser.add_argument("-o", "--outfile",
                        type=str,
                        default=sys.stdout,
                        help="Output file [stdout]")
    
    return parser.parse_args()


def main():
    args = parse_args()
    if args.input == "-":
        args.input = sys.stdin
    else:
        args.input = open(args.input)

    kmer_cov = count_kmer_coverage(args.input, args.kmer_length)
    print_output(kmer_cov, args.outfile)

    args.input.close()


if __name__ == "__main__":
    main()
