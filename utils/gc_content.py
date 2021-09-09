#!/usr/bin/env python3
"""
Calculate gc content from a FASTA, FASTQ or BAM file.
"""

import argparse
import os.path
import sys
import Bio.SeqIO
import pysam

def calc_gc(sequence):
    """Given a DNA sequence, returns the GC content."""
    if type(sequence) != str:
        raise ValueError("Input sequence must be a string")
    
    length = len(sequence)
    gc = 0

    for char in sequence.upper():
        if char == "C" or char == "G":
            gc += 1
    
    return gc / length


def get_fastx_seqs(filename, format):
    """Parses a FASTA/FASTQ file and returns a list of tuples (name, sequence)."""
    if not (format == "fasta" or format == "fastq"):
        raise ValueError("Format must be 'fasta' or 'fastq")
    
    seqs = []
    with open(filename, "r") as fh:
        for entry in Bio.SeqIO.parse(fh, format):
            seqs.append((entry.id, str(entry.seq)))
    
    return seqs


def get_aln_seqs(filename):
    """Parses a SAM/BAM file and return a list of tuples (name, sequence)."""
    seqs = []
    aln_file = pysam.AlignmentFile(filename, check_sq=False)

    for read in aln_file.fetch(until_eof=True):
        seqs.append((read.query_name, read.query_sequence))
    
    return seqs


def parse_args():
    parser = argparse.ArgumentParser(description="Calculate GC content from FASTA, FASTQ or BAM/SAM files.")
    parser.add_argument("input",
                        type=str,
                        help="Input file (use '-' for stdin)")
    parser.add_argument("filetype",
                        type=str,
                        choices=["fasta", "fastq", "sam", "bam"],
                        help="Input file type")
    parser.add_argument("-p", "--per_seq",
                        action="store_true",
                        help="Calculate GC content per sequence and print in tabular (.tsv) format")
    
    return parser.parse_args()


def main():
    args = parse_args()
    
    # Get input file
    if args.input == "-":
        args.input = sys.stdin
    
    else:
        if not os.path.exists(args.input):
            print(f"gc_content.py: error: file {args.input} does not exist")
            sys.exit(1)

    # Get sequences
    if args.filetype == "fasta" or args.filetype == "fastq":
        seqs = get_fastx_seqs(args.input, args.filetype)
    else:
        seqs = get_aln_seqs(args.input)
    
    if not args.per_seq:
        entire_seq = "".join([seq for (name, seq) in seqs])
        gc = calc_gc(entire_seq)
        print(gc)
    else:
        for name, seq in seqs:
            gc = calc_gc(seq)
            print(name, gc, sep="\t")


if __name__ == "__main__":
    main()
