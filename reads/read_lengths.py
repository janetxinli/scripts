#!/usr/bin/env python3
"""Print a tsv file of read lengths or a histogram of read lengths."""

import sys
import argparse
import numpy as np

def get_readlengths(read_file, outfile=None, print_lengths=True):
    """Calculate and either print read lengths or return lengths as a list."""
    lengths = []
    with open(read_file, "rt") as long_reads:
        if not long_reads.isatty():
            if print_lengths:
                print("length", file=outfile, flush=True)
            first = long_reads.readline()
            if first[0] == "@":  # Fastq file
                for i, line in enumerate(long_reads):
                    if i % 4 == 0:
                        cur_length = len(line.strip())
                        if print_lengths:
                            print(cur_length, file=outfile, flush=True)
                        else:
                            lengths.append(cur_length)
            elif first[0] == ">":  # Fasta file
                for i, line in enumerate(long_reads):
                    if i % 2 == 0:
                        cur_length = len(line.strip())
                        if print_lengths:
                            print(cur_length, file=outfile, flush=True)
                        else:
                            lengths.append(cur_length)
            if not print_lengths:
                return lengths
        else:
            print("read_lengths: error: reads must be piped from stdin or given with the -r argument",
                file=sys.stderr, flush=True)
            sys.exit(1)


def print_histogram(read_lengths, bin_width, outfile):
    """Print a histogram of read lengths and width of bin_width to outfile."""
    max_length = max(read_lengths)
    hist_max = (bin_width - (max_length % bin_width)) + max_length
    counts, bins = np.histogram(read_lengths, bins=hist_max//bin_width, range=(0, hist_max))
    print("length", "count", sep="\t", file=outfile)
    for i, b in enumerate(bins[1:]):
        print(int(b), counts[i], sep="\t", file=outfile)


def get_args():
    """Get command line arguments."""
    parser = argparse.ArgumentParser(description="Print read lengths or a read length histogram \
        for a long read file")
    parser.add_argument("-r", "--reads",
                        type=str,
                        default="-",
                        help="Read file [stdin]")
    parser.add_argument("-g", "--hist",
                        action="store_true",
                        default=False,
                        help="Print a histogram of read lengths")
    parser.add_argument("-b", "--bin_width",
                        type=int,
                        default=1000,
                        help="Desired bin width for histogram [1000]")
    parser.add_argument("-o", "--outfile",
                        type=str,
                        default=sys.stdout,
                        help="Output file [stdout]")
    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    if args.reads == "-":
        args.reads = "/dev/stdin"
    if args.hist:
        lengths = get_readlengths(args.reads, print_lengths=False)
        print_histogram(lengths, args.bin_width, args.outfile)
    else:
        get_readlengths(args.reads, outfile=args.outfile)
