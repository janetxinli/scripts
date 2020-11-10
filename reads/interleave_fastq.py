#!/usr/bin/env python3
"""Interleave paired end fastq files."""

import sys
import argparse
import gzip
import glob

def interleave_reads(read1, read2, outfile):
    """Combined reads from r1 and r2 in interleaved format and print to outfile."""
    if outfile != sys.stdout:
        outfile = open(outfile, "w+")
    with gzip.open(read1, "rt") as r1, gzip.open(read2, "rt") as r2:
        for i, line in enumerate(r1):
                outfile.write(line)
                if i % 4 == 3:
                    for _ in range(4):
                        outfile.write(r2.readline())


def get_args():
    """Parse arguments from command line."""
    parser = argparse.ArgumentParser(description="Interleave two paired end fastq read files")
    parser.add_argument("r1",
                        type=str,
                        help="Read 1 file (gzip)")
    parser.add_argument("r2",
                        type=str,
                        help="Read 2 file (gzip)")
    parser.add_argument("-o", "--outfile",
                        type=str,
                        default=sys.stdout,
                        help="Output file for interleaved reads to be printed [stdout]")
    return parser.parse_args()


def main():
    args = get_args()
    args.r1 = glob.glob(args.r1)[0]
    args.r2 = glob.glob(args.r2)[0]
    interleave_reads(args.r1, args.r2, args.outfile)


if __name__ == "__main__":
    main()
