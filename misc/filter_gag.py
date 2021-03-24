#!/usr/bin/env python3
"""
Given a set of valid IDs, filter GAG output to keep only features with these IDs.
"""

import os
import re
import sys
import argparse

def get_ids(infile):
    """Read in valid IDs and return a set containing these values."""
    valid_ids = set()
    valid_genes = set()
    gene_re = "(.+)-mRNA-[\d+]"
    with open(infile, "r") as fh:
        if not fh.isatty():  # Check for input from stdin
            for line in fh:
                line = line.strip()
                valid_ids.add(line)
                valid_genes.add(re.search(gene_re, line)[1])
        else:
            print("error: input IDs must be piped from stdin or passed as a positional argument",
                  file=sys.stderr)
            sys.exit(1)

    return valid_ids, valid_genes


def filter_fasta(first_header, infh, valid_ids, outfh):
    """Filter a fasta file."""
    cur_id = first_header.strip()[1:].split(" ")[0]
    print_seq = cur_id in valid_ids
    for i, line in enumerate(infh):
        if i % 2 == 0:  # Seq file
            if print_seq:
                print(">{}".format(cur_id), line.strip(), sep="\n", file=outfh)
        else:
            cur_id = line.strip()[1:].split(" ")[0]
            print_seq = cur_id in valid_ids


def filter_gff(infh, valid_ids, valid_genes, outfh):
    """Filter a gff file."""
    id_re = "ID=([\w\.-]+)[;|:].+"  # fix this to match end of line
    for line in infh:
        line = line.strip().split("\t")
        feature_type = line[2]
        info = line[8]
        if feature_type == "gene":
            cur_gene = info.split("=")[1]
            if cur_gene in valid_genes:
                print(*line, sep="\t", file=outfh)
        else:
            cur_id = re.search(id_re, info)[1]
            if cur_id in valid_ids:
                print(*line, sep="\t", file=outfh)


def run_filter(infile, valid_ids, valid_genes, outfile):
    """Given a .gff or .fasta file and a set of valid IDs, filter out invalid sequences/features."""
    if outfile is None:
        outfh = sys.stdout
    else:
        outfh = open(outfile, "w")

    with open(infile, "r") as fh:
        line = fh.readline()
        if line.startswith(">"):  # fasta input
            filter_fasta(line, fh, valid_ids, outfh)
        elif line.startswith("#"):  # gff input
            print(line.strip(), file=outfh)
            filter_gff(fh, valid_ids, valid_genes, outfh)
        else:
            print("error: input file content doesn't match expected format",
                  file=sys.stderr)
            sys.exit(1)

    outfh.close()


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser("Filter GAG output files given a set of valid IDs to keep")
    parser.add_argument("ids",
                        nargs="?",
                        type=str,
                        default="/dev/stdin", 
                        help="Line-separated IDs to keep")
    parser.add_argument("-i", "--input",
                        type=str,
                        required=True,
                        help="GAG output file for processing (.fasta or .gff)")
    parser.add_argument("-o", "--outfile",
                        type=str,
                        default=None,
                        help="Output file [stdout]")
    return parser.parse_args()


def main():
    args = parse_args()
    if not (args.input.endswith(".fasta") or args.input.endswith(".gff")):
        print("error: input GAG files for filtering must be .gff or .fasta",
              file=sys.stderr)
        sys.exit(1)
    valid_ids, valid_genes = get_ids(args.ids)
    run_filter(args.input, valid_ids, valid_genes, args.outfile)


if __name__ == "__main__":
    main()
