#!/usr/bin/env python3

import sys
import re
import argparse
from collections import namedtuple

Interval = namedtuple("Interval", ["info", "start", "end"])

def load_repeats(gff):
    """Load repeats into a dictionary."""
    repeats = {}  # scaf -> [Interval]
    with open(gff, "r") as fh:
        for line in fh:

            if line.startswith("#"):
                continue
            line = line.strip().split("\t")
            cur_scaf = line[0]
            if cur_scaf not in repeats:
                repeats[cur_scaf] = []
            start = int(line[3])
            end = int(line[4])
            info = line[8]
            repeats[cur_scaf].append(Interval(info, start, end))


    # Sort intervals
    for scaf in repeats:
        repeats[scaf].sort(key=lambda x: x.start)
    
    return repeats

def load_gene_names(gene_names="-"):
    """Returns a set of gene names of interest."""
    genes_to_search = set()

    if gene_names == "-":
        if sys.stdin.isatty():
            print("distance_to_repeats.py: error: gene names must be piped from "
                "stdin or passed as a positional argument")
            sys.exit(1)
        fh = sys.stdin
    else:
        fh = open(gene_names, "r")

    for line in fh:
        genes_to_search.add(line.strip())

    fh.close()
    return genes_to_search


def load_orthogroup_genes(gene_names, gff):
    """Load information of mRNAs in gene_names from annotation gff."""
    gene_info = {}  # scaf -> [Interval]
    id_re = "ID=(JB-\w+-R[A-Z]+)\;.+"
    with open(gff, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            line = line.strip().split("\t")
            if line[2] == "mRNA":
                gene_id = re.match(id_re, line[8])[1]
                if gene_id in gene_names:
                    scaf = line[0]
                    start = int(line[3])
                    end = int(line[4])
                    if scaf not in gene_info:
                        gene_info[scaf] = []
                    
                    gene_info[scaf].append(Interval(gene_id, start, end))

        # Sort genes by start pos
        for scaf in gene_info:
            gene_info[scaf].sort(key=lambda x: x.start)

    return gene_info

def print_distances(genes, repeats):
    """Finds distances between genes and repeats."""
    repeat_re = "\"(.+)\""
    print("scaffold\tgene\tclosest_dist\tclosest_rep\toverlap")
    for scaf in genes:
        if scaf in repeats:
            i = 0
            repeat_intervals = repeats[scaf]
            for gene in genes[scaf]:
                num_reps = len(repeat_intervals)
                gene_start = gene.start
                gene_end = gene.end
                found = False
                while i < num_reps and not found:
                    if repeat_intervals[i].start >= gene_end:
                        found = True
                    else:
                        i += 1
                
                left = repeat_intervals[i-1] if i > 0 else None
                left_dist = gene_start - left.end if left is not None else "NA"
                left_info = re.search(repeat_re, left.info)[1] if left is not None else "NA"
                right = repeat_intervals[i] if i < num_reps else None
                right_dist = right.start - gene_end if right is not None else "NA"
                right_info = re.search(repeat_re, right.info)[1] if right is not None else "NA"
                
                if left_dist == "NA":
                    closest_dist = right_dist
                    closest_info = right_info
                    overlap = right_dist < 0
                elif right_dist == "NA":
                    closest_dist = left_dist
                    closest_info = left_info
                    overlap = left_dist < 0
                else:
                    if abs(left_dist) < abs(right_dist):
                        closest_dist = left_dist
                        closest_info = left_info
                        overlap = left_dist < 0
                    else:
                        closest_dist = right_dist
                        closest_info = right_info
                        overlap = right_dist < 0
                
                print(scaf, gene.info, abs(closest_dist), closest_info, overlap, sep="\t")


def parse_args():
    """Parse the command line arguments."""
    parser = argparse.ArgumentParser(description="Calculate the average distance "
                                     "of genes to repeat regions")
    parser.add_argument("repeats",
                        type=str,
                        help="RepeatModeler GFF file")
    parser.add_argument("annotation",
                        type=str,
                        help="Annotation GFF file")
    parser.add_argument("genes",
                        nargs="?",
                        type=str,
                        default="-",
                        help="Genes names to search in annotation gff file (one gene per line) [stdin]")
    return parser.parse_args()

def main():
    args = parse_args()
    genes_to_search = load_gene_names(args.genes)
    genes = load_orthogroup_genes(genes_to_search, args.annotation)
    repeats = load_repeats(args.repeats)
    print_distances(genes, repeats)

if __name__ == "__main__":
    main()
