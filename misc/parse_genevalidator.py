#!/usr/bin/env python3
"""
Parse a GeneValidator output results.csv file.
Expected input format:
AnalysisNumber,GVScore,Identifier,NumberOfHits,LengthCluster,LengthRank,GeneMerge,Duplication,MissingExtraSequences
1,0,protein|genemark-scaf0-processed-gene-0.15-mRNA-1,2,Not enough evidence,Not enough evidence,Not enough evidence,Not enough evidence,Not enough evidence
2,0,protein|augustus-scaf0-processed-gene-0.23-mRNA-1,0,Not enough evidence,Not enough evidence,Not enough evidence,Not enough evidence,Not enough evidence
3,64,protein|augustus-scaf0-processed-gene-0.21-mRNA-1,18,382 [202 - 243],22%,0.0,1.0,77% conserved; 39% extra; 7% missing.
"""

import sys
import argparse
from collections import defaultdict

def parse_gv(csv_file):
    """
    Parse results from results.csv file. Returns a dictionary of prediction
    problems -> count
    """
    summary = defaultdict(int)
    with open(csv_file, "r") as fh:
        fh.readline()  # Read header file
        for i, line in enumerate(fh):
            line = line.strip().split(",")
            num_hits = int(line[3])

            # Number of BLAST hits
            if num_hits < 5:
                summary["too_few_hits"] += 1
                continue
            
            length_info = line[5]
            gene_merge = float(line[6])
            duplication = float(line[7])
            
            # Length rank
            if "too long" in length_info:
                summary["too_long"] += 1
            elif "too short" in length_info:
                summary["too_short"] += 1
            
            # Gene merge (if 0.4 < score < 1.2, query may contain multiple genes)
            if gene_merge > 0.4 and gene_merge < 1.2:
                summary["gene_merge"] += 1
            
            # Duplication (p-value <= 0.05 suggests duplication)
            if duplication <= 0.05:
                summary["duplication"] += 1
    
    summary["total"] = i + 1
    return summary

def print_tsv(summary, outfile, name="GeneValidator"):
    """Prints summary in tsv format to outfile."""
    if outfile is None:
        outfh = sys.stdout
    else:
        outfh = open(outfile, "w")
    
    cols = ["name", "too_few_hits", "too_long", "too_short", "gene_merge", "duplication", "total"]
    print(*cols, sep="\t", file=outfh)
    print(name, end="\t", file=outfh)
    print(*[summary[c] for c in cols[1:]], sep="\t", file=outfh)
    
    outfh.close()


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Parse a GeneValidator output csv file")
    parser.add_argument("csv",
                        type=str,
                        help="GeneValidator output results.csv file")
    parser.add_argument("-n", "--name",
                        type=str,
                        default=None,
                        help="Name of GV run [same as csv file]")
    parser.add_argument("-o", "--outfile",
                        type=str,
                        default=None,
                        help="Output file for results [stdout]")
    return parser.parse_args()

def main():
    args = parse_args()
    if args.name is None:
        args.name = args.csv
    
    summary = parse_gv(args.csv)
    print_tsv(summary, args.outfile, args.name)

if __name__ == "__main__":
    main()