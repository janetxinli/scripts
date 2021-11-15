#!/usr/bin/env python3

import argparse
import sys
import pandas as pd
from gprofiler import GProfiler

def load_genes(genes):
    gene_list = []
    if genes == "-":
        if not sys.stdin.isatty():
            fh = sys.stdin
        else:
            print("run_gprofiler.py: error: no genes detected from stdin")
            sys.exit(1)
    else:
        fh = open(genes, "r")
    
    for line in fh:
        gene_list.append(line.strip())
    
    if genes != "-":
        fh.close()
    
    return gene_list

def run_gp(genes, organism, sorted, correction_method, cutoff, all):
    """Run gProfiler analysis."""
    gp = GProfiler(return_dataframe=True)
    res = gp.profile(genes, organism, ordered=sorted, user_threshold=cutoff, 
               all_results=all, significance_threshold_method=correction_method).drop("parents", axis=1)
    
    return res

def parse_args():
    parser = argparse.ArgumentParser("Perform enrichment analysis with gProfiler and print results as tsv")
    parser.add_argument("genes",
                        type=str,
                        help="Gene list file, one gene per line (use '-' for stdin)")
    parser.add_argument("organism",
                        type=str,
                        help="gProfiler organism")
    parser.add_argument("-s", "--sorted",
                        action="store_true",
                        help="Runs ordered query; assumes input list is sorted by significance")
    parser.add_argument("-m", "--corr_method",
                        choices=["fdr", "bonferroni", "g_SCS"],
                        default="fdr")
    parser.add_argument("-o", "--output",
                        type=str,
                        default=None,
                        help="Output file name [stdout]")
    parser.add_argument("-c", "--cutoff",
                        type=float,
                        default=0.05,
                        help="Significance cutoff [0.05]")
    parser.add_argument("-a", "--all",
                        action="store_true",
                        help="Print results for all terms, not just significant results")
    
    return parser.parse_args()


def main():
    args = parse_args()
    genes = load_genes(args.genes)
    gp = run_gp(genes, args.organism, args.sorted, args.corr_method, args.cutoff, args.all)

    outfile = sys.stdout if args.output is None else args.output
    gp.to_csv(outfile, sep="\t", index=False)

if __name__ == "__main__":
    main()
