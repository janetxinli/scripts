#!/usr/bin/env python3
"""Perform overrepresentation analysis on Pfam domains."""

import sys
import argparse
import numpy as np
import pandas as pd
from collections import Counter
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
from orthofinder import get_pfam_desc

def get_enriched(all, selection, name, method, cutoff):
    all_counts = Counter(all)
    sel_counts = Counter(selection)
    df = pd.DataFrame({"all": pd.Series(all_counts), name: pd.Series(sel_counts)}).fillna(0)

    # Hypergeometric test
    M = df["all"].sum()
    N = df[name].sum()
    df["p_value"] = df.apply(lambda x: hypergeom(M, df["all"], N).sf(df[name]-1), axis=1)

    # Multiple test correction
    corr = "fdr_bh" if method == "bh" else method
    df[corr] = multipletests(df["p_value"], method=corr)[1]
    df = df.sort_values(corr)
    df["significant"] = df[corr] <= cutoff

    # Add pfam domain and description columns
    df = df.reset_index().rename(columns={"index": "pfam_domain"})
    df.insert(1, "description", df["pfam_domain"].map(lambda x: get_pfam_desc))

    return df


def parse_args():
    """Parse the command line arguments."""
    parser = argparse.ArgumentParser(description="Perform Pfam enrichment analysis")
    parser.add_argument("all",
                        type=str,
                        help="List of all Pfam domains in genes being selected from")
    parser.add_argument("selection",
                        type=str,
                        help="List of Pfam domains in a group of selected genes")
    parser.add_argument("-n", "--name",
                        default="selection",
                        help="Name of selected gene group [selection]")
    parser.add_argument("-m", "--method",
                        choices=["bh", "bonferroni"],
                        default="fdr",
                        help="Multiple test correction method [bh (Benjamini-Hochberg FDR)]")
    parser.add_argument("-c", "--cutoff",
                        default="0.05",
                        help="Significance cutoff for multiple test corrected p-value [0.05]")
    parser.add_argument("-p", "--print_all",
                        action="store_true",
                        help="Print all Pfam domains tested, not just significantly enriched")
    parser.add_argument("-o", "--outfile",
                        default=sys.stdout,
                        help="Output file [stdout]")

    return parser.parse_args()

def main():
    args = parse_args()
    
    enriched = get_enriched(args.all, args.selection, args.name,
                            args.method, args.cutoff)
    
    if args.print_all:
        enriched.write_csv(args.outfile, sep="\t", index=False)
    else:
        enriched.loc[enriched["significant"]].write_csv(args.outfile, sep="\t", index=False)

if __name__ == "__main__":
    main()
