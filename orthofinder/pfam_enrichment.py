#!/usr/bin/env python3
"""Perform overrepresentation analysis on Pfam domains."""

import sys
import argparse
import pandas as pd
from collections import Counter
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests
from orthofinder import get_pfam_desc
from gff import functional_info

def get_pfam(filename):
    """Get all pfam domains and associated labels."""
    df = pd.read_csv(filename, sep="\t")
    pfam = df["pfam_domains"].dropna().map(lambda x: x.split(",")).explode()

    return pfam

def get_enriched(all_genes, selection, name, method, cutoff, print_all):
    """Get enrichment for pfam domains."""
    all_counts = Counter(all_genes)
    sel_counts = Counter(selection)
    df = pd.DataFrame({"all": pd.Series(all_counts), name: pd.Series(sel_counts)}).fillna(0)

    # Hypergeometric test
    M = df["all"].sum()
    N = df[name].sum()
    df["p_value"] = df.apply(lambda x: hypergeom(M, x["all"], N).sf(x[name]-1), axis=1)

    # Multiple test correction
    corr = "fdr_bh" if method == "bh" else method
    df[corr] = multipletests(df["p_value"], method=corr)[1]
    df = df.sort_values(corr)
    df["significant"] = df[corr] <= cutoff

    # Add pfam domain and description columns
    df = df.reset_index().rename(columns={"index": "pfam_domain"})
    
    if not print_all:
        df = df.loc[df["significant"]]
   
    df.insert(1, "description", df["pfam_domain"].map(get_pfam_desc))
    return df


def parse_args():
    """Parse the command line arguments."""
    parser = argparse.ArgumentParser(description="Perform Pfam enrichment analysis")
    parser.add_argument("selection",
                        type=str,
                        help="Functional info for selected genes (output from orthogroup_func_info.py)")
    parser.add_argument("gff",
                        nargs="+",
                        type=str,
                        help="GFF files representing background for enrichment analysis")
    parser.add_argument("-n", "--name",
                        type=str,
                        default="selection",
                        help="Name of selected gene group [selection]")
    parser.add_argument("-m", "--method",
                        type=str,
                        choices=["bh", "bonferroni"],
                        default="bh",
                        help="Multiple test correction method [bh (Benjamini-Hochberg FDR)]")
    parser.add_argument("-c", "--cutoff",
                        type=float,
                        default="0.05",
                        help="Significance cutoff for multiple test corrected p-value [0.05]")
    parser.add_argument("-p", "--print_all",
                        action="store_true",
                        help="Print all Pfam domains tested, not just significantly enriched")
    parser.add_argument("-o", "--outfile",
                        type=str,
                        default=sys.stdout,
                        help="Output file [stdout]")

    return parser.parse_args()

def main():
    args = parse_args()

    all_genes = []
    for gff in args.gff:
        func = functional_info(gff, feature="gene")
        for v in func:
            pfam = func[v][1]
            all_genes.extend(pfam)
    
    sel_genes = get_pfam(args.selection)

    enriched = get_enriched(all_genes, sel_genes, pfam_descs,
                            args.method, args.cutoff, args.print_all)
    
    enriched.to_csv(args.outfile, sep="\t", index=False)
 
if __name__ == "__main__":
    main()
