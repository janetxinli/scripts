#!/usr/bin/env python3
"""Score genes according to gProfiler GO enrichment p-values."""

import sys
import numpy as np
import pandas as pd

def split_gene_terms(df):
    """Split enrichment terms into a list."""
    df = df.loc[df["go_terms"].notnull()]
    df["go_terms"] = df["go_terms"].map(lambda x: x.split(","))

    return df

def get_sig(df):
    """Return a dictionary of GO term -> gProfiler enrichment p-value (fdr)."""
    return dict(zip(df["term_id"], -np.log10(df["p_value"])))

def calculate_scores(df, sig):
    """Calculate enrichment scores for genes."""
    df["score"] = df["go_terms"].map(lambda x: sum([sig[i] if i in sig else 0 for i in x]))
    df.sort_values("score", ascending=False, inplace=True)
    df["go_terms"] = df["go_terms"].map(lambda x: ",".join(x))

    return df.loc[df["score"] > 0]

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <functional info> <enrichment results>")
        sys.exit(1)
    
    func = pd.read_csv(sys.argv[1], sep="\t")
    enrichment = pd.read_csv(sys.argv[2], sep="\t")

    func = split_gene_terms(func)
    sig = get_sig(enrichment)

    scores = calculate_scores(func, sig)
    scores.to_csv(sys.stdout, sep="\t", index=False)
    
