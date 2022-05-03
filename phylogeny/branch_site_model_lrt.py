#!/usr/bin/env python3

import argparse
import glob
import os
import sys

import pandas as pd

from Bio.Phylo.PAML import codeml
from scipy.stats import chi2
from statsmodels.stats.multitest import fdrcorrection

def load_models(null, alt):
    """Load CodeML results for branch-site models from separate files."""
    null = codeml.read(null)
    alt = codeml.read(alt)
    lnL0 = null["NSsites"][2]["lnL"]  # Null hypothesis
    lnL1 = alt["NSsites"][2]["lnL"]  # Alternative hypothesis

    site_class_2a = alt["NSsites"][2]["parameters"]["site classes"][2]
    site_class_2b = alt["NSsites"][2]["parameters"]["site classes"][3]

    return lnL0, lnL1, site_class_2a, site_class_2b


def calculate_lrt(lnL0, lnL1, df):
    """Calculate likelihood ratio test p-value from two log-likelihood values."""
    delta = 2 * (lnL1 - lnL0)
    return chi2.sf(delta, df)


def parse_args():
    parser = argparse.ArgumentParser(description="Calculate LRTs for CodeML branch-site model")
    parser.add_argument("null_dir", help="Directory for null test (M0)")
    parser.add_argument("alt_dir", help="Directory for alternative test (M2a)")
    parser.add_argument("-o", "--outfile", default=None, help="Output file name [stdout]")

    return parser.parse_args()


def main():
    args = parse_args()

    if args.outfile is None:
        args.outfile = sys.stdout

    runs = []
    likelihood_0 = []
    likelihood_1 = []
    p_vals = []
    site_classes_2a = []
    site_classes_2b = []

    # Load results from all runs
    for run in glob.glob(f"{args.null_dir}/*.cds"):
        basename = os.path.basename(run)
        null = f"{run}/mlc"
        alt = f"{args.alt_dir}/{basename}/mlc"
        try:
            lnL0, lnL1, site_class_2a, site_class_2b = load_models(null, alt)
            p_value = calculate_lrt(lnL0, lnL1, 1)
            runs.append(run)
            likelihood_0.append(lnL0)
            likelihood_1.append(lnL1)
            p_vals.append(p_value)
            site_classes_2a.append(site_class_2a["proportion"])
            site_classes_2b.append(site_class_2b["proportion"])
        except:
            print(f"Error parsing {run}", flush=True, file=sys.stderr)
    
    # Correct p-values
    _, p_adj = fdrcorrection(p_vals)

    df = pd.DataFrame(data={
        "run": runs,
        "lnL0": likelihood_0,
        "lnL1": likelihood_1,
        "p_value": p_vals,
        "fdr": p_adj,
        "site_class_2a": site_classes_2a,
        "site_class_2b": site_classes_2b
    })

    df.to_csv(args.outfile, sep="\t", index=False)
    

if __name__ == "__main__":
    main()
