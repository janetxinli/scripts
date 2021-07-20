#!/usr/bin/env python3
"""Get the occurrence of Pfam domains in gProfileR output."""

import argparse
from collections import Counter
import sys
import pandas as pd
from orthofinder import get_pfam_desc

def get_sorted_occurrences(pfam):
    """Count the occurrences of each pfam domain."""
    all_pfam = []
    
    with open(pfam, "r") as fh:
        if not fh.isatty():  # Check for input from stdin
            for line in fh:
                all_pfam.append(line.strip())
        else:
            print("error: input must be piped from stdin or passed as a position argument",
                  file=sys.stderr)
            sys.exit(1)
    
    pfam_counts = Counter(all_pfam)
    df = pd.DataFrame.from_dict(pfam_counts, orient="index", columns=["count"]).reset_index().rename({"index": "pfam_id"}, axis=1)
    df["pfam_desc"] = df["pfam_id"].map(get_pfam_desc)
    df.sort_values("count", ascending=False, inplace=True)

    return df


def parse_args():
    """Parse the command line arguments."""
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("pfam",
                        type=str,
                        nargs="?",
                        default="/dev/stdin",
                        help="File containing line-separated pfam domains [stdin]")
    parser.add_argument("-o", "--outfile",
                        type=str,
                        default="pfam_counts.tsv")
    
    return parser.parse_args()

def main():
    args = parse_args()
    
    pfam_counts = get_sorted_occurrences(args.pfam)
    pfam_counts.to_csv(args.outfile, sep="\t", index=False)

if __name__ == "__main__":
    main()
