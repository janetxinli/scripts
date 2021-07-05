#!/usr/bin/env python3
"""
Create a lower left similarity matrix for creating a neighbour-joining tree 
in MEGAX from an all-vs-all Mash dist table output file.
"""

import argparse
import sys

def check_file(tsv):
    """Check that tsv file matches expected format."""
    with open(tsv, "r") as fh:
        cols = fh.readline().strip().split("\t")
        if not cols[0] == "#query":
            print(f"mash_mega_matrix.py: error: Mash output {tsv} not in expected format."
                    "Please provide table output")
            sys.exit(1)
        n_samples = len(cols) - 1
        n_rows = 0
        for _ in fh:
            n_rows += 1
        if not n_rows == n_samples:
            print("mash_mega_matrix.py: error: Number of columns and rows don't match")
            sys.exit(1)

def load_distances(tsv):
    """Load distances."""
    distances = dict()
    with open(tsv, "r") as fh:
        samples = fh.readline().strip().split("\t")[1:]
        for i, row in enumerate(fh):
            row = row.strip().split("\t")
            cur_sample = row[0]
            cur_distances = row[1:]
            if cur_sample not in distances:
                distances[cur_sample] = dict()
            j = 0
            while j < i:
                distances[cur_sample][samples[j]] = float(cur_distances[j]) * 100
                j += 1
            
    return distances

def print_mega(distances):
    """Print lower left distance matrix for MEGA."""
    # Print header
    samples = distances.keys()
    print("#mega", "!Title;", "!Description;", sep="\n")
    for s in samples:
        print(f"#{s}")
    print("")
    
    # Print distance matrix
    all_rows = []
    for s1 in samples:
        row = []
        for s2 in samples:
            if s1 != s2:
                if s2 in distances[s1]:
                    row.append(distances[s1][s2])
                else:
                    row.append("")
            else:
                row.append("")
        
        all_rows.append(row)
    
    for row in all_rows:
        print(*row, sep=" ")


def parse_args():
    """Parse the command line arguments."""
    parser = argparse.ArgumentParser(description="Print a MEGAX lower left similarity matrix")
    parser.add_argument("tsv",
                        type=str,
                        help="Mash dist table tsv file")
    parser.add_argument("-o", "--outfile",
                        type=str,
                        default=sys.stdout,
                        help="Ouptut file name [stdout]")
    return parser.parse_args()


def main():
    args = parse_args()
    check_file(args.tsv)
    distances = load_distances(args.tsv)
    print_mega(distances)


if __name__ == "__main__":
    main()