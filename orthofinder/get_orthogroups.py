#!/usr/bin/env python3

import sys
import argparse
from orthofinder import *

def load_orthogroups(tsv):
    """Load all orthogroups from file."""
    orthogroups = {}  # OG -> [names of genes in each species]

    with open(tsv, "r") as fh:
        header = fh.readline().strip().split("\t")
        total_species = len(header) - 3  # 3 extra columns in N0.tsv file

        for line in fh:
            line = line.strip().split("\t")
            ortho_genes = line[3:]
            og = line[1]
            orthogroups[og] = ortho_genes
    
    return orthogroups


def parse_args():
    """Get the command line arguments."""
    parser = argparse.ArgumentParser(description="Print core or accessory orthogroups")
    parser.add_argument("tsv",
                        type=str,
                        help="OrthoFinder PhylogeneticHierarchicalOrthogroups/N0.tsv result file")
    parser.add_argument("-t", "--orth_type",
                        type=str,
                        default="core_single",
                        choices=["core_single", "core_all", "accessory", "singleton"],
                        help="Orthogroup type of interest [core_single]")
    parser.add_argument("-o", "--outfile",
                        type=str,
                        default=None,
                        help="Output file [stdout]")
    
    return parser.parse_args()

def main():
    args = parse_args()
    tests = {
        "core_single": is_single_core,
        "core_all": is_core,
        "accessory": is_accessory,
        "singleton": is_singleton
    }
    
    og_to_print = {k for k, v in get_orthogroup_genes(args.tsv).items() if tests[args.orth_type](v)}
    orthogroups = load_orthogroups(args.tsv)

    outfile = sys.stdout if args.outfile is None else open(args.outfile, "w+")

    for o in orthogroups:
        if o in og_to_print:
            print(o, *orthogroups[o], sep="\t", file=outfile)
    
    if not args.outfile is None:
        outfile.close()


if __name__ == "__main__":
    main()