#!/usr/bin/env python3
"""
Map orthogroups to GO terms.
"""

import argparse
import sys
from orthofinder import load_orthogroups, load_gff_go

def print_orthogroup_go(orthogroups, go_terms):
    """Prints orthogroup GO terms."""
    print("HOG", "OG", "GO terms", sep="\t")
    for og in orthogroups:
        og_mrnas = orthogroups[og]
        og_go_terms = {}  # GO term -> occurrences
        for mrna in og_mrnas:
            if mrna != "" and mrna in go_terms:
                for go in go_terms[mrna]:
                    if go not in og_go_terms:
                        og_go_terms[go] = 0
                    
                    og_go_terms[go] += 1
        
        if len(og_go_terms) > 0:
            max_count = max(og_go_terms.values())
            og_go = [k for k, v in og_go_terms.items() if v == max_count]
            go_to_print = ",".join(og_go)
        else:
            go_to_print = ""
        
        print(og[0], og[1], go_to_print, sep="\t")

def parse_args():
    """Parse the command line arguments."""
    parser = argparse.ArgumentParser(description="Map orthogroups to most common GO terms of "
                                     "orthogroup genes")
    parser.add_argument("tsv",
                        type=str,
                        help="OrthoFinder N0.tsv file")
    parser.add_argument("gff",
                        nargs="*",
                        type=str,
                        help="GFF files for OrthoFinder run "
                             "(in the same order as N0.tsv columns)")
    
    return parser.parse_args()

def main():
    args = parse_args()
    orthogroups, species = load_orthogroups(args.tsv)
    if len(species) != len(args.gff):
        print("orthogroup_go_terms.py: error: species in N0.tsv and gff files "
              "provided is unequal")
        sys.exit(1)
    
    # Load go terms in each gff file
    all_go_terms = {}
    for gff in args.gff:
        all_go_terms.update(load_gff_go(gff))
    
    for og in orthogroups:
        genes = orthogroups[og]
    
    print_orthogroup_go(orthogroups, all_go_terms)


if __name__ == "__main__":
    main()