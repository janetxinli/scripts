#!/usr/bin/env python3
"""
Map orthogroups to GO terms and Pfam domains from a MAKER gff file.
Selects the most common GO and Pfam annotations.
"""

import argparse
import sys
from orthofinder import load_orthogroups
from gff import functional_info

def get_most_common(terms):
    """
    Given a dict mapping terms to occurrences, returns
    the most common terms.
    """
    max_count = max(terms.values())
    return [k for k, v in terms.items() if v == max_count]

def print_orthogroup_func(orthogroups, func_info):
    """Prints orthogroup GO terms."""
    print("OG", "go_terms", "pfam_domains", sep="\t")
    for og in orthogroups:
        og_mrnas = orthogroups[og]
        og_go_terms = {}  # GO term -> occurrences
        og_pfam_doms = {}  # Pfam domain -> occurrences
        for mrna in og_mrnas:
            if mrna != "":
                mrnas = mrna.split(", ")
                for m in mrnas:
                    if m in func_info:
                        for go in func_info[m][0]:
                            if go not in og_go_terms:
                                og_go_terms[go] = 0
                    
                            og_go_terms[go] += 1
                
                        for pfam in func_info[m][1]:
                            if pfam not in og_pfam_doms:
                                og_pfam_doms[pfam] = 0
                    
                            og_pfam_doms[pfam] += 1 

        if len(og_go_terms) > 0:
            max_count = max(og_go_terms.values())
            og_go = [k for k, v in og_go_terms.items() if v == max_count]
            go_to_print = ",".join(og_go)
        else:
            go_to_print = ""
        
        go_to_print = ",".join(og_go_terms) if len(og_go_terms) > 0 else ""
        pfam_to_print = ",".join(og_pfam_doms) if len(og_pfam_doms) > 0 else ""

        print(og, go_to_print, pfam_to_print, sep="\t")

def parse_args():
    """Parse the command line arguments."""
    parser = argparse.ArgumentParser(description="Map orthogroups to most common "
                                     "GO terms of orthogroup genes")
    parser.add_argument("tsv",
                        type=str,
                        help="OrthoFinder Orthogroups.tsv file")
    parser.add_argument("gff",
                        nargs="*",
                        type=str,
                        help="GFF files for OrthoFinder run "
                             "(in the same order as Orthogroups.tsv columns)")
    
    return parser.parse_args()

def main():
    args = parse_args()
    orthogroups, species = load_orthogroups(args.tsv)
    
    if len(species) != len(args.gff):
        print("orthogroup_go_terms.py: error: species in Orthogroups.tsv and gff files "
              "provided is unequal")
        sys.exit(1)
    
    # Load go terms in each gff file
    all_func_info = {}
    for gff in args.gff:
        all_func_info.update(functional_info(gff))
   
    print_orthogroup_func(orthogroups, all_func_info)


if __name__ == "__main__":
    main()

