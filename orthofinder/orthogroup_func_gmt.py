#!/usr/bin/env python3
"""
Creates a GMT file for HOGs from orthogroup_func_info.py output.
"""

import argparse
import sys
from orthofinder import get_go_label

def load_hog_func(tsv):
    """
    Loads functional info for HOGs. Returns a tuple of hogs, go_terms.
    hogs is a dict of HOG -> [associated GO terms]; go_terms is a set of all
    unique go terms present.
    """
    hogs = {}
    go_terms = set()

    with open(tsv, "r") as fh:
        fh.readline()  # Skip header
        for line in fh:
            line = line.strip().split("\t")
            if len(line) < 3 or line[2] == "":
                continue
            cur_go_terms = line[2].split(",")
            hogs[line[0]] = cur_go_terms
            go_terms.update(cur_go_terms)
    
    return hogs, go_terms

def search_go_terms(go_terms):
    """Searches for GO descriptions."""
    id_to_desc = {}
    
    for go in go_terms:
        desc = get_go_label(go)
        if desc == go:
            desc == "NA"
        
        id_to_desc[go] = desc
    
    return id_to_desc

def print_gmt(hogs, go_terms, outfile):
    go_to_hogs = {}  # go_term -> [HOGs]

    outfh = open(outfile, "w+") if outfile is not None else sys.stdout

    for hog in hogs:
        for go in hogs[hog]:
            if go not in go_to_hogs:
                go_to_hogs[go] = []
            
            go_to_hogs[go].append(hog)
    
    for go in go_to_hogs:
        print(go, go_terms[go], *go_to_hogs[go], sep="\t", flush=True, file=outfh)

    if outfh != sys.stdout:
        outfh.close()

def parse_args():
    parser = argparse.ArgumentParser(description="Convert orthogroups to GMT")
    parser.add_argument("tsv",
                        type=str,
                        help="orthogroup_func_info.py tsv file containing functional info for ogs")
    parser.add_argument("-o", "--outfile",
                        type=str,
                        default=None,
                        help="Output file [stdout]")
    
    return parser.parse_args()

def main():
    args = parse_args()
    hogs, go_terms = load_hog_func(args.tsv)
    go_terms = search_go_terms(go_terms)
    print_gmt(hogs, go_terms, args.outfile)

if __name__ == "__main__":
    main()
