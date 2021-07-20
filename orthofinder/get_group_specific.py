#!/usr/bin/env python3
"""
Get orthogroups specific to a group of species. Group-specific OGs are 
present in 1+ copy each species of the group, and not present in species
not in the group.

Requires output from orthogroup_func_info.py
"""

import sys
import argparse
from xml.etree import ElementTree
import numpy as np
import pandas as pd
import requests
from gff import mrna_note

def get_header(tsv):
    """Get a list of species in OrthoFinder NO.tsv file."""
    with open(tsv, "r") as fh:
        header = fh.readline().strip().split("\t")
        
        return header


def get_go_label(go_id):
    """Get the Gene Ontology label given a GO id."""
    res = requests.get(f"http://api.geneontology.org/api/ontology/term/{go_id}")
    
    return res.json()["label"]


def get_pfam_desc(acc):
    """Get the Pfam description given an accession term."""
    res = requests.get(f"https://pfam.xfam.org/family/{acc}?output=xml")
    
    try:
        tree = ElementTree.fromstring(res.content)
        desc = tree[0][0].text.strip()
    except ElementTree.ParseError:
        desc = "" 
    
    return desc


def clean_dfs(hog, func, indices):
    """
    Merge orthogroups with functional information, and find
    group-specific orthogroups.
    """
    drop_cols = [i for i in range(3, hog.shape[1]) if i not in indices]
    arr = np.ones(hog.shape[0], dtype=bool)
    
    for i in range(3, hog.shape[1]):
        if i in indices:
            arr = arr & ~hog.iloc[:, i].isnull()
        else:
            arr = arr & hog.iloc[:, i].isnull()
    
    group = hog[arr]
    group.drop(group.columns[tuple([drop_cols])], axis=1, inplace=True)
    group_func = group.merge(func, on=["HOG", "OG"], how="left")

    return group_func

def get_func_info(df, colname, new_colname, info_func,):
    """Get functional information and merge into orthogroup df."""
    df_copy = df.copy()
    unique_terms = df[colname].dropna().map(lambda x: x.split(",")).explode().unique()
    term_names = {}

    # Map accession to term name/description
    for acc in unique_terms:
        term_names[acc] = info_func(acc)
    
    idx = df_copy.columns.get_loc(colname)
    df_copy \
        .insert(idx + 1, new_colname, df_copy[colname] \
        .map(lambda x: ",".join([term_names[i] for i in x.split(",")]), na_action="ignore"))

    return df_copy


def get_mrna_note(df, species, gff):
    mrna_to_note = mrna_note(gff)

    idx = df.columns.get_loc(species)
    new_col = species + "_gene_notes"
    df.insert(idx + 1, new_col, df[species].map(lambda x: ",".join([mrna_to_note[i] for i in x.split(", ") if i in mrna_to_note])))
    
    return df

def parse_args():
    parser = argparse.ArgumentParser(description="Get group-specific orthogroups")
    parser.add_argument("tsv",
                        type=str,
                        help="OrthoFinder N0.tsv file")
    parser.add_argument("func",
                        type=str,
                        help="Orthogroup functional info file")
    parser.add_argument("group",
                        type=str,
                        help="Group prefix for output files")
    parser.add_argument("-s", "--species",
                        nargs="*",
                        type=str,
                        help="Species belonging to group of interest")
    parser.add_argument("-g", "--gff",
                        nargs="*",
                        type=str,
                        help="GFF files for species (must be given in the same "
                             "order as the --species argument")
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        help="Print verbose output")


    return parser.parse_args()

def main():
    args = parse_args()

    if (args.verbose):
        print("getting header from N0.tsv", file=sys.stdout)
    
    header = get_header(args.tsv)

    if len(args.species) != len(args.gff):
        print("get_group_specific.py: error: mismatching number of species and "
              "gff files provided")
        sys.exit(1)

    # Get indices of species in group
    indices = [i for i, sp in enumerate(header) if sp in args.species]

    if len(indices) != len(args.species):
        print(f"get_group_specific.py: error: not all species found "
              "in {args.tsv}")
        sys.exit(1)
    
    if (args.verbose):
        print("reading hierarchical orthogroups into memory", file=sys.stdout)
    
    hog = pd.read_csv(args.tsv, sep="\t")

    if (args.verbose):
        print("reading functional information into memory", file=sys.stdout)
    
    func = pd.read_csv(args.func, sep="\t")
    
    if (args.verbose):
        print("finding group-specific orthogroups", file=sys.stdout)
    
    group_func = clean_dfs(hog, func, indices)

    # Print all core orthogroups to file
    group_func.drop(["go_terms", "pfam_domains"], axis=1) \
              .to_csv(f"{args.group}_core_orthogroups.tsv", sep="\t", index=False)

    if (args.verbose):
        print("getting GO term labels", file=sys.stdout)
    
    og_func = get_func_info(group_func, "go_terms", "go_labels", get_go_label)

    if (args.verbose):
        print("getting Pfam accession descriptions", file=sys.stdout)
    
    og_func = get_func_info(og_func, "pfam_domains", "pfam_descs", get_pfam_desc)

    if (args.verbose):
        print("adding gene information for core orthogroups", file=sys.stdout)
    
    for i, sp in enumerate(args.species):
        if (args.verbose):
            print(f"current species: {sp}", file=sys.stdout)
        
        og_func = get_mrna_note(og_func, sp, args.gff[i])

    og_func.to_csv(f"{args.group}_core_orthogroups.func.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()
