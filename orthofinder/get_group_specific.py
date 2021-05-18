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
    tree = ElementTree.fromstring(res.content)
    
    return tree[0][0].text.strip()


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
    group.drop(group.columns[tuple([drop_cols])], axis=1, inplace=True)  # Tuple here?
    group_func = group.merge(func, on=["HOG", "OG"], how="left")

    return group_func

def get_func_info(df, colname, info_func, dropcol):
    """Get functional information and merge into orthogroup df."""
    df_copy = df.copy()
    df_copy[colname] = df_copy[colname].dropna().map(lambda x: x.split(","))
    df_copy = df_copy.explode(colname)

    # Explode comma-separated column
    terms = df_copy[colname].dropna().unique()
    func_terms = pd.DataFrame(terms, columns=[colname])
    func_terms["label"] = func_terms[colname].map(info_func)
    df_copy = df_copy.merge(func_terms, on=colname).drop(dropcol, axis=1)

    return df_copy

def get_mrna_note(df, species, gff):
    mrna_to_note = mrna_note(gff)

    # Split and unnest gene column
    df[species] = df[species].map(lambda x: x.split(",")).explode(species)
    idx = df.columns.get_loc(species)
    new_col = species + "_gene_note"
    df.insert(idx + 1, new_col, df[species].map(lambda x: mrna_to_note[x] if x in mrna_to_note else ""))
    
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

    return parser.parse_args()

def main():
    args = parse_args()
    header = get_header(args.tsv)

    if len(args.species) != len(args.gff):
        print("get_group_specific.py: error: mismatching number of species and "
              "gff files provided")
        sys.exit(1)

    # Get indices of species of interest
    indices = [i for i, sp in enumerate(header) if sp in args.species]
    
    if len(indices) != len(args.species):
        print(f"get_group_specific.py: error: not all species found "
              "in {args.tsv}")
        sys.exit(1)
    
    hog = pd.read_csv(args.tsv, sep="\t")
    func = pd.read_csv(args.func, sep="\t")
    group_func = clean_dfs(hog, func, indices)

    print("getting GO term labels")
    go = get_func_info(group_func, "go_terms", get_go_label, "pfam_domains")

    # print("getting Pfam accession descriptions")
    # pfam = get_func_info(group_func, "pfam_domains", get_pfam_desc, "go_terms")

    for i, sp in enumerate(args.species):
        go = get_mrna_note(go, sp, args.gff[i])
        # pfam = get_mrna_note(pfam, sp, args.gff[i])

    go.to_csv(f"{args.group}_orthogroups.go_terms.tsv", sep="\t", index=False)
    # pfam.to_csv(f"{args.group}_orthogroups.pfam_domains.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()
