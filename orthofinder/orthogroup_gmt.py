#!/usr/bin/env python3

import sys
import re
import argparse
import requests
import urllib.parse
from orthofinder import load_orthogroups

def load_gff(gff):
    """Load GO terms for genes in a gff file."""
    gene_go_terms = {}
    id_re = "ID=([^;]+)"
    go_re = "(GO:\d+)"
    with open(gff, "r") as fh:
        for line in fh:
            if not line.startswith("#"):
                line = line.strip().split("\t")
                if line[2] == "mRNA":
                    info = line[8]
                    mrna_id = re.match(id_re, info)[1]
                    go_terms = re.findall(go_re, info)
                    if len(go_terms) > 0:
                        for term in go_terms:
                            gene_go_terms[mrna_id] = go_terms
    
    return gene_go_terms

def get_names_to_files(map_file):
    """Read files from map."""
    names_to_files = {}
    with open(map_file, "r") as fh:
        for line in fh:
            line = line.strip().split("\t")
            if len(line) != 2:
                print("error: map file must have one file "
                      "per name in tsv format", sys.stderr)
                sys.exit(1)
            names_to_files[line[0]] = line[1]
    
    return names_to_files

def get_go_label(go_id):
    """Get associated label with GO ID."""
    base_url = "http://api.geneontology.org/api/ontology/term"
    encoded = urllib.parse.quote_plus(go_id)
    res = requests.get(f"{base_url}/{encoded}")
    data = res.json()
    return data["label"] if "label" in data else "NA"

def print_gmt(species, orthogroups, all_mrnas):
    """Print info in gmt format."""
    go_to_genes = {}
    for orth in orthogroups:
        for i, genes in enumerate(orthogroups[orth]):
            cur_sp = species[i]
            if genes[0] == "":
                continue
            else:
                for gene in genes:
                    if gene in all_mrnas[cur_sp]:
                        go_terms = all_mrnas[cur_sp][gene]
                        for go in go_terms:
                            if go not in go_to_genes:
                                go_to_genes[go] = []
                        
                            go_to_genes[go].append(gene)
    
    for go in go_to_genes:
        go_label = get_go_label(go)
        print(go, go_label, *go_to_genes[go], sep="\t")

def parse_args():
    """Get the command line arguments."""
    parser = argparse.ArgumentParser(description="Convert orthogroups to GMT")
    parser.add_argument("tsv",
                        type=str,
                        help="OrthoFinder N0.tsv file containing all orthogroups")
    parser.add_argument("map",
                        type=str,
                        help="File mapping OrthoFinder tsv columns to "
                              "corresponding gff files")
    
    return parser.parse_args()

def main():
    args = parse_args()
    names_to_files = get_names_to_files(args.map)
    orthogroups, species = load_orthogroups(args.tsv)
    orthogroups = {k: [i.split(", ") for i in v] for k, v in orthogroups.items()}

    for s in species:
        if s not in names_to_files:
            print("error: {} not found in map file".format(s),
                file=sys.stderr)
            sys.exit(1)
    
    all_mrnas = {}
    for name in names_to_files:
        mrnas = load_gff(names_to_files[name])
        all_mrnas[name] = mrnas
    
    print_gmt(species, orthogroups, all_mrnas)


if __name__ == "__main__":
    main()