#!/usr/bin/env python3

import sys
import re
import argparse
import requests
import urllib.parse
from orthofinder import load_orthogroups
from gff import ID_RE, GO_RE

def load_gff_go(gff):
    """Load GO terms for genes in a gff file."""
    gene_go_terms = {}  # gene ID -> [associated go terms]

    with open(gff, "r") as fh:
        for line in fh:
            if not line.startswith("#"):
                line = line.strip().split("\t")
                if line[2] == "gene":
                    info = line[8]
                    gene_id = re.match(ID_RE, info)[1]
                    go_terms = re.findall(GO_RE, info)
                    if len(go_terms) > 0:
                        for _ in go_terms:
                            gene_go_terms[gene_id] = go_terms
    
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

def print_gmt(species, orthogroups, all_genes):
    """Print info in gmt format."""
    go_to_genes = {}
    for orth in orthogroups:
        for i, mrnas in enumerate(orthogroups[orth]):
            cur_sp = species[i]
            if mrnas[0] == "":
                continue
            else:
                for m in mrnas:
                    gene = re.search("(.+)-R[A-Z]+", m)[1]  # Convert mRNA ID to gene ID
                    if gene in all_genes[cur_sp]:
                        go_terms = all_genes[cur_sp][gene]
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
    orthogroups = {k: [i.split(", ") for i in v] for k, v in orthogroups.items()}  # Split csv mRNAs into a list

    for s in species:
        if s not in names_to_files:
            print("error: {} not found in map file".format(s),
                file=sys.stderr)
            sys.exit(1)
    
    all_genes = {}
    for name in names_to_files:
        genes = load_gff_go(names_to_files[name])
        all_genes[name] = genes
    
    print_gmt(species, orthogroups, all_genes)


if __name__ == "__main__":
    main()