#!/usr/bin/env python3

import sys
import argparse
from orthofinder import is_core, is_single_core, is_accessory, is_singleton, get_orthogroup_genes

def load_orthogroups(tsv):
    """Load all orthogroups from file."""
    orthogroups = {}  # OG -> [names of genes in each species]

    with open(tsv, "r") as fh:
        header = fh.readline().strip().split("\t")
        species = header[1:]

        for line in fh:
            line = line.strip().split("\t")
            ortho_genes = line[1:]
            og = line[0]
            orthogroups[og] = ortho_genes
    
    return orthogroups, species

def print_table(orthogroups, og_to_print, species, outfile=None):
    """Print orthogroups in table format."""

    outfh = sys.stdout if outfile is None else open(outfile, "w+")
    print("OG", *species, sep="\t", file=outfh)
    
    for o in orthogroups:
        if o in og_to_print:
            print(o, *orthogroups[o], sep="\t", file=outfh)
    
    if not outfile is None:
        outfh.close()

def print_list(orthogroups, og_to_print, outfile=None):
    """Print orthogroup genes in list format."""
    outfh = sys.stdout if outfile is None else open(outfile, "w+")
    orthogroups = [[i.split(", ") for i in v] for k, v in orthogroups.items() if k in og_to_print]
    for og in orthogroups:
        for genes in og:
            if genes[0] != "":
                print(*genes, sep="\n", file=outfh)
    
    if not outfile is None:
        outfh.close()

def is_var_core(orthogroup):
    return is_core(orthogroup) and not is_single_core(orthogroup)

def parse_args():
    """Get the command line arguments."""
    parser = argparse.ArgumentParser(description="Print core or accessory orthogroups")
    parser.add_argument("tsv",
                        type=str,
                        help="OrthoFinder Orthogroups.tsv result file")
    parser.add_argument("-t", "--orth_type",
                        type=str,
                        default="core_single",
                        choices=["core_single", "core_all", "core_var", "accessory", "singleton"],
                        help="Orthogroup type of interest [core_single]")
    parser.add_argument("-o", "--outfile",
                        type=str,
                        default=None,
                        help="Output file [stdout]")
    parser.add_argument("-f", "--format",
                        default="table",
                        choices=["table", "t", "list", "l"],
                        help="Output format: table (t) of orthogroups by species "
                             "or a simple list (l) of genes")
    
    return parser.parse_args()

def main():
    args = parse_args()
    tests = {
        "core_single": is_single_core,
        "core_all": is_core,
        "core_var": is_var_core,
        "accessory": is_accessory,
        "singleton": is_singleton
    }
    
    og_to_print = {k for k, v in get_orthogroup_genes(args.tsv).items() if tests[args.orth_type](v)}

    orthogroups, species = load_orthogroups(args.tsv)

    if args.format == "table" or args.format == "t":
        print_table(orthogroups, og_to_print, species, args.outfile)
    
    else:
        print_list(orthogroups, og_to_print, args.outfile)


if __name__ == "__main__":
    main()