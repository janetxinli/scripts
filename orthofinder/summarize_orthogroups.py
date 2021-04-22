#!/usr/bin/env python3

import sys
import argparse

def is_core(orthogroup):
    """Tests whether an orthogroup is present in all species."""
    for sp in orthogroup:
        if sp == 0:
            return False
    
    return True

def is_single_core(orthogroup):
    """
    Tests whether an orthogroup is present in a single copy
    in all species.
    """
    for sp in orthogroup:
        if sp != 1:
            return False
    
    return True

def is_accessory(orthogroup):
    """Tests whether 2+ genes have an orthogroup."""
    contains = 0
    for sp in orthogroup:
        if sp > 0:
            contains += 1
    
    return contains >= 2

def is_singleton(orthogroup):
    return orthogroup.count(0) == len(orthogroup) - 1

def read_orthogroups(tsv):
    """
    Read orthogroups file. Returns a dictionary mapping orthogroups
    to number of genes in each species.
    """
    orthogroups = {}  # orthogroup --> [num genes in each orthogroup]
    
    with open(tsv, "r") as fh:
        header = fh.readline().strip().split("\t")
        total_species = len(header) - 3
        
        for line in fh:
            line = line.strip().split("\t")
            og = line[1]
            og_info = []
            species_in_line = len(line) - 3
            for sp in line[3:]:
                if sp.split(", ")[0] == "":
                    num_genes = 0
                else:
                    num_genes = len(sp.split(", "))
                og_info.append(num_genes)
            
            if species_in_line < total_species:
                for s in range(species_in_line, total_species):
                    og_info.append(0)
            
            orthogroups[og] = og_info
    
    return orthogroups


def summarize_orthogroups(orthogroups, outfile=None):
    """Print summary of orthogroups in tsv format."""
    if outfile is None:
        outfh = sys.stdout
    else:
        outfh = open(outfile, "w+")
    
    core = 0
    single_core = 0
    acc = 0
    singleton = 0

    for og in orthogroups:
        if is_core(orthogroups[og]):
            core += 1
            if is_single_core(orthogroups[og]):
                single_core += 1
        elif is_accessory(orthogroups[og]):
            acc += 1
        elif is_singleton(orthogroups[og]):
            singleton += 1
    
    print("type\tcategory\tcount", file=outfh)
    print(f"core\tall\t{core}", file=outfh)
    print(f"core\tsingle_copy\t{single_core}", file=outfh)
    print(f"accessory\tall\t{acc}", file=outfh)
    print(f"singleton\tall\t{singleton}", file=outfh)

    outfh.close()
        

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Summarize OrthoFinder results")
    parser.add_argument("tsv",
                        type=str,
                        help="OrthoFinder PhylogeneticHierarchicalOrthogroups/N0.tsv result file")
    parser.add_argument("-o", "--outfile",
                        type=str,
                        default=None,
                        help="Output results file")
    
    return parser.parse_args()


def main():
    args = parse_args()
    orthogroups = read_orthogroups(args.tsv)
    summarize_orthogroups(orthogroups, args.outfile)


if __name__ == "__main__":
    main()