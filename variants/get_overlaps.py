#!/usr/bin/env python3
"""
Get overlaps of variants for an upset plot.
"""

import argparse
import sys

def get_variants(filename, variants=None):
    """
    Load and keep track of variants in a dictionary. Variants are named as
    chrom_pos_ref_alt.
    """
    if variants is None:
        variants = {}  # {variant -> {filename -> 1 (if variant exists in file)}}
    
    with open(filename, "r") as fh:
        for line in fh:
            if line.startswith("#"):  # Skip headers
                continue
            
            line = line.strip().split("\t")
            chrom = line[0]
            pos = line[1]
            ref = line[3]
            alt = line[4]
            cur_var = f"{chrom}_{pos}_{ref}_{alt}"

            if cur_var not in variants:
                variants[cur_var] = {}
            
            variants[cur_var][filename] = 1
    
    return variants


def print_variant_sets(variants, filenames, outfile):
    """Print a tab-separated table of variant sets."""
    outfile = sys.stdout if outfile is None else open(outfile, "w+")

    print("variant", *filenames, sep="\t", file=outfile, flush=True)
    for variant in variants:
        exists = [variants[variant][f] if f in variants[variant] else 0 for f in filenames]
        print(variant, *exists, sep="\t", file=outfile, flush=True)
    
    if not outfile == sys.stdout:
        outfile.close()


def parse_args():
    parser = argparse.ArgumentParser(description="Get overlapping sets of variants")
    parser.add_argument("files",
                        type=str,
                        nargs="+",
                        help="VCF files")
    parser.add_argument("-o", "--outfile",
                        type=str,
                        default=None,
                        help="Output file [stdout]")
    
    return parser.parse_args()


def main():
    args = parse_args()
    variants = {}
    for f in args.files:
        variants = get_variants(f, variants)
    
    print_variant_sets(variants, args.files, args.outfile)

if __name__ == "__main__":
    main()
