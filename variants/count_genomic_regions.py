#!/usr/bin/env python3
"""
Count the number of variants in exons, introns and intergenic regions.
"""

import argparse
from collections import namedtuple
from intervaltree import IntervalTree

Variant = namedtuple("Variant", ["contig", "pos", "indel"])

def add_intervals(exons, intervals):
        """
        Given a gene's start and end position, list of exons' start and end positions
        and an IntervalTree, adds introns and exons to the tree. Exons have a data
        attribute value of "exon", and introns have a value of "intron".
        """
        # Sort exons by start position
        exons.sort(key=lambda x: x[0])

        i = 0
        while i < (len(exons) - 1):
            current_exon = exons[i]
            next_exon = exons[i+1]

            # Add current exon
            if current_exon[0] != current_exon[1]:
                intervals.addi(current_exon[0], current_exon[1], "exon")
            # Add next intron
            intervals.addi(current_exon[1]+1, next_exon[0]-1, "intron")
            i += 1

        # Add final exon
        current_exon = exons[i]    
        intervals.addi(current_exon[0], current_exon[1], "exon")


def load_genic_regions(gff):
    """
    Load genic regions (introns and exons) into an IntervalTree. Returns a
    dictionary of contig_id -> gene_intervals
    """
    genic_regions = {}
    current_contig = None
    current_intervals = IntervalTree()
    current_gene = None
    current_exons = []
    
    with open(gff, "r") as fh:
        for line in fh:
            if not line.startswith("#"):
                line = line.strip().split("\t")
                
                if line[0] != current_contig:  # At a new contig
                    if not current_contig is None:
                        # Add final gene to current_intervals
                        if len(current_exons):
                            add_intervals(current_exons, current_intervals)

                        # Keep track of prev contig gene regions
                        genic_regions[current_contig] = current_intervals
                    
                    current_contig = line[0]
                    current_intervals = IntervalTree()
                    current_gene = None
                    current_exons = []
                
                if line[2] == "gene":  # At a new gene
                    if not current_gene is None:
                        add_intervals(current_exons, current_intervals)

                    current_gene = (int(line[3]), int(line[4]))  # start, end
                    current_exons = []
                
                if line[2] == "exon":  # At a new exon of the current gene
                    current_exons.append((int(line[3]), int(line[4])))
        
        # Add final contig (plus final gene) to current_intervals
        if len(current_exons):
            add_intervals(current_exons, current_intervals)
        genic_regions[current_contig] = current_intervals
    
    return genic_regions

def load_variants(vcf):
    """Load variants from vcf. Returns a list of Variants."""
    variants = []
    with open(vcf, "r") as fh:
        for line in fh:
            if not line.startswith("#"):
                line = line.strip().split("\t")
                contig = line[0]
                pos = int(line[1])
                indel = line[7].startswith("INDEL")
                variants.append(Variant(contig, pos, indel))
    
    return variants

def count_variant_locations(variants, genic_regions):
    """
    Count the number of occurrences of SNPs and indels in genic and
    intergenic regions. Returns a dict of variant_type -> position -> count.
    """
    variant_counts = {
        "snp": {
            "exon": 0,
            "intron": 0,
            "intergenic": 0
            },
        "indel": {
            "exon": 0,
            "intron": 0,
            "intergenic": 0
            }
        }
    
    for v in variants:
        key = "indel" if v.indel else "snp"
        query = genic_regions[v.contig].at(v.pos)
        
        if len(query) == 0:
            position = "intergenic"
        else:
            position= query.pop().data
        
        variant_counts[key][position] += 1

    return variant_counts


def parse_args():
    parser = argparse.ArgumentParser(description="Count the number of variants in exons, introns and intergenic regions")
    parser.add_argument("vcf",
                        type=str,
                        help="VCF file containing variants")
    parser.add_argument("gff",
                        type=str,
                        help="GFF file for reference genome")
    
    return parser.parse_args()

def main():
    args = parse_args()

    # Load genomic regions in from GFF
    genic_regions = load_genic_regions(args.gff)

    # Load variants
    variants = load_variants(args.vcf)

    variant_location_counts = count_variant_locations(variants, genic_regions)
    print("location", "snp", "indel", sep="\t")
    locations = ["exon", "intron", "intergenic"]

    for loc in locations:
        print(loc, variant_location_counts["snp"][loc], variant_location_counts["indel"][loc], sep="\t")

if __name__ == "__main__":
    main()
