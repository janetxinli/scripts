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
    attribute value of "exon", and introns have a value of "intron". Returns length
    of exons and introns.
    """
    exon_lengths = 0
    intron_lengths = 0

    # Sort exons by start position
    exons.sort(key=lambda x: x[0])

    i = 0
    while i < (len(exons) - 1):
        current_exon = exons[i]
        next_exon = exons[i+1]

        # Add current exon
        if current_exon[0] != current_exon[1]:
            intervals.addi(current_exon[0], current_exon[1], "exon")
            exon_lengths += abs(current_exon[1] - current_exon[0])
        # Add next intron
        intervals.addi(current_exon[1]+1, next_exon[0]-1, "intron")
        intron_lengths += abs((next_exon[0]-1) - (current_exon[1]+1))
        i += 1

    # Add final exon
    current_exon = exons[i]
    intervals.addi(current_exon[0], current_exon[1], "exon")
    exon_lengths += abs(current_exon[1] - current_exon[0])

    return exon_lengths, intron_lengths

def load_genic_regions(gff):
    """
    Load genic regions (introns and exons) into an IntervalTree. Returns a
    tuple of genic_regions (contig_id -> gene_intervals), feature_lengths
    (contig_id -> {type (intron, exon, intergenic) -> length})  TODO change this description
    """
    genic_regions = {}
    feature_lengths = {
        "total": 0,
        "exon": 0,
        "intron": 0
    }

    current_contig = None
    current_intervals = IntervalTree()
    current_gene = None
    current_exons = []

    current_exon_lengths = 0
    current_intron_lengths = 0
    
    with open(gff, "r") as fh:
        for line in fh:
            if not line.startswith("#"):
                line = line.strip().split("\t")
                
                if line[2] == "region":  # Add new contig length
                    feature_lengths["total"] += int(line[4])
                
                if line[0] != current_contig:  # At a new contig
                    if not current_contig is None:
                        # Add final gene to current_intervals
                        if len(current_exons):
                            exon_lengths, intron_lengths = add_intervals(current_exons, current_intervals)
                            current_exon_lengths += exon_lengths
                            current_intron_lengths += intron_lengths

                        # Keep track of prev contig gene regions
                        genic_regions[current_contig] = current_intervals
                        feature_lengths["exon"] += current_exon_lengths
                        feature_lengths["intron"] += current_intron_lengths
                    
                    current_contig = line[0]
                    current_intervals = IntervalTree()
                    current_gene = None
                    current_exons = []

                    current_exon_lengths = 0
                    current_intron_lengths = 0
                
                if line[2] == "gene":  # At a new gene
                    if not current_gene is None:
                        exon_lengths, intron_lengths = add_intervals(current_exons, current_intervals)
                        current_exon_lengths += exon_lengths
                        current_intron_lengths += intron_lengths

                    current_gene = (int(line[3]), int(line[4]))  # start, end
                    current_exons = []
                
                if line[2] == "exon":  # At a new exon of the current gene
                    current_exons.append((int(line[3]), int(line[4])))
        
        # Add final contig (plus final gene) to current_intervals
        if len(current_exons):
            exon_lengths, intron_lengths = add_intervals(current_exons, current_intervals)
            current_exon_lengths += exon_lengths
            current_intron_lengths += intron_lengths

        genic_regions[current_contig] = current_intervals
        feature_lengths["exon"] += current_exon_lengths
        feature_lengths["intron"] += current_intron_lengths
    
    # Tally up total genic/intergenic length
    total_feature_lengths = {
        "genic": {
            "total": 0,
            "exon": 0,
            "intron": 0,
        },
        "intergenic": {
            "total": 0
        }
    }

    genic = feature_lengths["exon"] + feature_lengths["intron"]
    intergenic = feature_lengths["total"] - genic
    total_feature_lengths["intergenic"]["total"] = intergenic
    total_feature_lengths["genic"]["total"] = genic
    total_feature_lengths["genic"]["exon"] = feature_lengths["exon"]
    total_feature_lengths["genic"]["intron"] = feature_lengths["intron"]

    return genic_regions, total_feature_lengths

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
    intergenic regions. Returns a dict of variant_type -> position -> count.  TODO fix this output note
    """
    variant_counts = {
        "snp": {
            "intergenic": {
                "total": 0
            },
            "genic": {
                "exon": 0,
                "intron": 0,
            }
        },
        "indel": {
            "intergenic": {
                "total": 0
            },
            "genic": {
                "exon": 0,
                "intron": 0,
            }
        }
    }
    
    for v in variants:
        key = "indel" if v.indel else "snp"
        query = genic_regions[v.contig].at(v.pos)
        
        if len(query) == 0:
            region = "intergenic"
            type = "total"
        else:
            region = "genic"
            type = query.pop().data
        
        variant_counts[key][region][type] += 1
    
    # Sum up genic variants
    for k in variant_counts:
        variant_counts[k]["genic"]["total"] = variant_counts[k]["genic"]["exon"] + variant_counts[k]["genic"]["intron"]

    return variant_counts

def get_snps_per_kb(snp_count, bp_length):
    return snp_count / (bp_length / 1000)

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
    genic_regions, total_feature_lengths = load_genic_regions(args.gff)

    # Load variants
    variants = load_variants(args.vcf)

    variant_location_counts = count_variant_locations(variants, genic_regions)
    print("region", "type", "snp", "indel", "total_length", "snps_per_kb", sep="\t")
    
    # Define output table structure
    regions = {
        "genic": ["total", "exon", "intron"],
        "intergenic": ["total"]
    }

    for r in regions:
        for type in regions[r]:
            length = total_feature_lengths[r][type]
            snps_per_kb = get_snps_per_kb(variant_location_counts["snp"][r][type], length)
            print(r, type, variant_location_counts["snp"][r][type], 
                  variant_location_counts["indel"][r][type], length, snps_per_kb, sep="\t")


if __name__ == "__main__":
    main()
