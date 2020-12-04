#!/usr/bin/env python3
"""
Count variant types (SNP, indel), zygosity and effects in SnpSift extractFields output
Expected format:
CHROM   POS     ID      REF     ALT     INDEL   ANN[*].EFFECT   ANN[*].IMPACT   GEN[*].GT
scaf5   422             C       T       false   downstream_gene_variant,intergenic_region       MODIFIER,MODIFIER        0/1
"""

import argparse

class VariantFile:
    """Defines an extractFields output file."""
    def __init__(self, filename):
        self.filename = filename
        self.indel_effects = {}  # {effect -> count}
        self.snp_effects = {}  # {effect -> count}
        self.genotype_effects = {}  # {genotype -> {effect -> count}}
        self.total_counts = {}
    
    def parse_variants(self):
        """Count the types and effects of variants in the output file."""
        with open(self.filename) as fh:
            fh.readline()
            for var in fh:
                var_info = var.strip().split("\t")
                indel = var_info[5] == "true"  # 6th field describes variant type
                geno = var_info[8]  # 9th field describes variant genotype
                if geno not in self.genotype_effects:
                    self.genotype_effects[geno] = dict()
                effects = var_info[6].split(",")
                effects_clean = set()
                for effect in effects:
                    if "&" in effect:
                        split = set(effect.split("&"))
                        for s in split:
                            effects_clean.add(s)
                    else:
                        effects_clean.add(effect)
                tracker = self.indel_effects if indel else self.snp_effects
                for e in effects_clean:
                    if e in tracker:
                        tracker[e] += 1
                    else:
                        tracker[e] = 1
                    if e in self.genotype_effects[geno]:
                        self.genotype_effects[geno][e] += 1
                    else:
                        self.genotype_effects[geno][e] = 1
                if geno not in self.total_counts:
                    self.total_counts[geno] = 1
                else:
                    self.total_counts[geno] += 1
                var_type = "INDEL" if indel else "SNP"
                if var_type not in self.total_counts:
                    self.total_counts[var_type] = 1
                else:
                    self.total_counts[var_type] += 1
    
    def print_variant_types(self, outfh):
        """Print variant types and effects in tsv format."""
        var_types = {"SNP": self.snp_effects, "INDEL": self.indel_effects}
        for vt in var_types:
            for effect in var_types[vt]:
                print(f"{vt}\t{effect}\t{var_types[vt][effect]}\t{self.filename}", file=outfh)

    def print_variant_genos(self, outfh):
        """Print variant genotypes and effects in tsv format."""
        for gt in self.genotype_effects:
            for effect in self.genotype_effects[gt]:
                print(f"{gt}\t{effect}\t{self.genotype_effects[gt][effect]}\t{self.filename}", file=outfh)
    
    def print_total_counts(self, outfh):
        """Print total counts of SNP, INDEL, genotypes in tsv format."""
        for var_type in sorted(self.total_counts.keys()):
            print(f"{var_type}\t{self.total_counts[var_type]}\t{self.filename}", file=outfh)

def main():
    """Read input files and print summaries to stdout."""
    parser = argparse.ArgumentParser(description="Print counts of variants")
    parser.add_argument("tsv", nargs="+",
            help="SnpSift extractFields tsv output file(s)")
    parser.add_argument("-p", "--outprefix",
            type=str,
            default="summary",
            help="output file prefix [variants]") 
    args = parser.parse_args()
    var_type_filename = args.outprefix + ".types.summary.tsv"
    var_gt_filename = args.outprefix + ".genotypes.summary.tsv"
    total_counts_filename = args.outprefix + ".total_counts.summary.tsv"

    with open(var_type_filename, "w+") as var_types, open(var_gt_filename, "w+") as var_gt, open(total_counts_filename, "w+") as total_counts: 
        print("var_type", "var_effect", "count", "assembly", sep="\t", file=var_types)
        print("var_gt", "var_effect", "count", "assembly", sep="\t", file=var_gt)
        print("var", "count", "assembly", sep="\t", file=total_counts)
        for f in args.tsv:
            cur_file = VariantFile(f)
            cur_file.parse_variants()
            cur_file.print_variant_types(var_types)
            cur_file.print_variant_genos(var_gt)
            cur_file.print_total_counts(total_counts)

if __name__ == "__main__":
    main()
