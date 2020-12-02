#!/usr/bin/env python3
"""
Count variant types (SNP, indel), zygosity and effects in SnpSift extractFields output
Expected format:
CHROM   POS     ID      REF     ALT     INDEL   ANN[*].EFFECT   ANN[*].IMPACT   GEN[*].GT
scaf5   422             C       T       false   downstream_gene_variant,intergenic_region       MODIFIER,MODIFIER        0/1
"""

import sys
import argparse

class VariantFile:
    """Defines a extractFields output file."""
    def __init__(self, filename):
        self.filename = filename
        self.indel_effects = {}  # effect -> (total_count, heterozygous, homozygous)
        self.snp_effects = {}  # effect -> (total_count, heterozygous, homozygous)
    
    def parse_variants(self):
        """Count the types and effects of variants in the output file."""
        with open(self.filename) as fh:
            fh.readline()
            for var in fh:
                var_info = var.strip().split("\t")
                indel = var_info[5] == "true"  # 6th field describes variant type
                geno = 2 if var_info[8] == "1/1" else 1  # homozygous = 2; heterozygous = 1
                effects = var_info[6].split(",")
                tracker = self.indel_effects if indel else self.snp_effects
                for e in effects:
                    if e in tracker:
                        tracker[e][0] += 1
                        tracker[e][geno] += 1
                    else:
                        if geno == 1:  # heterozyogus
                            tracker[e] = [1, 1, 0]
                        else:  # homozygous
                            tracker[e] = [1, 0, 1]
    
    def print_variants(self):
        """
        Print variant types and effects in tsv format to stdout, eg:
        var_type\tvar_effect\tcount\tassembly\n
        """
        self.parse_variants()
        var_types = {"SNP": self.snp_effects, "INDEL": self.indel_effects}
        for vt in var_types:
            for effect in var_types[vt]:
                print(f"{vt}\t{effect}\t{var_types[vt][effect][0]}\t{self.filename}")    


def main():
    """Read input files and print summaries to stdout."""
    if len(sys.argv) < 2:
        print("Variant files must be provided as arguments")
        sys.exit(1)
    variant_files = sys.argv[1:]
    print("var_type", "var_effect", "count", "assembly", sep="\t")
    for f in variant_files:
        VariantFile(f).print_variants()


if __name__ == "__main__":
    main()
