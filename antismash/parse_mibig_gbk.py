#!/usr/bin/env python3

import sys
from Bio import SeqIO

def parse_mibig_gbk(genbank):
    """
    Parse a MIBiG cluster GenBank file. Returns a dictionary
    describing the cluster: protein id -> gene info
    """
    cluster = {}
    with open(genbank, "r") as fh:
        for record in SeqIO.parse(fh, "gb"):
            for feature in record.features:
                if feature.type == "CDS" and "protein_id" in feature.qualifiers:
                    protein_id = feature.qualifiers["protein_id"][0]
                    gene_info = {}
                    gene_info["product"] = feature.qualifiers["product"][0]
                    if "gene_kind" in feature.qualifiers:
                        gene_info["gene_kind"] = feature.qualifiers["gene_kind"][0]
                    if "gene_functions" in feature.qualifiers:
                        gene_info["gene_functions"] = feature.qualifiers["gene_functions"]
                    cluster[protein_id] = gene_info
    return cluster

if __name__ == "__main__":
    gb = sys.argv[1]
    print(parse_mibig_gbk(gb))
