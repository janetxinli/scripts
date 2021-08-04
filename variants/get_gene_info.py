#!/usr/bin/env python3
"""
Get gene names and AED/eAED score for snpEff gene output.
Also cleans file to have one impact/effect per line.
"""

import argparse
import os
import re
import sys
import pandas as pd
from collections import namedtuple

Gene = namedtuple("Gene", ["name", "aed", "eaed"])

def load_genes(fasta):
    """
    Reads gene names and AED/eAED score from a MAKER fasta 
    file. Returns a dictionary of gene_id -> Gene (name, aed, eaed)
    """
    genes = dict()
    
    with open(fasta, "r") as fh:
        for line in fh:
            if line.startswith(">"):  # fasta header
                line = line.strip()
                gene_id = re.match(">([^\s]+)", line)[1]
                name = re.search("Name:\"(.+)\"", line)[1]
                aed = re.search("AED:(\d*\.\d*)", line)[1]
                eaed = re.search("eAED:(\d*\.\d*)", line)[1]
                genes[gene_id] = Gene(name, aed, eaed)
    
    return genes


def get_gene_info(x, gene_info):
    """Get gene info for a row in snpEff gene variants file."""
    gene = gene_info[x.TranscriptId]
    return pd.Series([gene.name, gene.aed, gene.eaed])


def add_gene_info(vargenes, gene_info, outfile):
    """Add gene info and pivot table, then write to outfile."""
    if outfile is None:
        outfile = sys.stdout
    
    df = pd.read_csv(vargenes, sep="\t", skiprows=1)
    df[["maker_gene_name", "AED", "eAED"]] = df.apply(lambda x: get_gene_info(x, gene_info), axis=1)

    # Pivot variants_impact_* columns
    impact_cols = [i for i in df.columns if i.startswith("variants_impact")]
    id_vars = list(set(df.columns) - set(impact_cols))
    df = pd.melt(df, id_vars=id_vars, value_vars=impact_cols, var_name="variant_impact_type", value_name="variant_impact_count")
    df["variant_impact_type"] = df["variant_impact_type"].map(lambda x: x.replace("variants_impact_", ""))
    
    # Pivot variants_effect_* columns
    effect_cols = [i for i in df.columns if i.startswith("variants_effect")]
    id_vars = list(set(df.columns) - set(effect_cols))
    df = pd.melt(df, id_vars=id_vars, value_vars=effect_cols, var_name="variant_effect_type", value_name="variant_effect_count")
    df["variant_effect_type"] = df["variant_effect_type"].map(lambda x: x.replace("variants_effect_", ""))

    df = df.loc[(df["variant_impact_count"] != 0) & (df["variant_effect_count"] != 0)]
    
    df[[
        "#GeneName", "GeneId", "TranscriptId", "BioType", "variant_impact_type",
        "variant_impact_count", "variant_effect_type", "variant_effect_count",
        "AED", "eAED", "maker_gene_name"
        ]].to_csv(outfile, sep="\t", index=False)


def parse_args():
    parser = argparse.ArgumentParser(description="Get gene names and annotation scores for snpEff gene output")
    
    parser.add_argument("vargenes",
                        type=str,
                        help="snpEff gene summary output txt file")
    parser.add_argument("fasta",
                        type=str,
                        help="MAKER Fasta file for variant calling reference "
                             "containing gene names and scores in headers")
    parser.add_argument("-o", "--outfile",
                        type=str,
                        default=None,
                        help="Output file name [stdout]")
    
    return parser.parse_args()


def main():
    args = parse_args()

    if not args.outfile is None:
        if os.path.exists(args.outfile):
            print(f"Error: output file {args.outfile} already exists. Exiting.")
    
    gene_info = load_genes(args.fasta)
    add_gene_info(args.vargenes, gene_info, args.outfile)


if __name__ == "__main__":
    main()
