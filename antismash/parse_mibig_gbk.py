#!/usr/bin/env python3

import sys
from Bio import SeqIO

def get_cluster_record(gb, id):
    """
    Get specified cluster record from GenBank file from its ID.
    Returns None if no record is present with the name id in the
    GenBank file.
    """
    full = SeqIO.parse(gb, "gb")
    for record in full:
        if record.name == id:
            return record
    
    return None


def get_cluster_info(record):
    """
    Parse cluster info from a GenBank record. Returns a dictionary
    describing the cluster: protein id -> gene info
    """
    cluster = {}
    for feature in record.features:
        if feature.type == "CDS":
            gene_info = {}
            if "protein_id" in feature.qualifiers:
                protein_id = feature.qualifiers["protein_id"][0]
            elif "locus_tag" in feature.qualifiers:
                protein_id = feature.qualifiers["locus_tag"][0]
            elif "gene" in feature.qualifiers:
                protein_id = feature.qualifiers["gene"][0]
            if "product" in feature.qualifiers:
                gene_info["product"] = feature.qualifiers["product"][0]
            if "gene_kind" in feature.qualifiers:
                gene_info["gene_kind"] = feature.qualifiers["gene_kind"][0]
            if "gene_functions" in feature.qualifiers:
                gene_info["gene_functions"] = feature.qualifiers["gene_functions"]                        
            cluster[protein_id] = gene_info
    
    return cluster