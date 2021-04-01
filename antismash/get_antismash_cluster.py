#!/usr/bin/env python3

import sys
import argparse
import json
from antismash import parse_mibig_gbk

def print_cluster_info(ref, hits, print_name=False):
    """Summarize AntiSMASH hits and misses."""
    matches = set()
    for pairing in hits:
        prot_id = pairing[2]["locus_tag"] 
        if prot_id in ref:
            matches.add(pairing[2]["locus_tag"])
        else:
            print("error: AntiSMASH output cluster {} not found in ref".format(prot_id),
                file=sys.stderr)
            sys.exit(1)
        match_id = pairing[0]
        product = ref[prot_id]["product"] if "product" in ref[prot_id] else "NA"
        kind = ref[prot_id]["gene_kind"] if "gene_kind" in ref[prot_id] else "NA"
        fn = ",".join(ref[prot_id]["gene_functions"]) if "gene_functions" in ref[prot_id] else "NA"
        hit = 1
        cov = pairing[2]["perc_coverage"]
        identity = pairing[2]["perc_ident"]
        name = print_name if print_name else ""
        print(prot_id, match_id, product, kind,
              fn, hit, cov, identity, name, sep="\t")
    
    remaining = [i for i in ref if i not in matches]
    for r in remaining:
        product = ref[r]["product"] if "product" in ref[r] else "NA"
        kind = ref[r]["gene_kind"] if "gene_kind" in ref[r] else "NA"
        fn = ",".join(ref[r]["gene_functions"]) if "gene_functions" in ref[r] else "NA"
        name = print_name if print_name else ""
        print(r, "NA", product, kind, fn, 0, "NA", "NA", name, sep="\t")


def find_cluster(json_file, cluster_id):
    """
    Parse AntiSMASH JSON file and return cluster of interest.
    Return None if cluster of interest is not found.
    """
    with open(json_file, "r") as fh:
        content = json.load(fh)

        for i, record in enumerate(content["records"]):
            for feature in record["features"]:
                if feature["type"] == "region":
                    region_no = int(feature["qualifiers"]["region_number"][0])
                    knowncluster = record["modules"]["antismash.modules.clusterblast"]["knowncluster"]\
                        ["results"][region_no-1]
                    if knowncluster["total_hits"] > 0 and knowncluster["ranking"][0][0]["accession"] == cluster_id:
                        return knowncluster["ranking"][0][1]["pairings"]
        
        return None

def read_fof(fof):
    """Returns a list of files in fof."""
    files = []
    with open(fof, "r") as fh:
        for line in fh:
            files.append(line.strip())
    
    return files

def run(json_file, ref, cluster_id, print_name=False):
    """Run analysis for a single JSON file."""
    antismash_cluster = find_cluster(json_file, cluster_id)
    print_cluster_info(ref, antismash_cluster, print_name)

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Pull out information about a cluster from "
                                     "AntiSMASH JSON output",
                                     usage="usage: get_antismash_cluster.py [-h] [-j [JSON] or -f [IN]] genbank cluster")
    parser.add_argument("genbank",
                        type=str,
                        help="AntiSMASH cluster GenBank file")
    parser.add_argument("cluster",
                        type=str,
                        help="ID of cluster of interest")
    parser.add_argument("-j", "--json",
                        type=str,
                        nargs="?",
                        help="AntiSMASH JSON output. If this option is selected, "
                             "-f must not be")
    parser.add_argument("-f", "--fof",
                        type=str,
                        nargs="?",
                        help="Input file containing several AntiSMASH JSON files to parse. "
                             "Files must be line-separated. If this option is selected, "
                             "-j must not be")
    return parser.parse_args()


def main():
    args = parse_args()
    if args.json and args.fof:
        print("error: provide either a single AntiSMASH JSON output file or a file of files")
    
    cluster_ref = parse_mibig_gbk(args.genbank)

    if args.json:
        print("protein_id\tmatch_id\tproduct\tgene_kind\tgene_fn\thit\tperc_cov\tperc_identity")
        run(args.json, cluster_ref, args.cluster)
    elif args.fof:
        files = read_fof(args.fof)
        print("protein_id\tmatch_id\tproduct\tgene_kind\tgene_fn\thit\tperc_cov\tperc_identity\tname")
        for f in files:
            run(f, cluster_ref, args.cluster, print_name=f)

    

if __name__ == "__main__":
    main()
