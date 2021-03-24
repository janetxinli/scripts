#!/usr/bin/env python3

import sys
import argparse
import json
from misc import parse_mibig_gbk

def print_cluster_info(ref, hits):
    """Summarize AntiSMASH hits and misses."""
    print("protein_id\tmatch_id\tproduct\tgene_kind\tgene_fn\thit\tperc_cov\tperc_identity")
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
        product = ref[prot_id]["product"]
        kind = ref[prot_id]["gene_kind"] if "gene_kind" in ref[prot_id] else "NA"
        fn = ",".join(ref[prot_id]["gene_functions"]) if "gene_functions" in ref[prot_id] else "NA"
        hit = 1
        cov = pairing[2]["perc_coverage"]
        identity = pairing[2]["perc_ident"]
        print(prot_id, match_id, product, kind,
              fn, hit, cov, identity, sep="\t")
    
    remaining = [i for i in ref if i not in matches]
    for r in remaining:
        product = ref[r]["product"]
        kind = ref[r]["gene_kind"] if "gene_kind" in ref[r] else "NA"
        fn = ref[r]["gene_functions"] if "gene_functions" in ref[r] else "NA"
        print(r, "NA", product, kind, fn, 0, "NA", "NA", sep="\t")


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

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Pull out information about a cluster from "
                                     "AntiSMASH JSON output")
    parser.add_argument("json",
                        type=str,
                        help="AntiSMASH JSON output file")
    parser.add_argument("genbank",
                        type=str,
                        help="AntiSMASH cluster GenBank file")
    parser.add_argument("cluster",
                        type=str,
                        help="ID of cluster of interest")
    return parser.parse_args()


def main():
    args = parse_args()
    cluster_ref = parse_mibig_gbk(args.genbank)
    antismash_cluster = find_cluster(args.json, args.cluster)
    if antismash_cluster is None:
        print("error: cluster not found in AntiSMASH JSON {}".format(args.json), file=sys.stderr)
        sys.exit(1)
    print_cluster_info(cluster_ref, antismash_cluster)
    

if __name__ == "__main__":
    main()