#!/usr/bin/env python3
"""
Parse antiSMASH JSON output (regions and KnownClusterBlast best hits) and print to TSV.
Output fields:
sample: file name/sample named
scaffold: scaffold ID that region belongs to
region_no: region number
start: start position of region on scaffold
end: end position of region on scaffold
best_hit_acc: MiBIG accession best KnownClusterBlast hit (if available)
best_hit_desc: MiBIG description of best KnownCluseterBlast hit (if available)
best_hit_sim: Similarity (%) of cluster to best KnownClusterBlast hit (if available) 
"""

__author__ = "janetxinli"

import os
import json
import sys
import argparse

class AntismashJson:
    def __init__(self, filename, name):
        self.fh = open(filename)
        self.data = json.load(self.fh)["records"]
        self.basename = os.path.basename(filename)
        self.name = name
    
    def get_scaffold(self, idx):
        """Returns the scaffold ID for a given record index."""
        return self.data[idx]["id"]
    
    def get_location(self, idx):
        """
        Returns the location of a region along a scaffold for a given
        record index as a tuple of (start, stop).
        """
        for f in self.data[idx]["features"]:
            if f["type"] == "region":
                loc = f["location"].strip("[]").split(":")
                return loc[0], loc[1]

    def get_type(self, idx, region_no):
        """
        Returns the cluster type for a given record. If multiple types are 
        defined, they are separated by "-".
        """
        for f in self.data[idx]["features"]:
            if f["type"] == "region" and f["qualifiers"]["region_number"][0] == str(region_no):
                clust_type = ";".join(f["qualifiers"]["product"])
                return clust_type
    
    def get_known_cluster(self, idx, region_no):
        """
        Returns the best hit (KnownClusterBlast) for a given record index
        as a tuple of (accession, description, cluster type, similarity). If no hits are found,
        returns a tuple of empty strings.
        """
        sub_dict = self.data[idx]["modules"]["antismash.modules.clusterblast"] \
            ["knowncluster"]["results"][region_no - 1]
        if sub_dict["total_hits"] > 0:
            acc = sub_dict["ranking"][0][0]["accession"]
            desc = sub_dict["ranking"][0][0]["description"]
            #num_prots = len(sub_dict["ranking"][0][0]["proteins"])
            hits = sub_dict["ranking"][0][1]["hits"]
            #similarity = int((hits / num_prots) * 100)
            return acc, desc, hits
        else:
            return "", "", ""
    
    def get_best_hit(self, idx, region_no):
        """
        Returns the top hit (ClusterBlast) for a given record index
        as a tuple of (accession, description, cluster type, similarity). If no hits are found,
        returns a tuple of empty strings.
        """
        sub_dict = self.data[idx]["modules"]["antismash.modules.clusterblast"] \
            ["general"]["results"][region_no - 1]["ranking"]
        if len(sub_dict) > 0:
            acc = sub_dict[0][0]["accession"] + "_" + sub_dict[0][0]["cluster_label"]
            desc = sub_dict[0][0]["description"]
            #num_prots = len(sub_dict[0][0]["proteins"])
            hits = sub_dict[0][1]["hits"]
            #similarity = int((hits / num_prots) * 100)
            return acc, desc, hits
        else:
            return "", "", ""

    def parse_and_print(self, outfh):
        """Parses JSON and prints to file handle."""
        print_name = self.name if (self.name is not None) else self.basename
        
        for i, record in enumerate(self.data):
            for f in record["features"]:
                if f["type"] == "region":
                    region_no = int(f["qualifiers"]["region_number"][0])
                    scaf = self.get_scaffold(i)
                    start, end = self.get_location(i)
                    clust_type = self.get_type(i, region_no)
                    acc, desc, sim = self.get_known_cluster(i, region_no)
                    if (acc, desc, sim) != ("", "", ""):
                        print(print_name, scaf, region_no, start, end, "1",
                              acc, desc, clust_type, sim, sep="\t", file=outfh)
                    else:
                        acc, desc, sim = self.get_best_hit(i, region_no)
                        print(print_name, scaf, region_no, start, end, "0",
                              acc, desc, clust_type, sim, sep="\t", file=outfh)
    
    def close(self):
        """Close input file."""
        self.fh.close()


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Parse antiSMASH json output and summarize regions")
    parser.add_argument("files",
                        type=str,
                        help="Files to be summarized (separated by spaces and surrounded by quotation "
                             "marks if more than one, e.g. 'file1.json file2.json file2.json')")
    parser.add_argument("-n", "--names",
                        type=str,
                        default=None,
                        help="Matching names for files given (also separated by spaces and surrounded "
                             "by quotation marks) (optional)")
    parser.add_argument("-o", "--outfile",
                        default=None,
                        help="Output file for printing results [stdout]")
    return parser.parse_args()


def main():
    args = parse_args()
    args.files = args.files.split(" ")

    if args.names is not None:
        args.names = args.names.split(" ")
        if len(args.files) != len(args.names):
            print("error: number of files and names must be equal")
            sys.exit(1)
    else:
        args.names = [None] * len(args.files)
    
    if args.outfile is None:
        outfh = sys.stdout
    else:
        outfh = open(args.outfile, "w+")
    
    print("sample\tscaffold\tregion_no\tstart\tend\tknown_cluster\taccession\t"
          "description\tcluster_type\thits", file=outfh)
    for i, json_file in enumerate(args.files):
        content = AntismashJson(json_file, args.names[i])
        content.parse_and_print(args.outfile)
        content.close()
    outfh.close()

if __name__ == "__main__":
    main()
