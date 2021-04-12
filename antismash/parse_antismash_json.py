#!/usr/bin/env python3
"""
Parse antiSMASH JSON output (regions and KnownClusterBlast best hits) and print to TSV.
Takes a single JSON or a file of files in the following format:
name1\tfile1
name2\tfile2
name4\tfile3
or
file1
file2
file3
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
    
    def get_location(self, idx, region_no):
        """
        Returns the location of a region along a scaffold for a given
        record index as a tuple of (start, stop).
        """
        for f in self.data[idx]["features"]:
            if f["type"] == "region" and int(f["qualifiers"]["region_number"][0]) == region_no:
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
            hits = sub_dict["ranking"][0][1]["hits"]
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
            hits = sub_dict[0][1]["hits"]
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
                    start, end = self.get_location(i, region_no)
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


def get_file_info(filename):
    """
    Reads a file of names and json files and returns a
    tuple of names, files
    """
    names = []
    files = []
    with open(filename, "r") as fh:
        for line in fh:
            line = line.strip().split("\t")
            if len(line) == 2:
                names.append(line[0])
                files.append(line[1])
            elif len(line) == 1:  # Assuming no name is provided
                names.append(None)
                files.append(line[0])
            else:
                print("error: file of files {} formatted incorrectly".format(filename))
                sys.exit(1)
    
    return names, files


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Parse antiSMASH json output and summarize regions")
    parser.add_argument("-j", "--json",
                        type=str,
                        default=None,
                        help="antiSMASH JSON file to be analysed")
    parser.add_argument("-n", "--name",
                        type=str,
                        default=None,
                        help="Matching name for JSON file given")
    parser.add_argument("-f", "--fof",
                        type=str,
                        default=None,
                        help="File of antismash JSON files and optional names")
    parser.add_argument("-o", "--outfile",
                        default=None,
                        help="Output file for printing results [stdout]")
    return parser.parse_args()


def main():
    args = parse_args()

    if args.json and args.fof:
        print("error: provide either a single json file or file of files")
        sys.exit(1)
    
    if args.fof:
        names, files = get_file_info(args.fof)
    else:
        names = [args.name]
        files = [args.json]
    
    if args.outfile is None:
        outfh = sys.stdout
    else:
        outfh = open(args.outfile, "w+")
    
    print("sample\tscaffold\tregion_no\tstart\tend\tknown_cluster\taccession\t"
          "description\tcluster_type\thits", file=outfh)
    for i, json_file in enumerate(files):
        content = AntismashJson(json_file, names[i])
        content.parse_and_print(args.outfile)
        content.close()
    outfh.close()

if __name__ == "__main__":
    main()
