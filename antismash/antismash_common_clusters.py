#!/usr/bin/env python3
"""
Find common and group-specific antiSMASH clusters. Requires summary .tsv file (output
from parse_antismash_json.py).

For group-specific antiSMASH clusters, provide a .txt file defining groups in the
following format:
group_name\tsample1,sample2,sample3

Sample names must be present in summary .tsv file.
"""

import sys
import argparse

def find_common(known_clusters, outfile):
    """Find common clusters in dictionary of sets (intersection)."""
    with open(outfile, "w") as outfh:
        common_clusters = [i for i in set.intersection(*[i for i in known_clusters.values()])]
        print("accession\tdescription\tknownclusterblast", file=outfh)
        if len(common_clusters) > 0:
            for cluster in common_clusters:
                print(cluster[0], cluster[1], cluster[2], sep="\t", file=outfh)


def find_group_specific(known_clusters, groups, outfile):
    """Find clusters that are unique to a group."""
    with open(outfile, "w") as outfh:
        print("group\taccession\tdescription\tknownclusterblast", file=outfh)
        group_names = [i for i in groups.keys()]
        for n in group_names:
            group_clusters = []
            non_group = [g for g in group_names if g != n]
            non_group_clust = []
            for non in non_group:
                for sample in groups[non]:
                    non_group_clust.append(known_clusters[sample])
            non_group_clust_union = set.union(*non_group_clust)
            for sample in groups[n]:
                sample_specific_clust = set()
                sample_clust = known_clusters[sample]
                for s in sample_clust:
                    if s not in non_group_clust_union:
                        sample_specific_clust.add(s)
                    group_clusters.append(sample_specific_clust)
            group_specific_clusters = set.intersection(*group_clusters)
            for cluster in group_specific_clusters:
                print(n, cluster[0], cluster[1], cluster[2], sep="\t", file=outfh)

def find_clusters(filename):
    """
    Read antiSMASH tsv summary and load clusters into a dictionary.
    Returns dict{sample -> {(accession, description, knowncluster)}}
    """
    known_clusters = dict()

    with open(filename, "r") as fh:
        fh.readline()
        for line in fh:
            line = line.strip().split("\t")
            cur_sample = line[0]
            if len(line) > 7:
                known = line[5]
                acc = line[6]
                desc = line[7].replace(" ", "_")
                clust_type = line[8]
                if acc != "":
                    if cur_sample not in known_clusters:
                        known_clusters[cur_sample] = set()
                    known_clusters[cur_sample].add((acc, desc, known))

    return known_clusters

def parse_groups(txt):
    """
    Parse text file containing groups. Returns a dictionary mapping
    group -> [samples].
    """
    groups = {}
    with open(txt, "r") as fh:
        for line in fh:
            line = line.strip().split("\t")
            if len(line) != 2:
                print("error: group text file {} formatted incorrectly".format(txt))
                sys.exit(1)
            groups[line[0]] = line[1].split(",")
    
    return groups

def parse_args():
    """Parse the command line arguments."""
    parser = argparse.ArgumentParser(description="Find common and group-specific "
                                     "clusters in an antismash tsv summary file")
    parser.add_argument("summary",
                        type=str,
                        help="Summary tsv file")
    parser.add_argument("-p", "--out_prefix",
                        type=str,
                        default=None,
                        help="Output file prefix")
    parser.add_argument("-g", "--groups",
                        type=str,
                        default=None,
                        help="Text file containing groups (leave blank to only search "
                             "common clusters)")
    return parser.parse_args()

def main():
    args = parse_args()

    if args.out_prefix is None:
        args.out_prefix = ""
    else:
        args.out_prefix = args.out_prefix + "."

    # Load all clusters
    clusters = find_clusters(args.summary)

    # Find common clusters
    common_clust_file = args.out_prefix + "common_clusters.tsv"
    common_clusters = find_common(clusters, common_clust_file)

    # Find group specific clusters
    if not args.groups is None:
        sample_groups = parse_groups(args.groups)
        group_clust_file = args.out_prefix + "group_specific_clusters.tsv"
        find_group_specific(clusters, sample_groups, group_clust_file)


if __name__ == "__main__":
    main()
