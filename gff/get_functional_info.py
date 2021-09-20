#!/usr/bin/env python3
"""
Create table of functional info from a GFF file for easier parsing.
"""

import argparse
import re
from gff import *
from orthofinder import get_pfam_desc, get_go_label

class GeneFeatureInfo:
    def __init__(self, info):
        """Given a info field of a GFF file, extract information about the feature."""
        self.info = info
        self.id = GeneFeatureInfo.get_id(self.info)
        self.pfam = GeneFeatureInfo.get_pfam(self.info)
        self.go_terms = GeneFeatureInfo.get_go_terms(self.info)
        self.note = GeneFeatureInfo.get_note(self.info)
        self.parent = GeneFeatureInfo.get_parent(self.info)
    
    @staticmethod
    def get_id(text):
        search = re.search(ID_RE, text)
        return search[1] if search else None
    
    @staticmethod
    def get_parent(text):
        search = re.search(PARENT_RE, text)
        return search[1] if search else None
    
    @staticmethod
    def get_pfam(text):
        pfam = re.findall(PFAM_RE, text)
        return pfam
    
    @staticmethod
    def get_go_terms(text):
        go_terms = re.findall(GO_RE, text)
        return go_terms
    
    @staticmethod
    def get_note(text):
        search = re.search(NOTE_RE, text)
        return search[1] if search else None


def extract_and_lookup(gff, feature, element):
    """Extract info and perform lookup of element ID (GO term of Pfam domain)."""
    element_methods = {
        "pfam": (GeneFeatureInfo.get_pfam, get_pfam_desc),
        "go": (GeneFeatureInfo.get_go_terms, get_go_label),
    }

    info = dict()  # id -> [[accession, description], ...]
    elem_ids = set()
    elem_id_to_name = dict()

    with open(gff, "r") as fh:
        for line in fh:
            if line[0] == "#":
                continue
            line = line.strip().split("\t")
            if line[2] == feature:
                id = GeneFeatureInfo.get_id(line[8])
                element_vals = element_methods[element][0](line[8])

                if len(element_vals):
                    info[id] = []
                    for e in element_vals:
                        elem_ids.add(e)
                        info[id].append([e])  # Only add accession at first
    
    # Look up names of element accession
    for e in elem_ids:
        elem_id_to_name[e] = element_methods[element][1](e)
    
    # Add element name to dict
    for id in info:
        for i, elem in enumerate(info[id]):
            elem_name = elem_id_to_name[elem[0]]
            info[id][i].append(elem_name)
    
    return info


def extract(gff, feature, element):
    """Extract a given element from a GFF feature."""
    element_methods = {
        "note": GeneFeatureInfo.get_note,
        "parent": GeneFeatureInfo.get_parent
    }
    
    info = dict()

    with open(gff, "r") as fh:
        for line in fh:
            if line[0] == "#":  # Skip header lines
                continue
            line = line.strip().split("\t")
            if line[2] == feature:
                id = GeneFeatureInfo.get_id(line[8])
                element_val = element_methods[element](line[8])
                if not element_val is None:
                    info[id] = [element_val]
    
    return info


def parse_args():
    """Parse the command line arguments."""
    parser = argparse.ArgumentParser(description="Extract GFF functional info elements into a table")
    parser.add_argument("gff",
                        type=str,
                        help="MAKER GFF file with functional info")
    parser.add_argument("-f", "--feature",
                        type=str,
                        default="mRNA",
                        help="Feature type to be extracted from GFF")
    parser.add_argument("-e", "--element",
                        type=str,
                        choices=["pfam", "go", "note", "parent"],
                        default="note",
                        help="Element to extract from information field (along with ID) [note]")
    
    return parser.parse_args()

def main():
    args = parse_args()
    if args.element == "pfam" or args.element == "go":
        info = extract_and_lookup(args.gff, args.feature, args.element)
        print(args.feature, args.element, "description", sep="\t")
        for id in info:
            for elem in info[id]:
                print(id, *elem, sep="\t")
    else:
        info = extract(args.gff, args.feature, args.element)
        print(args.feature, args.element)
        for i in info:
            print(i, *info[i], sep="\t")

if __name__ == "__main__":
    main()
