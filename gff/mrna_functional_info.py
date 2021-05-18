#!/usr/bin/env python3
"""Get functional info (GO and Pfam) from a MAKER gff file."""

import re
from gff import ID_RE, GO_RE, PFAM_RE

def mrna_functional_info(gff):
    """
    Get functional info (GO terms and Pfam domains) from
    a MAKER gff file.
    """
    func = {}  # mrna_id -> ([go terms], [pfam domains])

    with open(gff, "r") as fh:
        for line in fh:
            if not line.startswith("#"):
                line = line.strip().split("\t")
                if line[2] == "mRNA":
                    info = line[8]
                    mrna_id = re.search(ID_RE, info)[1]
                    go_terms = re.findall(GO_RE, info)
                    pfam_doms = re.findall(PFAM_RE, info)
                    func[mrna_id] = (go_terms, pfam_doms)
    
    return func