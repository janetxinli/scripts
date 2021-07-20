#!/usr/bin/env python3
"""Get mRNA functional info (GO and Pfam) from a MAKER gff file."""

import re
from gff import ID_RE, GO_RE, PFAM_RE

def functional_info(gff, feature="mRNA"):
    """
    Get functional info (GO terms and Pfam domains) from
    a MAKER gff file.
    """
    func = {}  # feature_id -> ([go terms], [pfam domains])

    with open(gff, "r") as fh:
        for line in fh:
            if not line.startswith("#"):
                line = line.strip().split("\t")
                if line[2] == feature:
                    info = line[8]
                    id = re.search(ID_RE, info)[1]
                    go_terms = re.findall(GO_RE, info)
                    pfam_doms = re.findall(PFAM_RE, info)
                    func[id] = (go_terms, pfam_doms)
    
    return func