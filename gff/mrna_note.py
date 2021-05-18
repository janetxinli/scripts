#!/usr/bin/env python3

import re
from gff import ID_RE, NOTE_RE

def mrna_note(gff):
    """
    Get functional info (GO terms and Pfam domains) from
    a MAKER gff file.
    """
    notes = {}  # mrna_id -> note

    with open(gff, "r") as fh:
        for line in fh:
            if not line.startswith("#"):
                line = line.strip().split("\t")
                if line[2] == "mRNA":
                    info = line[8]
                    mrna_id = re.search(ID_RE, info)[1]
                    note = re.search(NOTE_RE, info)[1]
                    notes[mrna_id] = note

    return notes
