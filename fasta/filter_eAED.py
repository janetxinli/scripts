#!/usr/bin/env python3
"""
Filter sequences in a fasta file by eAED score.
Expected scored fasta header format (MAKER):
>augustus_masked-scaf926-processed-gene-0.3-mRNA-1 protein AED:1.00 eAED:1.00 QI:0|0|0|0|1|1|5|0|1308
Expected fasta header format (GAG):
>protein|genemark-scaf0-processed-gene-0.15-mRNA-1 ID=genemark-scaf0-processed-gene-0.15-mRNA-1|Parent=genemark-scaf0-processed-gene-0.15|Name=
"""

import sys
import argparse

def maker_id(fasta_header):
    """Returns the sequence ID in a MAKER fasta header."""
    return fasta_header.split(" ")[0].strip(">")


def gag_id(fasta_header):
    """Returns the sequence ID in a GAG fasta header."""
    return fasta_header.split(" ")[0].split("|")[1]


def eAED_score(fasta_header):
    """
    Returns the eAED score in a fasta header. Raises RuntimeError if header
    doesn't include eAED score.
    """
    if "eAED" in fasta_header:
        return float(fasta_header.strip().split(" ")[3].split(":")[1])
    else:
        raise RuntimeError("Fasta header does not contain eAED score")


def load_scores(scored_fasta):
    """Load eAED scores and headers from MAKER output."""
    header_scores = {}
    with open(scored_fasta, "r") as fh:
        for line in fh:
            #print(line)
            if line[0] == ">":  # Header
                hid = maker_id(line)
                #print(hid)
                header_scores[maker_id(line)] = eAED_score(line)
    return header_scores


def filter_seqs(fasta, filter, header_scores=None):
    """Filter sequences with eAED equal or greater than a given value."""
    if header_scores is not None:
        htype = "gag"
    else:
        htype = "maker"
    #print("header score type is {}".format(htype))
    with open(fasta, "r") as fh:
        for line in fh:
            if line[0] == ">":
                if htype == "gag":
                    header = gag_id(line)
                    print_seq = header_scores[header] < filter
                else:
                    try:
                        score = eAED_score(line)
                    except RuntimeError:
                        print("error: input fasta file {} does not contain eAED scores".format(fasta), file=sys.stderr)
                        sys.exit(1)
                    print_seq = score < filter
                if print_seq:
                    print(line.strip())
            else:
                if print_seq:
                    print(line.strip())


def get_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Filter out sequences by eAED")
    parser.add_argument("input",
                        nargs="?",
                        type=str,
                        default="/dev/stdin",
                        help="Input fasta file [stdin]")
    parser.add_argument("-f", "--filter",
                        type=float,
                        default=1.0,
                        help="Filter sequences with eAED equal to or greater than this value [1]")
    parser.add_argument("-s", "--scored_seqs",
                        type=str,
                        default=None,
                        help="MAKER Fasta file with scored sequences and matching headers (alternative to filtering input directly)")
    return parser.parse_args()


def main():
    args = get_args()
    if args.scored_seqs:
        #print("scored seqs")
        header_scores = load_scores(args.scored_seqs)
        #print(header_scores)
        filter_seqs(args.input, args.filter, header_scores=header_scores)
    else:
        filter_seqs(args.input, args.filter)


if __name__ == "__main__":
    main()
