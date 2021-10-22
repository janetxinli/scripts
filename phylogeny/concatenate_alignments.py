#!/usr/bin/env python3

import argparse
import sys
from Bio import Align, AlignIO

def combine_fasta_aln(filenames, ids):
    """
    Combine a list of FASTA alignments into a single MultipleSeqAlignment object.
    Returns a tuple of MultipleSeqAlignment, alignment_lengths, filenames
    """
    alignment_lengths = []
    files = []

    # Combine alignments
    with open(filenames, "r") as fh:
        for i, line in enumerate(fh):
            line = line.strip()
            files.append(line)
            if i == 0:
                alns = AlignIO.read(line, "fasta")
                
                # Check that correct number of IDs were provided
                if len(alns) != len(ids):
                    print(f"{sys.argv[0]}: error: number of IDs and alignments do not match")
                    sys.exit(1)
                
                alignment_lengths.append(alns.get_alignment_length())
            else:
                cur_aln = AlignIO.read(line, "fasta")
                alns += cur_aln
                alignment_lengths.append(cur_aln.get_alignment_length())
    
    # Add IDs
    for i, id in enumerate(ids):
        alns[i].id = id
    
    return alns, alignment_lengths, files

def print_partitions(alignment_lengths, filenames, partition):
    if len(alignment_lengths) != len(filenames):
        raise ValueError("Alignment lengths must equal number of alignment files")
    
    with open(partition, "w+") as fh:
        end = 0
        for i, filename in enumerate(filenames):
            start = end + 1
            end = start + alignment_lengths[i] - 1
            print(f"DNA, {filename}={start}-{end}", file=fh, flush=True)

def parse_args():
    parser = argparse.ArgumentParser(description="Combine several FASTA DNA alignments into relaxed "
                                     "PHYLIP format and create a partition file")
    parser.add_argument("filenames",
                        type=str,
                        help="Newline-separated list of FASTA files to be concatenated")
    parser.add_argument("phylip",
                        type=str,
                        help="Name of output PHYLIP file")
    parser.add_argument("partition",
                        type=str,
                        help="Name of output partition file")
    parser.add_argument("-i", "--ids",
                        nargs="+",
                        required=True,
                        help="Space-separated list of IDs for output PHYLIP file. Must be "
                             "provided in the same order as FASTA alignments")

    return parser.parse_args()


def main():
    args = parse_args()
    concat_aln, aln_lengths, filenames = combine_fasta_aln(args.filenames, args.ids)
    
    # Write concatenated alignment to file
    AlignIO.write(concat_aln, args.phylip, "phylip-relaxed")

    # Write partitions to file
    print_partitions(aln_lengths, filenames, args.partition)

if __name__ == "__main__":
    main()
