#!/usr/bin/env python3

import os
import sys
import argparse

def convert_file(in_fh, outfile, column_headers, no_header):
    """Takes as input a tsv file handle and formats the input for Jira."""
    if outfile == None:
        out = sys.stdout
    else:
        out = open(outfile, "w+")
    header_split = in_fh.readline().strip().split("\t")
    if no_header:
        jira_header = "|" + "|".join(header_split) + "|"
    else:
        jira_header = "||" + "||".join(header_split) + "||"
    print(jira_header, file=out)
    for line in in_fh:
        jira_line = ""
        if column_headers:
            jira_line += "|"
        line_split = line.strip().split("\t")
        jira_line += "|" + "|".join(line_split) + "|"
        print(jira_line, file=out)
    out.close()

def main():
    parser = argparse.ArgumentParser(description="Convert a tsv file to Jira table format")
    parser.add_argument("-i", "--input", 
                        type=str,
                        nargs="?",
                        default=None,
                        help="Input tsv file to be converted [stdin]")
    parser.add_argument("-o", "--output",
                        type=str, default=None,
                        help="Output file [stdout]")
    parser.add_argument("-c", "--columnheaders",
                        default=False, action="store_true",
                        help="Create column headers [false]")
    parser.add_argument("-n", "--noheader",
                        action="store_true",
                        default=False,
                        help="Do not print a header row")
    args = parser.parse_args()

    if args.input == None:
        args.input = "/dev/stdin"

    with open(args.input) as fh:
        if not fh.isatty():
            convert_file(fh, args.output,
                args.columnheaders, args.noheader)
        else:
            print("to_jira.py: error: missing input, pass a file as an argument or from stdin")
            sys.exit(1)

if __name__ == "__main__":
    main()