#!/usr/bin/env python3

import sys
import argparse

def convert_file(infile, outfile, column_headers, no_header):
    """Takes as input a tsv file handle and formats the input for Jira."""
    with open(infile, "r") as fh:
        if outfile == None:
            out = sys.stdout
        else:
            out = open(outfile, "w")
        header_split = fh.readline().strip().split("\t")
        if no_header:
            jira_header = "|" + "|".join(header_split) + "|"
        else:
            jira_header = "||" + "||".join(header_split) + "||"
        print(jira_header, file=out)
        for line in fh:
            jira_line = ""
            if column_headers:
                jira_line += "|"
            line_split = line.strip().split("\t")
            jira_line += "|" + "|".join(line_split) + "|"
            print(jira_line, file=out)
        out.close()

def main():
    parser = argparse.ArgumentParser(description="Convert a tsv file to Jira table format.")
    parser.add_argument("files", 
                        type=str, nargs="*",
                        help="Input tsv file to be converted.")
    parser.add_argument("-o", "--output",
                        type=str, default=None,
                        help="Output file (leave blank to print to stdout) [stdout].")
    parser.add_argument("-c", "--columnheaders",
                        default=False, action="store_true",
                        help="Create column headers [false].")
    parser.add_argument("-s", "--stdin",
                        action="store_true",
                        default=False,
                        help="Get input from stdin")
    parser.add_argument("-n", "--noheader",
                        action="store_true",
                        default=False,
                        help="Do not print a header row")
    args = parser.parse_args()

    if args.stdin:
        if len(args.files) > 1:
            print("to_jira: error: only one file can be taken if getting input from stdin")
            sys.exit(1)
        args.files.append("/dev/stdin")

    for f in args.files:
        convert_file(f, args.output,
            args.columnheaders, args.noheader)

if __name__ == "__main__":
    main()