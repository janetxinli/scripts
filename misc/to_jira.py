#!/usr/bin/env python3

import os
import sys
import argparse

def format_number(number):
    """Formats a number, rounding it to 2 decimal places and adding commas."""
    try:
        if number.isnumeric():
            formatted = "{:,}".format(int(number))
        else:
            formatted = "{:,.2f}".format(float(number))
    except ValueError:
        formatted = number
    finally:
        return formatted

def convert_file(in_fh, outfile, column_headers, no_header, no_rounding, delim):
    """Takes as input a tsv file handle and formats the input for Jira."""
    if outfile == None:
        out = sys.stdout
    else:
        out = open(outfile, "w+")
    header_split = [i for i in in_fh.readline().strip().split(delim)]
    if not no_rounding:
        header_split = [format_number(i) for i in header_split]
    if no_header and not column_headers:
        jira_header = "|" + "|".join(header_split) + "|"
    elif no_header and column_headers:
        jira_header = "||" + "|".join(header_split) + "|"
    else:
        jira_header = "||" + "||".join(header_split) + "||"
    print(jira_header, file=out)
    for line in in_fh:
        jira_line = ""
        if column_headers:
            jira_line += "|"
        line_split = [i for i in line.strip().split(delim)]
        if not no_rounding:
            line_split = [format_number(i) for i in line_split]
        jira_line += "|" + "|".join(line_split) + "|"
        print(jira_line, file=out)
    out.close()

def main():
    parser = argparse.ArgumentParser(description="Convert a file to Jira table format")
    parser.add_argument("file", 
                        type=str,
                        nargs="?",
                        default="/dev/stdin",
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
    parser.add_argument("-r", "--no_rounding",
                        action="store_true",
                        default=False,
                        help="Do not round floating point numbers")
    parser.add_argument("-d", "--delim",
                        default="\t",
                        help="Delimiter used in input file [tab]")
    args = parser.parse_args()

    with open(args.file) as fh:
        if not fh.isatty():
            convert_file(fh, args.output,
                args.columnheaders, args.noheader,
                args.no_rounding, args.delim)
        else:
            print("to_jira.py: error: missing input, pass a file as an argument or from stdin")
            sys.exit(1)

if __name__ == "__main__":
    main()
