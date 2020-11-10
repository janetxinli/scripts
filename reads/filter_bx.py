#!/usr/bin/env python
"""
Filter, partition or subsample barcoded reads by a given read multiplicity range. Requires 
a longranger basic-processed reads file and the unprocessed reads from supernova mkfastq. 
Optionally subsamples reads to a certain coverage.
"""

import sys
import gzip
import os
import math
import time
import re
import argparse

def get_bx(read_file, min_reads, max_reads):
    """
    Opens read_file and gets the headers of reads with barcodes that have
    at least min_reads and at most max_reads reads. Returns a tuple of dict, num headers
    where dict = bx -> set(headers)
    """
    prev_bx = None
    all_headers = {}
    cur_headers = set()
    cur_reads = 0
    header_count = 0
    with open(read_file) as fh:
        if not fh.isatty():
            add_header = True
            for i, line in enumerate(fh):
                if i % 4 == 0:  # Header line
                    line_split = line.split(" ")
                    if len(line_split) > 1:  # Read has a bx
                        bx = line_split[1].strip()
                        header = line_split[0]
                        if bx != prev_bx:  # Reached the next bx
                            if cur_reads >= min_reads and cur_reads <= max_reads:
                                all_headers[prev_bx] = cur_headers
                                header_count += len(cur_headers)
                            add_header = True
                            cur_headers = {header}
                            cur_reads = 1
                            prev_bx = bx
                        else:
                            cur_reads += 1
                            if add_header:
                                if cur_reads > max_reads:
                                    add_header = False
                                    cur_headers = set()
                                else:
                                    cur_headers.add(header)  # Store valid header
            # Add the last headers if valid
            if cur_reads >= min_reads and cur_reads <= max_reads:
                all_headers[prev_bx] = cur_headers
                header_count += len(cur_headers)
            return all_headers, header_count
        else:
            print("filter_bx.py: error: barcoded reads must be piped from stdin or "
                "provided with the -b parameter", file=sys.stderr, flush=True)
            sys.exit(1)


def get_reads(read_file, valid_headers):
    """Given a list of valid headers, return a dictionary of header -> other read lines."""
    file_content = {}
    with gzip.open(read_file, "rt") as reads:
        add_lines = False
        for i, line in enumerate(reads):
            if i % 4 == 0:
                cur_header = line.split(" ")[0]
                if cur_header in valid_headers:
                    add_lines = True
                    file_content[cur_header] = [line.strip()]
                else:
                    add_lines = False
            else:
                if add_lines == True:
                    file_content[cur_header].append(line.strip())
    return file_content


def filter_reads(read_1, read_2, new_1, new_2, valid_headers):
    """
    Filters reads in paired-end read files that have a valid multiplicity
    (present in valid_headers).
    """
    files = {read_1: new_1, read_2: new_2}
    for read_file in files:
        new_file = gzip.open(files[read_file], "wt")
        with gzip.open(read_file, "rt") as fh:
            for i, line in enumerate(fh):
                if i % 4 == 0:
                    if line.strip().split(" ")[0] in valid_headers:
                        add_read = True
                        print(line.strip(), file=new_file, flush=True)
                    else:
                        add_read = False
                else:
                    if add_read:
                        print(line.strip(), file=new_file, flush=True)
        new_file.close()


def partition_reads(r1_content, r2_content, new_r1_files, new_r2_files, valid_bx, num_reads):
    """
    Partition a given read file into new files with a given number of reads,
    filtering for reads present in valid_headers (with a certain bx multiplicity).
    """
    f = 0
    cur_new_r1 = gzip.open(new_r1_files[f], "wt")
    cur_new_r2 = gzip.open(new_r2_files[f], "wt")
    cur_reads = 0
    for bx in valid_bx:
        if cur_reads >= num_reads:
                cur_new_r1.close()
                cur_new_r2.close()
                f += 1
                try:
                    cur_new_r1 = gzip.open(new_r1_files[f], "wt")
                    cur_new_r2 = gzip.open(new_r2_files[f], "wt")
                except IndexError:
                    return
                cur_reads = 0
        for h in valid_bx[bx]:
            print("\n".join(r1_content[h]), file=cur_new_r1)
            print("\n".join(r2_content[h]), file=cur_new_r2)
            cur_reads += 1
    # Close final files
    cur_new_r1.close()
    cur_new_r2.close()


def get_range(multiplicity_str):
    """Read multiplicity range from string format min-max (inclusive)."""
    mult_split = [int(i) for i in multiplicity_str.split("-")]
    if len(mult_split) != 2:
        print("filter_bx.py: barcode multiplicity range must be provided in the format min-max",
            file=sys.stderr, flush=True)
        sys.exit(1)
    min_mult, max_mult = mult_split[0], mult_split[1]
    if min_mult >= max_mult:
        print("filter_bx.py: max of multiplicity range provided is not larger than the min",
            file=sys.stderr, flush=True)
        sys.exit(1)
    return mult_split[0], mult_split[1]


def ensure_writable(file_name):
    """Ensure file_name is writable."""
    if os.access(file_name, os.W_OK):
        return
    if os.access(file_name, os.F_OK):
        print("filter_bx.py: error: file '{0}' exists and is not writable".format(file_name),
            file=sys.stderr, flush=True)
        sys.exit(1)


def get_int(size_string):
    """
    Convert genome size from string to integer. Allows string to be formatted as an 
    explicit integer or in scientific notation.
    """
    try:
        return int(float(size_string))
    except ValueError:
        print("filter_bx.py: error: genome size must be given as an integer or in scientific notation (e.g. 3e9)",
            file=sys.stderr, flush=True)
        sys.exit(1)


def get_new_file_names(read_prefix, read_suffix, min_mult, max_mult, num_files):
    """Return a list of new read file names for partitioned reads."""
    files = []
    for part in range(1, num_files + 1):
        if min_mult == -1 and max_mult == 9e12:  # Keeping all reads
            filter_name = "_all"
        else:
            filter_name = "{0}-{1}".format(min_mult, max_mult)
        files.append(read_prefix + "_filterbx{0}_partition{1}".format(filter_name, part) + read_suffix)
    return files


def get_read_info(r1, r2):
    """Return a tuple of read prefix and suffixes for r1 and r2 file names."""
    regex = "^(.+)(_S\d+_L00\d+_R[12]_001.fastq.gz)"
    search1 = re.search(regex, r1)
    search2 = re.search(regex, r2)
    read_prefix = search1.group(1)
    if read_prefix != search2.group(1):
        raise ValueError("Read files have different prefixes")
    else:
        read1_suffix = search1.group(2)
        read2_suffix = search2.group(2)
        return read_prefix, read1_suffix, read2_suffix


def run_filter(r1, r2, min_mult, max_mult, read_prefix, r1_suffix, r2_suffix, valid_headers):
    """Filter reads by barcode multiplicity."""
    if min_mult == -1 and max_mult == 9e12:
        filter_name = "_all"
    else:
        filter_name = "{0}-{1}".format(min_mult, max_mult)
    new_r1 = read_prefix + "_filterbx{0}".format(filter_name) + r1_suffix
    new_r2 = read_prefix + "_filterbx{0}".format(filter_name) + r2_suffix
    print("filtering reads...", file=sys.stderr, flush=True)
    filter_reads(r1, r2, new_r1, new_r2, valid_headers)
    print("reads filtered into new files:", new_r1, new_r2, sep="\n", flush=True)


def run_partition(r1, r2, new_r1_files, new_r2_files, valid_bx, valid_headers, reads_per_file):
    """Partition read files to a certain number of reads."""
    print("reading in r1...", file=sys.stderr, flush=True)
    r1_content = get_reads(r1, valid_headers)
    print("reading in r2...", file=sys.stderr, flush=True)
    r2_content = get_reads(r2, valid_headers)
    print("partitioning reads...", file=sys.stderr, flush=True)
    partition_reads(r1_content, r2_content, new_r1_files, new_r2_files, valid_bx, reads_per_file)

def get_args():
    """Parse the command line arguments."""
    parser = argparse.ArgumentParser(description="Filter or partition reads in a fastq file by \
                            barcode multiplicity or read coverage")
    parser.add_argument("mode",
                        choices=["filter", "partition", "subsample"],
                        default="filter")
    parser.add_argument("r1",
                        type=str,
                        help="Gzipped file containing read 1")
    parser.add_argument("r2",
                        type=str,
                        help="Gzipped file containing read 2")
    parser.add_argument("-m", "--multiplicity",
                        type=str,
                        required=True,
                        help="Bx multiplicity range for reads to be included in output (e.g. 10-150)")
    parser.add_argument("-c", "--coverage",
                        type=int,
                        default=None,
                        help="Desired read coverage for partitioned reads")
    parser.add_argument("-g", "--genome",
                        type=str,
                        default=None,
                        help="Haploid genome size; required to calculate read number in partition mode")
    parser.add_argument("-l", "--read_length",
                        type=int,
                        default=150,
                        help="Length of reads")
    parser.add_argument("-n", "--num_reads",
                        type=str,
                        default=None,
                        help="Number of reads to include in each partition")
    parser.add_argument("-b", "--barcoded_reads",
                        type=str,
                        default="-",
                        help="File containing barcoded reads [stdin]")
    return parser.parse_args()


def main():
    """Filter input reads."""
    start = time.asctime()
    print("filter_bx.py: started at: {0}".format(start), file=sys.stderr,
        flush=True)
    args = get_args()

    if args.coverage and args.num_reads:
        print("filter_bx.py: error: coverage and number of reads both provided. please provide one",
            file=sys.stderr, flush=True)
        sys.exit(1)

    if args.coverage:
        if not args.genome:
            print("filter_bx.py: error: haploid genome size must be required to partition reads by coverage", 
                file=sys.stderr, flush=True)
            sys.exit(1)
        args.genome = get_int(args.genome)

    try:
        read_prefix, read1_suffix, read2_suffix = get_read_info(args.r1, args.r2)
    except ValueError:
        sys.exit(1)

    if args.multiplicity == "all":
        min_mult = -1
        max_mult = 9e32  # Cheat
    else:
        min_mult, max_mult = get_range(args.multiplicity)

    print("finding valid headers...", file=sys.stderr, flush=True)
    if args.barcoded_reads == "-":
        args.barcoded_reads = "/dev/stdin"
    valid_bx, num_headers = get_bx(args.barcoded_reads, min_mult, max_mult)
    valid_headers = set(h for bx in valid_bx.values() for h in bx)

    if args.mode == "filter":  # Filter reads by bx multiplicity
        run_filter(args.r1, args.r2, min_mult, max_mult, read_prefix, read1_suffix, read2_suffix, valid_headers)

    elif args.mode == "partition":  # Filter reads by bx multiplicity and partition
        total_valid_reads = num_headers * 2
        if args.coverage:  # Partition by coverage
            total_cov = (total_valid_reads * args.read_length) // args.genome
            new_files_per_read = math.ceil((total_cov / args.coverage))
            reads_per_file = ((args.coverage * args.genome) // args.read_length) // 2
        elif args.num_reads:  # Partition by number
            args.num_reads = get_int(args.num_reads)
            reads_per_file = args.num_reads // 2
            new_files_per_read = math.ceil((total_valid_reads / 2) / reads_per_file)

        new_r1_files = get_new_file_names(read_prefix, read1_suffix, min_mult, max_mult, new_files_per_read)
        new_r2_files = get_new_file_names(read_prefix, read2_suffix, min_mult, max_mult, new_files_per_read)
        print("{0} valid reads".format(total_valid_reads), file=sys.stderr, flush=True)
        print("partitioning each read direction into {0} files with approx. {1} reads each".format(new_files_per_read, reads_per_file),
        file=sys.stderr, flush=True)
        run_partition(args.r1, args.r2, new_r1_files, new_r2_files, valid_bx, valid_headers, reads_per_file)
    
    elif args.mode == "subsample":  # Subsample reads
        if args.coverage:
            reads_per_file = ((args.coverage * args.genome) // args.read_length) // 2
        elif args.num_reads:
            reads_per_file = args.num_reads // 2
        
        new_r1_files = get_new_file_names(read_prefix, read1_suffix, min_mult, max_mult, 1)
        new_r2_files = get_new_file_names(read_prefix, read2_suffix, min_mult, max_mult, 1)
        print("{0} valid reads".format(total_valid_reads), file=sys.stderr, flush=True)
        print("subsampling {0} reads from each direction".format(reads_per_file),
            file=sys.stderr, flush=True)
        run_partition(args.r1, args.r2, new_r1_files, new_r2_files, valid_bx, valid_headers, reads_per_file)

    print("DONE!", file=sys.stderr, flush=True)
    end = time.asctime()
    print("filter_bx.py: ended at: {0}\nelapsed time: {1}".format(end, time.process_time()), file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
