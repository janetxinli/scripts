"""Print a histogram with a specified bin width given a list of values."""

import sys
import numpy as np

def print_histogram(values, binwidth, outfile=None):
    """
    Given a list of values and desired bin width, prints a histogram of values
    in tsv format. By default prints histogram to stdout.
    """
    max_val = np.max(values)
    # Calculate bin width and max
    upper_lim = max_val + (binwidth - (max_val % binwidth))
    bins = upper_lim // binwidth
    counts, bins = np.histogram(values, bins, range=(0, upper_lim))
    if outfile is None:
        of = sys.stdout
    else:
        of = open(outfile, "w+")
    for i, bin in enumerate(bins[1:]):
        print(int(bin), counts[i], sep="\t", file=of)
    of.close()