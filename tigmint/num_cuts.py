#!/usr/bin/env python3
"""
Counts the number of cuts made to contigs in a bed file.
Usage: python3 /projects/btl/scratch/janli/scripts/num_cuts.py <bedfile>
Author: Janet Li @janetxinli
"""

import sys

class LastCutContig:
	def __init__(self):
		self.cuts = 0	
	def addCut(self):
		self.cuts += 1

def count_breaktigs(fh):
	"""Given the file handle of a bed file, count the breaktigs in the file."""
	total_cuts = 0
	cur_breaktig = LastCutContig()
	for line in fh:
		cur_contig = line.strip().split(sep="\t")
		contig_name = cur_contig[3].split("-")
		if len(contig_name) > 1:
			if contig_name[1] == "1":  # First breaktig of a contig
				total_cuts += cur_breaktig.cuts
				cur_breaktig = LastCutContig()  # Start a new breaktig
			else:
				cur_breaktig.addCut()

	# Add last contig cuts to total
	total_cuts += cur_breaktig.cuts
	return total_cuts

def main():
	if len(sys.argv) != 2 or not sys.argv[1].endswith(".bed"):
		print("num_cuts: error: must provide bed file as an argument.")
		print("usage: %s <breaktigs bed file>" % sys.argv[0])
		sys.exit(1)
	
	bed = sys.argv[1]
	with open(bed, "r") as contigs:
		cuts = count_breaktigs(contigs)	
	print(cuts)

if __name__ == "__main__":
	main()
