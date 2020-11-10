#!/usr/bin/env python
"""Count the number of cuts made to contigs in a tigmint breaktig bed file."""

import sys

class Contig:
	"""Contig object that keeps track of number of cuts made."""
	def __init__(self):
		self.cuts = 0	
	def addCut(self):
		self.cuts += 1


def count_breaktigs(fh):
	"""Given the file handle of a bed file, count the breaktigs in the file."""
	total_cuts = 0
	cur_breaktig = Contig()
	for line in fh:
		cur_contig = line.strip().split(sep="\t")
		contig_name = cur_contig[3].split("-")
		if len(contig_name) > 1:
			if contig_name[1] == "1":  # First breaktig of a contig
				total_cuts += cur_breaktig.cuts
				cur_breaktig = Contig()  # Start a new breaktig
			else:
				cur_breaktig.addCut()

	# Add last contig cuts to total
	total_cuts += cur_breaktig.cuts
	return total_cuts


def main():
	if len(sys.argv) != 2 or not sys.argv[1].endswith(".bed"):
		print("num_cuts.py: error: must provide bed file as an argument", file=sys.stderr, flush=True)
		print("usage: %s <breaktigs bed file>" % sys.argv[0], file=sys.stderr, flush=True)
		sys.exit(1)
	
	bed = sys.argv[1]
	with open(bed, "r") as contigs:
		cuts = count_breaktigs(contigs)	
	print(cuts, file=sys.stdout)


if __name__ == "__main__":
	main()
