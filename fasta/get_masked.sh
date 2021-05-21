#!/bin/bash

if [[ $# -lt 1 ]]; then
	echo "Usage: $(basename $0) <genome fasta files>"
	exit 1
fi

fasta=$@

for f in ${fasta[@]}; do
	hard=$(grep -v "^>" ${f} | tr -d "actgATCGMRWSYKVHDBumrwsykvhdb" | tr -d "\n" | wc -c)
	soft=$(grep -v "^>" ${f} | tr -d "ATCGTMRWSYKVHDBNumrwsykvhdbn" | tr -d "\n" | wc -c)
	ambig=$(grep -v "^>" ${f} | tr -d "ATCGNactgn" | tr -d "\n" | wc -c)

	printf "${f}\t${hard}\t${soft}\t${ambig}\n"
done
