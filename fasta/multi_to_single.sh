#!/bin/bash
# Convert a multi-line fasta file to single-line

if [[ $# -lt 1 ]]; then
	echo "Usage: $(basename $0) <multiline fasta>"
	exit 1
fi

if [[ $1 =~ .gz$ ]]; then
	method="gunzip -c"
else
	method="cat"
fi

${method} ${1} | awk '/^>/{if(N>0) printf("\n"); ++N; printf("%s\n",$0);next;} {printf("%s",$0);}END{printf("\n");}'
