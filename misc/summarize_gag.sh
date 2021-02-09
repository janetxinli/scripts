#!/bin/bash

files=$@

af=${files[0]}

printf "assembly"
for f in ${files[@]}; do
	printf "\t${f}"
done
printf "\n"

for line in {3..39}; do
	row=$(sed "${line}q;d" ${af} | awk 'BEGIN {FS="  "} {print $1}')
	printf "${row}"
	for f in ${files[@]}; do
		val=$(sed "${line}q;d" ${f} | awk '{print $NF}')
		printf "\t${val}"
	done
	printf "\n"
done
