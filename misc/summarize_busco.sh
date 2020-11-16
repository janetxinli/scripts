#!/bin/bash
# Summarize a busco short summary file

if [[ $# -lt 1 ]]; then
	echo "Usage: $(basename $0) <busco short summary files>"
	exit 1
fi

b=$(basename $0)
summaries=$@; shift

printf "complete\tcomplete_single\tcomplete_dup\tfragmented\tmissing\ttotal\tcomplete_perc\tname\n"

for s in ${summaries[@]}; do
	f=$(basename ${s})
	if [[ "${f}" != "${b}" ]]; then
		for v in $(tail -n6 ${s} | awk '{print $1}'); do
			printf "$v\t"
		done
		printf "$(grep "C:" ${s} | awk -F : '{print $2}' | awk -F % '{print $1}')\t"
		printf "${f}\n"
	fi
done
