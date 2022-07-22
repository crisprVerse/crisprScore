#!/bin/bash
#$ -cwd

while read -r line; do
	g=$(echo $line | cut -d "," -f 1 )
	t=$(echo $line | cut -d "," -f 2 )
	RNAhybrid -c -s 3utr_human ${g} ${t}
done < $1