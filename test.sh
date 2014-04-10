#!/bin/bash

if [ $1 == "local" ]; then
	type=local
	matrix=pam120.txt
	gap_open=5
	gap_extend=1
	sequences=maguk-sequences.fasta
	size=1
elif [ $1 == "global" ]; then
	type=global
	matrix=blosum62.txt
	gap_open=5
	gap_extend=2
	sequences=PDZ-sequences.fasta
	size=5
fi

echo -e "CAREFUL : THIS program may take a while to finish as it tests all sequences against each other\n"

echo "type: " $type
echo "matrix: " $matrix
echo "opening gap penalty:" $gap_open
echo "extending gap penalty:" $gap_extend
echo -e "sequences file:" $sequences "\n"

for i in $(seq 0 $size); do
	for j in $(seq 0 $size); do
		if [ $i -ne $j ]; then
			echo -e "sequence #"$i " against sequence #"$j "\n"
			echo -e "$type\n$matrix\n$gap_open\n$gap_extend\n$sequences\n$i\n$sequences\n$j\n" | python Peche_A_L_Align.py > /dev/null
		fi
	done
done