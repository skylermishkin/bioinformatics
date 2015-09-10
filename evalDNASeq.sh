#!/bin/bash

#set up variables for the full file path and the directory for the results (same as the file)
if (($# == 0)); then
    echo "What is the full path to the FASTA DNA sequence?"
    read path
else
    path="$1"
fi
directory=$(dirname "${path}")

#touch result files
touch "$directory"/monomerCount.JSON
touch "$directory"/reverseComplement.txt
touch "$directory"/RNAFullSeq.txt
touch "$directory"/RNACodingSeq.txt
touch "$directory"/ProteinCodingSeq.txt
touch "$directory"/revRNAFullSeq.txt
touch "$directory"/revRNACodingSeq.txt
touch "$directory"/revProteinCodingSeq.txt

#run the python script, passing it the file path and directory
python evalDNASeq.py "$path" "$directory"

#run tests, check for errors in BIR
#todo

echo evalDNASeq has finished