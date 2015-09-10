#!/usr/bin/env python

from BIR import *
import json
import sys
import os

filePath = sys.argv[1]
directory = sys.argv[2]
data = parseFastaFile(filePath)

#assumes there is only one FASTA formatted sequence when unpacking
seq = ''
for key in data:
    seq = data[key]

print "Evaluating forward read..."

monomerCount = open(directory + '/monomerCount.JSON', 'a')
monomerCount.write(json.dumps(countMonomers(seq)))
monomerCount.close()

rComplement = open(directory + '/reverseComplement.txt', 'a')
revComp = reverseComplement(seq)
rComplement.write(revComp)
rComplement.close()

RNAFullSeq = open(directory + '/RNAFullSeq.txt', 'a')
fullRna = dna2Rna(seq)
RNAFullSeq.write(fullRna)
RNAFullSeq.close()

RNACodingSeq = open(directory + '/RNACodingSeq.txt', 'a')
codingRna = rnaCodingSeq(fullRna)
RNACodingSeq.write(codingRna)
RNACodingSeq.close()

proteinCodingSeq = open(directory + '/proteinCodingSeq.txt', 'a')
proteinCodingSeq.write(rna2Protein(codingRna))
proteinCodingSeq.close()

print "Evaluating reverse complement read..."

revRNAFullSeq = open(directory + '/revRNAFullSeq.txt', 'a')
revfullRna = dna2Rna(revComp)
revRNAFullSeq.write(revfullRna)
revRNAFullSeq.close()

revRNACodingSeq = open(directory + '/revRNACodingSeq.txt', 'a')
revcodingRna = rnaCodingSeq(revfullRna)
revRNACodingSeq.write(revcodingRna)
revRNACodingSeq.close()

revproteinCodingSeq = open(directory + '/revproteinCodingSeq.txt', 'a')
revproteinCodingSeq.write(rna2Protein(revcodingRna))
revproteinCodingSeq.close()