#!/usr/bin/env python

from BIR import *
import json
import sys
import os

filePath = sys.argv[1]
directory = sys.argv[2]
data = parseFastaFile(filePath)

seq = ''
for key in data:
    seq = data[key]

monomerCount = open(directory + '/monomerCount.JSON', 'a')
monomerCount.write(json.dumps(countMonomers(seq)))
monomerCount.close()

rComplement = open(directory + '/reverseComplement.txt', 'a')
rComplement.write(reverseComplement(seq))
rComplement.close()

# RNAFullSeq = open('D:/Python27/my_python_programs/DNAEval/RNAFullSeq.txt', 'a')
# RNAFullSeq.write()
# RNAFullSeq.close()

# RNACodingSeq = open('D:/Python27/my_python_programs/DNAEval/RNACodingSeq.txt', 'a')
# RNACodingSeq.write()
# RNACodingSeq.close()

# proteinFullSeq = open('D:/Python27/my_python_programs/DNAEval/proteinFullSeq.txt', 'a')
# proteinFullSeq.write()
# proteinFullSeq.close()

# proteinCodingSeq = open('D:/Python27/my_python_programs/DNAEval/proteinCodingSeq.txt', 'a')
# proteinCodingSeq.write()
# proteinCodingSeq.close()
