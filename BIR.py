#This repository includes the following functions:
#parseFastaFile(filename)
#countMonomers(seq)
#greatestContent(sequences, monomers)
#hammingDistance(seq1, seq2)
#indexSubstring(substring, string)
#reverseComplement(dna)
#consensusSequence(seq)
#dna2Rna(dna)
#indexStartCodon(mRna)
#indexStopCodon(mRna)
#rnaCodingSeq(mRna)
#rna2Codons(rna)
#codon2Aa(codon)
#rna2Protein(rna)

#---------------------------------------------------
#This function takes a DNA sequence and returns the reverse complement.
#The reverse complement is made by first taking the reverse, then
#using the complementary base for each ('G' for 'C', 'A' for 'T').
def reverseComplement(dna):
    revComp = ''

    revSeq = dna[::-1]
    for base in revSeq:
        if base == 'A':
            revComp += 'T'
        if base == 'T':
            revComp += 'A'
        if base == 'G':
            revComp += 'C'
        if base == 'C':
            revComp += 'G'
            
    return revComp
#---------------------------------------------------

#---------------------------------------------------
#This function takes a substring and a string, and returns an array
#of all positions of the substring in the string using 0-based count.
#If the substring is bigger than the string or no positions are found,
#it returns an empty array.
def indexSubstring(substring, string):
	array = []
	pos = 0

	for pos in range(len(string) - len(substring) + 1):
		if substring == string[pos : pos + len(substring)]:
			array.append(pos)
			
	return array
#---------------------------------------------------


#---------------------------------------------------
#This function takes two strings and returns the Hamming distance.
#The Hamming distance is defined as the number of mismatches
#between sequences.
def hammingDistance(seq1, seq2):
    dH = 0
    
    if not(len(seq1) == len(seq2)):
        print 'lengths not equal, errors may result    --hammingDist--'
    for pos in range(len(seq1)):
        if not(seq1[pos] == seq2[pos]):
            dH += 1
            
    return dH
#---------------------------------------------------


#---------------------------------------------------
#This function takes a file (in FASTA form) and returns a dictionary
#of tags and their corresponding sequences. filename is the directory.
def parseFastaFile(filename):
    dictionary = {}
        
    with open(filename, 'r') as r:
        for line in r:
            if line[0] == '>':
                tag = line[1:-1]
                dictionary[tag] = ''
                seq = ''
            else:
                dictionary[tag] += line.rstrip('\n')
                
    return dictionary
#---------------------------------------------------


#---------------------------------------------------
#This function takes a list of sequences and a list of monomers, and
#returns the sequence with the greatest %content of those monomers.
def greatestContent(sequences, monomers):
    contentList = []
    greatest = 0
    
    for seq in sequences:
        total = 0
        occurences = countDnaBases(seq)
        for monomer in monomers:
            total += occurences[monomer]
        contentList.append(float(total)/len(seq))

    for pos in range(len(contentList)):
        if contentList[pos] > greatest:
            greatest = contentList[pos]
            winner = sequences[pos]
            
    return winner
#---------------------------------------------------


#---------------------------------------------------
#This function takes a sequence and returns a dictionary of the number
#of occurences of each monomer. The dictionary uses each monomer as a key
#and their occurences as the values.
def countMonomers(seq):
    results = {}
    
    for pos in range(len(seq)):
        monomer = seq[pos]
        if (monomer in results):
            results[monomer] += 1
        else: results[monomer] = 1
    
    return results
#---------------------------------------------------


#---------------------------------------------------
#This function
def consensusSequence(seqs):
    consensusSeq = buildConcensus(buildProfile(seqs))
    
    return consensusSeq
            
#---------------------------------------------------


#---------------------------------------------------
#This function takes a list of sequences and returns the profile
#The profile is a dictionary of each type of character in the sequences,
#which each have value of an array. This array holds the number of occurences
#of that character at each position in the sequences.
def buildProfile(seqs):
    profile = {}
    
    for pos in range(len(seqs[0])):
        for key in profile:
            profile[key].append(0)
        for seq in seqs:
            monomer = seq[pos]
            if not(monomer in profile):
                profile[monomer] = []
                for i in range(pos+1):
                    profile[monomer].append(0)
            profile[monomer][pos] += 1
                
    return profile
#---------------------------------------------------

                
#---------------------------------------------------
#fix me
def buildConsensus(profile):
    consensus = ['']
    pos = 0
    
    winner = 0
    for key in profile:
        if profile[key][pos] > winner:
            winner = profile[key][pos]
            consensus[pos] = key
    pos += 1
    while (pos < len(profile[key])):
        winner = 0
        consensus.append('')
        for key in profile:
            if profile[key][pos] > winner:
                winner = profile[key][pos]
                consensus[pos] = key
        pos += 1
    consensusSeq = ''.join(consensus)

    return consensusSeq
#---------------------------------------------------


#---------------------------------------------------
#This function takes a DNA string and returns the corresponding RNA string.
def dna2Rna(dna):
    rna = ''
    for pos in range(len(dna)):
        if dna[pos] == 'T':
            rna += 'U'
        else: rna += dna[pos]
            
    return rna
#---------------------------------------------------


#---------------------------------------------------
#This function takes a mRNA sequence and returns the position of the first
#start codon ('AUG'). The position is where the start codon starts. If no start
#is found, None is returned.
def indexStartCodon(mRna):
    for pos in range(len(mRna)-2):
        if mRna[pos:pos+3] == 'AUG':
            start = pos
            return start
        
    print "no start codon found   --posStartCodon--"
    return None
#---------------------------------------------------


#---------------------------------------------------
#This function takes a mRNA sequence and returns the position of the first
#stop codon ('UAA', 'UAG', or 'UGA') that is in the initial reading frame.
#If no stop is found, None is returned.
def indexStopCodon(mRna):

    for codonPos in range(len(mRna)//3):
        if mRna[codonPos*3] == 'U':
            if mRna[codonPos*3+1] == 'A':
                if mRna[codonPos*3+2] == 'A' or mRna[codonPos*3+2] == 'G':
                    end = codonPos*3 + 3
                    return end
            if mRna[codonPos*3+1] == 'G':
                if mRna[codonPos*3+2] == 'A':
                    end = codonPos*3
                    return end
                
    return None
#---------------------------------------------------


#---------------------------------------------------
#This function takes an mRNA sequence and returns the slice of the
#sequence starting at the first found 'AUG', and ending at the next
#found 'UAA', 'UAG', 'UGA' or end. The function returns None if there
#is no start.
def rnaCodingSeq(mRna):	
    start = indexStartCodon(mRna)
    end = indexStopCodon(mRna[start:])+4
    if start == None:
        print "no start codon was found   --rnaCodingSeq--"
        return None
    if end == None:
        end = len(mRna)
        print "no stop codon was found    --rnaCodingSeq--"
        
    return mRna[start:end]
#---------------------------------------------------


#---------------------------------------------------
#This function takes a RNA sequence and returns a list of codons. Codons are
#three character strings. If the final codon is incomplete, it will not be
#included.
def rna2Codons(rna):
    codons = []

    if not (len(rna)%3 == 0):
        print "final codon will be incomplete, errors may result    --rna2Codons--"
    for j in range(len(rna)//3):    #for every 3 bases
        codon = ''
        for i in range(3):      #build 3 base codon
            if (3*j+i) < len(rna):
                codon += rna[3*j+i]
        codons.append(codon)
        
    return codons
#---------------------------------------------------

    
#---------------------------------------------------
#This function takes a codon string and returns the single letter code for the
#corresponding amino acid. Stop codons return 'stop'.
def codon2Aa(codon):
    dictionary = {'UUU' : 'F',
            'UUC' : 'F',
            'UUA' : 'L',
            'UUG' : 'L',
            'UCU' : 'S',
            'UCC' : 'S',
            'UCA' : 'S',
            'UCG' : 'S',
            'UAU' : 'Y',
            'UAC' : 'Y',
            'UAA' : 'stop',
            'UAG' : 'stop',
            'UGU' : 'C',
            'UGC' : 'C',
            'UGA' : 'stop',
            'UGG' : 'W',

            'CUU' : 'C',
            'CUC' : 'C',
            'CUA' : 'C',
            'CUG' : 'C',
            'CCU' : 'P',
            'CCC' : 'P',
            'CCA' : 'P',
            'CCG' : 'P',
            'CAU' : 'H',
            'CAC' : 'H',
            'CAA' : 'Q',
            'CAG' : 'Q',
            'CGU' : 'R',
            'CGC' : 'R',
            'CGA' : 'R',
            'CGG' : 'R',

            'AUU' : 'I',
            'AUC' : 'I',
            'AUA' : 'I',
            'AUG' : 'M',
            'ACU' : 'T',
            'ACC' : 'T',
            'ACA' : 'T',
            'ACG' : 'T',
            'AAU' : 'N',
            'AAC' : 'N',
            'AAA' : 'K',
            'AAG' : 'K',
            'AGU' : 'S',
            'AGC' : 'S',
            'AGA' : 'R',
            'AGG' : 'R',

            'GUU' : 'V',
            'GUC' : 'V',
            'GUA' : 'V',
            'GUG' : 'V',
            'GCU' : 'A',
            'GCC' : 'A',
            'GCA' : 'A',
            'GCG' : 'A',
            'GAU' : 'D',
            'GAC' : 'D',
            'GAA' : 'E',
            'GAG' : 'E',
            'GGU' : 'G',
            'GGC' : 'G',
            'GGA' : 'G',
            'GGG' : 'G'}
    aminoAcid = dictionary[codon]
    
    return aminoAcid
#-----------------------------------------------------


#-----------------------------------------------------
#This function takes a RNA sequence and returns the corresponding protein
#sequence. The protein ends at the first stop codon or end of the RNA.
def rna2Protein(rna):
    protein = ''

    codons = rna2Codons(rna)
    for codon in codons:
        if codon2Aa(codon) == 'stop':
            return protein
        else:
            protein += codon2Aa(codon)
            
    return protein
#------------------------------------------------------
##------------------------------------------------------------------------------------
