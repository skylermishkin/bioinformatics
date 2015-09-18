#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Thu Sep 17 13:36:39 2015

@author: Skyler Mishkin

This file contains functions related to exact and approximate sequence matching
algorithms. The outputs of all 'exact_matching_*Method*()' and 
'approximate_matching_*Method*()' functions are lists of indices for the text 
where the pattern was identified. Throughout this file, 't' is used to represent
the reference sequence which will be searched against, and 'p' is used to represent 
the pattern sequence which will be searched for.

Specifically, this file contains matching algorithms using:
    t (reference sequence) preprocessing:
        Kmer Indexing  |  class/method: KmerIndex()
        Subsequence Indexing  |  class/method: SubseqIndex()
    p (pattern sequence) preprocessing:
        Boyer-Moore  |  class/method: BoyerMoore()
    
Note: when dealing with massive sequences (or a slow processor) it is recommended 
to do the preprocessing first (perhaps saving to a file), and to then run the 
matching algorithm on the preprocessed data (rather than combining them into one go). 
If you want to try it all in one go, there are 'oneshot_exact_matching_...' and 
'oneshot_approximate_matching_...' functions for the indexing methods.
    
This file derives from the second module of the Coursera course: Algorithms for DNA 
Sequencing, from Johns Hopkins University; some code comes directly from course 
resources. Other code was discussed theoretically. And other code is original. I 
have made modifications throughout this file, and make no guarantees or 
assurances of any form.
"""
import bisect

#The following block contains an exact matching algorithm using a brute force method
#==============================================================================    
def exact_matching_naive(s, l):
    ''' This function uses brute force and is not suitable for long sequences.
    Its purpose is more demonstrational. '''
    occurrences = []
    aligns = 0  #to demonstrate how many alignments are checked (l - s + 1)
    comps = 0  #to demonstrate how many character comparisons are made
    if (len(l) < len(s)):
        print "the second sequence was longer than the first, switch arguments   --exact_matching_naive--"
        return False
    for i in range(len(l) - len(s) + 1):  # loop over alignments
        aligns += 1    
        match = True
        for j in range(len(s)):  # loop over characters
            comps += 1
            if l[i+j] != s[j]:  # compare characters
                match = False
                break
        if match:
            occurrences.append(i)  # all chars matched; record
    return occurrences#, aligns, comps
#============================================================================== 



#The following block contains code using a Kmer Indexing method
#==============================================================================     
class KmerIndex:
    ''' Create an index of all length k portions of t and their indices 
    in t. Binary search is used to query the index. E.g., KmerIndex('ACGT', 2) 
    extracts ('AC', 0), ('CG', 1), and ('GT', 2). '''
    def __init__(self, t, k):
        self.k = k
        self.index = []
        for i in range(len(t) - k + 1):
            self.index.append((t[i:i+k], i))
        self.index.sort()
    
    def query_prefix(self, p):
        kmer = p[:self.k]
        i = bisect.bisect_left(self.index, (kmer, -1))
        hits = []
        while i < len(self.index):
            if self.index[i][0] != kmer:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits


def exact_matching_KmerIndex(index, p, t):
    ''' Confirms the index hits are full matches.  '''
    k = index.k
    offsets = []
    for i in index.query_prefix(p):
        if p[k:] == t[i+k:i+len(p)]:
            offsets.append(i)
    return offsets
    
    
def approximate_matching_KmerIndex(index, p, t, mm=2):
    ''' Requires p to be evenly divisible by mm+1 and for that value to be 
    greater or equal to than k. These requirements are imposed by the use of a 
    'pigeon hole' principle dividing p into segments to query against the index. '''
    if (len(p) / (mm+1)) % 1 != 0:
        print 'p must be evenly divisible by mm+1, execution aborted.    --approximate_matching_KmerIndex--'
        return -1
    elif index.k > int(round(len(p) / (mm+1))):
        print 'k must be smaller than p / (mm+1), execution aborted.     --approximate_matching_KmerIndex--'
        return -1
        
    all_matches = set()
    count = 0
    segment_length = int(round(len(p) / (mm+1)))
    for num in range(mm+1):
        start = num*segment_length
        end = min((num+1)*segment_length, len(p))
        hits = index.query_prefix(p[start:end])
        count += len(hits)
        for i in hits:
            mismatches = 0
            for j in range(len(p)):  #scan p for mismatches, would be nice to skip kmer part
                if p[j] != t[i+j-start]:
                    mismatches += 1
                    if mismatches > mm:
                        break
            if mismatches <= mm:
                all_matches.add(i - start)
    print 'Hits:', count
    return list(all_matches)
        
        
def oneshot_approximate_matching_KmerIndex(p, t, k=8, mm=2):
    ''' Requires p to be evenly divisible by mm+1 and for that value to be 
    greater or equal to than k.These requirements are imposed by the use of a 
    'pigeon hole' principle dividing p into segments to query against the index. '''
    if (len(p) / (mm+1)) % 1 != 0:
        print 'p must be evenly divisible by mm+1, execution aborted.    --oneshot_approximate_matching_KmerIndex--'
        return -1
    elif k > int(round(len(p) / (mm+1))):
        print 'k must be smaller than p / (mm+1), execution aborted.     --oneshot_approximate_matching_KmerIndex--'
        return -1
        
    kmer_index = KmerIndex(t, k)
    return approximate_matching_KmerIndex(kmer_index, p, t, mm)
#============================================================================== 
    
    
    
#The following block contains code using a Subsequence Indexing method
#============================================================================== 
class SubseqIndex(object):
    """ Holds a subsequence index for a text T """
    
    def __init__(self, t, k, ival):
        """ Create index from all subsequences consisting of k characters
            spaced ival positions apart.  E.g., SubseqIndex("ATAT", 2, 2)
            extracts ("AA", 0) and ("TT", 1). """
        self.k = k  # num characters per subsequence extracted
        self.ival = ival  # space between them; 1=adjacent, 2=every other, etc
        self.index = []
        self.span = 1 + ival * (k - 1)
        for i in range(len(t) - self.span + 1):  # for each subseq
            self.index.append((t[i:i+self.span:ival], i))  # add (subseq, offset)
        self.index.sort()  # alphabetize by subseq
    
    def query_prefix(self, p):
        """ Return index hits for first subseq of p """
        subseq = p[:self.span:self.ival]  # query with first subseq
        i = bisect.bisect_left(self.index, (subseq, -1))  # binary search
        hits = []
        while i < len(self.index):  # collect matching index entries
            if self.index[i][0] != subseq:
                break
            hits.append(self.index[i][1])
            i += 1
        return hits 


def exact_matching_SubseqIndex(index, p, t):
    ''' Confirms the index hits are full matches.  '''
    offsets = []
    for i in index.query_prefix(p):
        if p == t[i:i+len(p)]:
            offsets.append(i)
    return offsets
    
    
def approximate_matching_SubseqIndex(index, p, t, mm=2, ival=3):
    ''' Requires p to be evenly divisible by mm+1 and for that value to be 
    greater or equal to than k. These requirements are imposed by the use of a 
    'pigeon hole' principle dividing p into segments to query against the index. '''
    if (len(p) / (mm+1)) % 1 != 0:
        print 'p must be evenly divisible by mm+1, execution aborted.    --approximate_matching_SubseqIndex--'
        return -1
    elif index.k > int(round(len(p) / (mm+1))):
        print 'k must be smaller than p / (mm+1), execution aborted.     --approximate_matching_SubseqIndex--'
        return -1
        
    all_matches = set()
    count = 0
    for num in range(mm+1):
        seg = p[num:]
        hits = index.query_prefix(seg)
        count += len(hits)
        for i in hits:
            mismatches = 0
            for j in range(len(p)):  #scan p for mismatches, would be nice to skip ivals
                if p[j] != t[i+j-num]:
                    mismatches += 1
                    if mismatches > mm:
                        break
            if mismatches <= mm:
                all_matches.add(i - num)
    print 'Index hits:', count
    return list(all_matches)
    
    
def oneshot_approximate_matching_SubseqIndex(p, t, k=8, mm=2, ival=3):
    ''' Requires p to be evenly divisible by mm+1 and for that value to be 
    greater or equal to than k. These requirements are imposed by the use of a 
    'pigeon hole' principle dividing p into segments to query against the index. '''
    if (len(p) / (mm+1)) % 1 != 0:
        print 'p must be evenly divisible by mm+1, execution aborted.    --oneshot_approximate_matching_SubseqIndex--'
        return -1
    elif k > int(round(len(p) / (mm+1))):
        print 'k must be smaller than p / (mm+1), execution aborted.     --oneshot_approximate_matching_SubseqIndex--'
        return -1
        
    subseq_index = SubseqIndex(t, k, ival)
    return approximate_matching_SubseqIndex(subseq_index, p, t, mm, ival)
#============================================================================== 
  
#The following block of code pertains to Boyer-Moore
#============================================================================== 
def z_array(s):
    """ Use Z algorithm (Gusfield theorem 1.4.1) to preprocess s """
    assert len(s) > 1
    z = [len(s)] + [0] * (len(s)-1)

    # Initial comparison of s[1:] with prefix
    for i in range(1, len(s)):
        if s[i] == s[i-1]:
            z[1] += 1
        else:
            break

    r, l = 0, 0
    if z[1] > 0:
        r, l = z[1], 1

    for k in range(2, len(s)):
        assert z[k] == 0
        if k > r:
            # Case 1
            for i in range(k, len(s)):
                if s[i] == s[i-k]:
                    z[k] += 1
                else:
                    break
            r, l = k + z[k] - 1, k
        else:
            # Case 2
            # Calculate length of beta
            nbeta = r - k + 1
            zkp = z[k - l]
            if nbeta > zkp:
                # Case 2a: zkp wins
                z[k] = zkp
            else:
                # Case 2b: Compare characters just past r
                nmatch = 0
                for i in range(r+1, len(s)):
                    if s[i] == s[i - k]:
                        nmatch += 1
                    else:
                        break
                l, r = k, r + nmatch
                z[k] = r - k + 1
    return z


def n_array(s):
    """ Compile the N array (Gusfield theorem 2.2.2) from the Z array """
    return z_array(s[::-1])[::-1]


def big_l_prime_array(p, n):
    """ Compile L' array (Gusfield theorem 2.2.2) using p and N array.
        L'[i] = largest index j less than n such that N[j] = |P[i:]| """
    lp = [0] * len(p)
    for j in range(len(p)-1):
        i = len(p) - n[j]
        if i < len(p):
            lp[i] = j + 1
    return lp


def big_l_array(p, lp):
    """ Compile L array (Gusfield theorem 2.2.2) using p and L' array.
        L[i] = largest index j less than n such that N[j] >= |P[i:]| """
    l = [0] * len(p)
    l[1] = lp[1]
    for i in range(2, len(p)):
        l[i] = max(l[i-1], lp[i])
    return l


def small_l_prime_array(n):
    """ Compile lp' array (Gusfield theorem 2.2.4) using N array. """
    small_lp = [0] * len(n)
    for i in range(len(n)):
        if n[i] == i+1:  # prefix matching a suffix
            small_lp[len(n)-i-1] = i+1
    for i in range(len(n)-2, -1, -1):  # "smear" them out to the left
        if small_lp[i] == 0:
            small_lp[i] = small_lp[i+1]
    return small_lp


def good_suffix_table(p):
    """ Return tables needed to apply good suffix rule. """
    n = n_array(p)
    lp = big_l_prime_array(p, n)
    return lp, big_l_array(p, lp), small_l_prime_array(n)


def good_suffix_mismatch(i, big_l_prime, small_l_prime):
    """ Given a mismatch at offset i, and given L/L' and l' arrays,
        return amount to shift as determined by good suffix rule. """
    length = len(big_l_prime)
    assert i < length
    if i == length - 1:
        return 0
    i += 1  # i points to leftmost matching position of P
    if big_l_prime[i] > 0:
        return length - big_l_prime[i]
    return length - small_l_prime[i]


def good_suffix_match(small_l_prime):
    """ Given a full match of P to T, return amount to shift as
        determined by good suffix rule. """
    return len(small_l_prime) - small_l_prime[1]


def dense_bad_char_tab(p, amap):
    """ Given pattern string and list with ordered alphabet characters, create
        and return a dense bad character table.  Table is indexed by offset
        then by character. """
    tab = []
    nxt = [0] * len(amap)
    for i in range(0, len(p)):
        c = p[i]
        assert c in amap
        tab.append(nxt[:])
        nxt[amap[c]] = i+1
    return tab


class BoyerMoore(object):
    """ Encapsulates pattern and associated Boyer-Moore preprocessing. """

    def __init__(self, p, alphabet='ACGT'):
        # Create map from alphabet characters to integers
        self.amap = {alphabet[i]: i for i in range(len(alphabet))}
        # Make bad character rule table
        self.bad_char = dense_bad_char_tab(p, self.amap)
        # Create good suffix rule table
        _, self.big_l, self.small_l_prime = good_suffix_table(p)

    def bad_character_rule(self, i, c):
        """ Return # skips given by bad character rule at offset i """
        assert c in self.amap
        assert i < len(self.bad_char)
        ci = self.amap[c]
        return i - (self.bad_char[i][ci]-1)

    def good_suffix_rule(self, i):
        """ Given a mismatch at offset i, return amount to shift
            as determined by (weak) good suffix rule. """
        length = len(self.big_l)
        assert i < length
        if i == length - 1:
            return 0
        i += 1  # i points to leftmost matching position of P
        if self.big_l[i] > 0:
            return length - self.big_l[i]
        return length - self.small_l_prime[i]

    def match_skip(self):
        """ Return amount to shift in case where P matches T """
        return len(self.small_l_prime) - self.small_l_prime[1]

def exact_matching_BoyerMoore(p, p_bm, t):
    """ Do Boyer-Moore matching. p=pattern, t=text,
        p_bm=BoyerMoore object for p """
    i = 0
    occurrences = []
    aligns = 0
    comps = 0
    while i < len(t) - len(p) + 1:
        aligns += 1
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            comps += 1
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences#, aligns, comps

def approximate_matching_BoyerMoore(p, t, mm=2):
    segment_length = int(round(len(p) / (mm+1)))
    all_matches = set()
    for i in range(mm+1):
        start = i*segment_length
        end = min((i+1)*segment_length, len(p))
        p_bm = BoyerMoore(p[start:end], alphabet='ACGT')
        matches = exact_matching_BoyerMoore(p[start:end], p_bm, t)
        
        for m in matches:
            if m < start or m-start+len(p) > len(t):
                continue
            
            mismatches = 0
            for j in range(0,start):
                if p[j] != t[m-start+j]:
                    mismatches += 1
                    if mismatches > mm:
                        break
            for j in range(end, len(p)):
                if not p[j] == t[m-start+j]:
                    mismatches += 1
                    if mismatches > mm:
                        break
            
            if mismatches <= mm:
                all_matches.add(m - start)
    return list(all_matches)
#============================================================================== 

 
 

#The remaining code was used to solve specific course assignments locally
#==============================================================================
# def getFastaSeq(filedir):
#     with open(filedir, 'r') as f:
#         for line in f:
#             if line[0] == '>':
#                 seq = ''
#             else:
#                 seq += line.rstrip('\n')
#                 
#     return seq  
# dir = 'C:\Users\Skyler\Downloads\chr1.GRCh38.excerpt.fasta'
# t = getFastaSeq(dir)
#==============================================================================

##QUESTION 1 & 2
#==============================================================================
# print 'QUESTION 1 & 2'
# p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
# pos, aligns, comps = exact_matching_naive(p, t)
# print 'Alignments:', aligns
# print 'Comparisons:', comps
#==============================================================================


##QUESTION 3
#==============================================================================
# print 'QUESTION 3'
# p = 'GGCGCGGTGGCTCACGCCTGTAATCCCAGCACTTTGGGAGGCCGAGG'
# p_bm = BoyerMoore(p)
# pos, aligns, comps = exact_matching_BoyerMoore(p, p_bm, t)
# print 'Alignments:', aligns
#==============================================================================


##QUESTION 4
#==============================================================================
# print 'QUESTION 4'
# p = 'GGCGCGGTGGCTCACGCCTGTAAT'
# matches = approximate_matching_KmerIndex(p, t, 8, 2)
# print len(matches)
#==============================================================================


##QUESTION 5 & 6
#==============================================================================
# print 'QUESTION 5 & 6'
# p = 'GGCGCGGTGGCTCACGCCTGTAAT'
# finds = approximate_matching_KmerIndex(p, t, 8, 2)
# print 'Using KmerIndex, found p:', len(finds), 'times'
#==============================================================================


##QUESTION 7
#==============================================================================
# print 'QUESTION 7'
# p = 'GGCGCGGTGGCTCACGCCTGTAAT'
# finds = approximate_matching_SubseqIndex(p, t, 8, 2)
# print 'Using SubseqIndex, found p:', len(finds), 'times'
#==============================================================================
