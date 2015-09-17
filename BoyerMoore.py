# -*- coding: utf-8 -*-
"""
Created on Wed Sep 16 11:46:23 2015

@author: Skyler
"""

import bisect

    
def bad_character_dict(p, alphabet):
    bc_dict = {}
    for e in alphabet:
        bc_dict[e] = [0] * len(p)
        for i in range(len(p)-1,-1,-1):
            offset = 0
            while offset <= i:
                if e != p[i - offset]:
                    if i - offset == 0:
                        bc_dict[e][i] = i + 1
                        break
                    offset += 1
                else:
                    bc_dict[e][i] = offset
                    break
    return bc_dict
    
def everySuffix(seq):
    suffixList = []
    for i in range(1, len(seq)):
        suffixList.append(seq[i:])
    suffixList.sort()
    return suffixList
    
def good_suffix_list(p):
    gs_list = []
    suffixList = everySuffix(p)
    for suffix in suffixList:
        found = False
        offset = 1
        while p[:-offset]:  #as long as there is some amount of p to check
            if 0 <= len(p) - len(suffix) - offset:  #when suffix fully fits on p
                substr = p[-len(suffix) - offset: -offset]
                suff = suffix
            else:  #when the suffix would hang past left end of p
                substr = p[:-offset]
                suff = suffix[len(suffix) + offset - len(p):]
            if substr == suff:
                found = True
                gs_list.append((suffix, offset))
                break
            offset += 1
        if not found:
            gs_list.append((suffix, len(p) - 3))
    return gs_list

class BoyerMoore:
    def __init__(self, p, alphabet='ACGT'):
        self.p = p
        self.alphabet = alphabet
        
        self.bc_dict = bad_character_dict(p)
        self.gs_list = good_suffix_list(p)
        
        def bad_character_rule(self, char, pos):
            return self.bc_dict['char'][pos] - 1
            
        def good_suffix_rule(self, pos):
            i = bisect.bisect_left(self.gs_list, (self.p[pos+1:], -1))
            return self.gs_list[i][1] - 1
            
        def match_skip(self):
            print 'TODO?'
    
def boyer_moore(p, p_bm, t):
    occurrences = []
    i = 0
    while i < len(t) - len(p) + 1:
        shift = 1
        mismatched = False
        for j in range(len(p)-1,-1,-1):
            if p[j] != t[i + j]:
                skip_bc = p_bm.bad_character_rule(j, t[i +j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences
