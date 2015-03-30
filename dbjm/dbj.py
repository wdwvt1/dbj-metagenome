#!/usr/bin/env python

from numpy import array, int32
from collections import defaultdict

class Nodes2():

    def __init__(self):
        self.nodes = {}

    def observe(self, kmer1, kmer2):
        '''Add kmer1 to self and incriment the count of kmer1.

        This function adds kmer1 if it isn't in self and incriments the count 
        of kmer one regardless to account for connection to kmer 2. 
        '''
        if kmer1 not in self.nodes:
            self.nodes[kmer1] = array([0, 0, 0, 0], dtype=int32)
        
        self.nodes[kmer1][nt_to_ind(kmer2[-1])] += 1

    def add_from_file(self, kmer, arr):
        '''Add a node when reading from file.

        Assumes node is not already in file. Will overwrite if node already in
        file.
        '''
        self.nodes[kmer] = arr

    def abundance(self, kmer):
        '''Calculate the abundance of a kmer in self.nodes.

        This function calculates the abundance of kmer as the sum of its ingoing
        and outgoing connection abundances.
        '''
        sum_outgoing = self.nodes[kmer].sum()
        inc_nodes = incoming_nodes(kmer)
        sum_incoming = sum([self.nodes[k][nt_to_ind(kmer[-2])] for k in
                            incoming_nodes])
        return sum_outgoing + sum_incoming

def nt_to_ind(nt):
    '''Return index of a nt in the scheme actg -> 0123.'''
    if nt == 'A':
        return 0
    elif nt == 'T':
        return 1
    elif nt == 'C':
        return 2
    elif nt == 'G':
        return 3

def chop(st, k, rc=False):
    '''Chop a sequence into kmers of length k. Optionally reverse complement.

    This function yields substrings of st of length k. If rc=True, this function
    will reverse compliment st before chopping.'''
    if rc:
        st = rc_dna(st)
    for i in range(0, len(st) - (k-1)):
        yield st[i:i+k]

def rc_dna(dna_st):
    '''Reverse complement dna_st according to canonical base pairing and 'N'.

    Parameters
    ----------
    dna_st : str
        Upper case DNA string containing {actg}.

    Returns
    -------
    str
    '''
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
    return ''.join([complement[i] for i in dna_st])[::-1]

def incoming_nodes(seq, alphabet='ACTG'):
    '''Return sequences who would have an outgoing edge to seq in dbj graph.

    Paramaters
    ----------
    seq : str
        Sequence of kmer. 
    alphabet : str
        Letters of the alphabet over which kmers are formed.

    Returns
    -------
    list of strs
    '''
    return [i + seq[:-1] for i in alphabet]

def outgoing_nodes(seq, alphabet='ACTG'):
    '''Return sequences who would have an incoming edge from seq in dbj graph.

    Paramaters
    ----------
    seq : str
        Sequence of kmer. 
    alphabet : str
        Letters of the alphabet over which kmers are formed.

    Returns
    -------
    list of strs
    '''
    return [seq[1:] + i for i in alphabet]

def kmers_within_r(kmer, r, alphabet='ACTG'):
    '''Return iterator that yields kmers within a distance r of kmer.'''
    pass

def iterate_through_kmers_within_r(kmer, r, alphabet='ACTG'):
    pass




