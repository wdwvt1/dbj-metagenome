#!/usr/bin/env python

from collections import defaultdict

class Nodes:

    def __init__(self):
        '''Initiate a Nodes object.'''
        def _f():
            return {'abundance': 0, 'reads': {'r1': [], 'r2': []}}
        self.nodes = defaultdict(_f)

    def addObservation(self, kmer, read_id, read_pair):
        '''Add observation to self.nodes

        This function adds a new node (kmer) to self if that node hasn't been 
        observed before. It adds a new observation to a known node otherwise.

        Parameters
        ----------
        kmer : str
            Sequence constituting the kmer. 
        read_id : numeric or str
            Reference to read id which is stored in Readlist class.
        read_pair : str
            One of 'r1' or 'r2'. Whether this kmer came from r1 or r2.. 
        '''
        self.nodes[kmer]['abundance'] += 1
        self.nodes[kmer]['reads'][read_pair].append(read_id)

    def __str__(self):
        return '<Nodes Object>'

    __repr__ = __str__


class ReadList:

    def __init__(self, expected_size):
        '''Create a list of reads with an expected size.'''
        self.reads = [None] * expected_size
        self.index = 0

    def addRead(self, read_id):
        '''Add a read.'''
        self.reads[self.index] = read_id
        self.index += 1

    def getLastIndex(self):
        '''Return the index of the last read.'''
        return self.index - 1

    def removeUnused(self):
        '''Remove entries which are not used because reads were discarded.'''
        self.reads = self.reads[:self.index]


def chop(st, k):
    '''produce kmers.'''
    for i in range(0, len(st) - (k-1)):
        yield st[i:i+k]

def rc_dna(dna_st):
    '''Reverse complement dna_st according to canonical base pairing.

    Parameters
    ----------
    dna_st : str
        Lower case DNA string containing {actg}.

    Returns
    -------
    str
    '''
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    return ''.join([complement[i] for i in dna_st])[::-1]
