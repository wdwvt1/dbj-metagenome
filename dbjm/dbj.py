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

    def addObservationFromLine(self, kmer, abundance, r1s, r2s):
        '''Add an observation to self.nodes when total abundance etc. is known.

        This function adds a new node and all the abundance and read presence
        information for this kmer. To be used when adding nodes from a txt file.

        Parameters
        ----------
        kmer : str
            Sequence constituting the kmer. 
        abundance : int
            Number of observations of this kmer. 
        r1s : list
            Read ids from read 1 that contained this kmer. 
        r2s : list
            Read ids from read 2 that contained this kmer.
        '''
        self.nodes[kmer]['abundance'] = abundance
        self.nodes[kmer]['reads']['r1'] = r1s
        self.nodes[kmer]['reads']['r2'] = r2s

    def removeLowCountNodes(self, threshold):
        '''Remove nodes whose abundance < threshold.'''
        n = self.nodes.keys() #can't change iterable length in for loop
        for k in n:
            if self.nodes[k]['abundance'] < threshold:
                self.nodes.pop(k)

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
    '''Reverse complement dna_st according to canonical base pairing and 'N'.

    Parameters
    ----------
    dna_st : str
        Lower case DNA string containing {actg}.

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

def count_edge_weight(node1_reads, node2_reads):
    '''count number of reads linking node1 and node2.'''
    e = 0 
    for read in set(node1_reads):
        if read in node2_reads:
            e +=1
    return e

def calculate_path_coverage(nodes, sequence, k, read_pair):
    '''Calculate coverage of sequence.'''
    kmers = list(chop(sequence, k))
    ews = []
    for ind in range(len(kmers) - 1):
        ew = count_edge_weight(nodes.nodes[kmers[ind]]['reads'][read_pair],
                               nodes.nodes[kmers[ind+1]]['reads'][read_pair])
        ews.append(ew)
    return ews







