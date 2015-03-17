#!/usr/bin/env python

from skbio.parse.sequences import parse_fastq
from dbj import Node, ReadList, chop



def create_components(fastq, expected_size, k, r1_or_r2):

    # create a readlist
    rl = ReadList(expected_size)

    # create node storage dict
    nodes = {}

    for read_id, seq, qual in parse_fastq(fastq):
        rl.addRead(read_id)
        read_index = rl.getLastIndex()

        for kmer in chop(seq, k):
            if kmer not in nodes: #node does not exist
                nodes[kmer] = Node(kmer)
            else: #node alread exists
                pass
            nodes[kmer].addObservation()
            nodes[kmer].addContainingRead(read_index, r1_or_r2)

    return nodes, rl


def discard_low_quality(qual_score, params):
    '''discard a low quality read.'''
    raise NotImplementedError 





