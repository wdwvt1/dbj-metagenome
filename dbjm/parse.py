#!/usr/bin/env python

from dbjm.dbj import Nodes, Nodes2, ReadList, chop
from itertools import islice
import subprocess as sub
from numpy import array, int32

def fastq_to_nodes2(fastq, k, rc=False):
    ''''''
    # create node storage dict
    N = Nodes2()

    with open(fastq, 'r') as f:
        while True:
            vals = list(islice(f, 4))
            if len(vals) != 4:
                break
            read_id = vals[0].strip()
            seq = vals[1].strip()
            kmers = list(chop(seq, k, rc))
            i = 0
            j = 1
            while j < len(kmers):
                N.observe(kmers[i], kmers[j])
                i += 1
                j += 1
    return N


def discard_low_quality(qual_score, params):
    '''discard a low quality read.'''
    raise NotImplementedError 

def node2_to_file(nodes, fp):
    '''store node2 to file.'''
    o = open(fp, 'w')
    header = '\t'.join(['kmer', 'A', 'T', 'C', 'G'])
    o.write(header+'\n')
    for kmer, edge_vals in nodes.nodes.iteritems():
        line = '\t'.join(map(str, [kmer] + list(edge_vals)))
        o.write(line+'\n')
    o.close()

def nodes_from_file(fp, threshold=None):
    '''Create Nodes object from filepath.'''
    n = Nodes2()
    o = open(fp)
    _ = o.readline()
    for line in o.readline()
        tmp = line.strip().split('\t')
        counts = map(int, tmp[1:]))
        if threshold is not None:
            if sum(counts) >= threshold:
                n.add_from_file(tmp, array(counts, dtype=int32))
        else:
            n.add_from_file(tmp, array(counts, dtype=int32))
    o.close()
    return n

def merge_nodes_files(fp1, fp2, out_fp):
    '''Merge two nodes objects and write to file.'''
    nodes1 = nodes_from_file(fp1)
    o = open(fp2)
    _ = o.readline()
    for line in o.readline()
        tmp = line.strip().split('\t')
        counts = map(int, tmp[1:]))
        if tmp[0] in nodes1.nodes:
            nodes1.nodes[tmp0] += counts
        else:
            nodes1.add_from_file(tmp[0], array(counts))
    o.close()
    node2_to_file(nodes1, out_fp)

