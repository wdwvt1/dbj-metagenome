#!/usr/bin/env python

from dbjm.dbj import Nodes, ReadList, chop
from itertools import islice
import subprocess as sub

def fastq_to_nodes(fastq, expected_num_seqs, k, r1_or_r2):
    '''Ingest a fastq file and convert it to the nodes of a DeBruijn graph. 

    Parameters
    ----------
    fastq : str
        Path to a fastq file. 
    expected_num_seqs : int
        The expected number of sequences that will be parsed. 
    k : int
        Kmer size to chop sequences into.
    r1_or_r2 : str
        'r1' or 'r2', representing which pair reads are from. 

    Returns
    -------
    nodes : Nodes object
        A Nodes object.
    rl : ReadList object
        Readlist object.
    '''
    # create a readlist
    rl = ReadList(expected_num_seqs)

    # create node storage dict
    N = Nodes()

    with open(fastq, 'r') as f:
        while True:
            vals = list(islice(f, 4))
            if len(vals) != 4:
                break
            read_id = vals[0].strip()
            seq = vals[1].strip()

            rl.addRead(read_id)
            read_index = rl.getLastIndex()
            for kmer in chop(seq, k):
                N.addObservation(kmer, read_index, r1_or_r2)

    return N, rl

def discard_low_quality(qual_score, params):
    '''discard a low quality read.'''
    raise NotImplementedError 

def nodes_to_lines(nodes):
    '''Return lines for writing to file containing information in nodes.'''
    header = '\t'.join(['kmer', 'abundance', 'r1', 'r2'])
    lines = [header]
    for node, vals in nodes.nodes.iteritems():
        r1 = ','.join(map(str, vals['reads']['r1']))
        r2 = ','.join(map(str, vals['reads']['r2']))
        line = '\t'.join([node, str(vals['abundance']), r1, r2])
        lines.append(line)
    return lines


def get_seq_from_read(read_number, readlist, reads_fp):
    '''return the sequence of a given read.'''
    cmd = "/usr/bin/grep -A 1 '%s' %s"
    seq_header = readlist.reads[read_number]
    seq = sub.Popen(cmd % (seq_header, reads_fp), stdout=sub.PIPE,
                    stderr=sub.PIPE, shell=True).communicate()[0].split('\n')[1]
    return seq





