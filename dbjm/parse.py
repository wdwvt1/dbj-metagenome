#!/usr/bin/env python

from dbjm.dbj import Nodes, Nodes2, ReadList, chop
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

def fastq_to_nodes2(fastq, k):
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


            kmers = list(chop(seq, k))
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

def nodes_to_lines_mem(nodes, fp):
    '''Write nodes to lines one at a time to avoid storing in memory.'''
    o = open(fp, 'w')
    header = '\t'.join(['kmer', 'abundance', 'r1', 'r2'])
    o.write(header+'\n')

    for node, vals in nodes.nodes.iteritems():
        r1 = ','.join(map(str, vals['reads']['r1']))
        r2 = ','.join(map(str, vals['reads']['r2']))
        line = '\t'.join([node, str(vals['abundance']), r1, r2])
        o.write(line+'\n')

    o.close()

def nodes_from_lines(fp, ignore_r1=False, ignore_r2=False):
    '''Return a Nodes object built from fp.'''
    n = Nodes()
    o = open(fp, 'U')
    _ = o.readline() # get rid of line0 header, can i do this with seek?
    for line in o:
        if ignore_r1:
            kmer, abd, r2s = line.strip().split('\t')
            n.addObservationFromLine(kmer, abd, [], r2s.strip().split(','))
        elif ignore_r2:
            kmer, abd, r1s = line.strip().split('\t')
            n.addObservationFromLine(kmer, abd, r1s.strip().split(','), [])
    o.close()
    return n


def get_seq_from_read(read_number, readlist, reads_fp):
    '''return the sequence of a given read.'''
    cmd = "/usr/bin/grep -A 1 '%s' %s"
    seq_header = readlist.reads[read_number]
    seq = sub.Popen(cmd % (seq_header, reads_fp), stdout=sub.PIPE,
                    stderr=sub.PIPE, shell=True).communicate()[0].split('\n')[1]
    return seq

def node2_to_file(nodes, fp):
    '''store node2 to file.'''
    o = open(fp, 'w')
    header = '\t'.join(['kmer', 'A', 'T', 'C', 'G'])
    o.write(header+'\n')
    for kmer, edge_vals in nodes.nodes.iteritems():
        line = '\t'.join(map(str, [kmer] + list(edge_vals)))
        o.write(line+'\n')
    o.close()





