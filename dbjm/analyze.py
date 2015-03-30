#!/usr/bin/env python

import matplotlib.pyplot as plt
import networkx as nx
from dbjm.dbj import outgoing_nodes, incoming_nodes, chop
from numpy import zeros, int32, int64, arange, hstack, array

def kmer_stats_from_file(fp, max_exp_count=100000):
    '''Return statistics about a list of kmers found in fp.

    This function returns the total number of kmers found, an array with bins
    containing the number of kmers at each abundance, any kmer counts above
    the max_exp_count, and the frequency of bases in the kmers.

    Parameters
    ----------
    fp : str
        Filepath to kmer file.
    max_exp_count : int
        The expected maximum abundance of any given kmer.
    '''
    num_kmers = 0
    atcg = zeros(4, dtype=int64)
    counts_bin = zeros(max_exp_count, dtype=int32)
    large_kmers = []
    o = open(fp)
    _ = o.readline()
    for line in o:
        tmp = line.strip().split('\t')
        count = sum(map(int, tmp[1:]))
        try:
            counts_bin[count] += 1
        except IndexError:
            large_kmers.append(count)
        num_kmers += 1
        atcg += [tmp[0].count(i) for i in 'ATCG']
    return counts_bin, num_kmers, large_kmers, atcg

def kmer_histogram(counts_bins, large_kmers=None, x_ub=None):
    '''Plot a histogram of kmers.'''
    xs = arange(counts_bins.shape[0])
    ys = counts_bins
    if large_kmers is not None:
        tmp = sorted(list(set(large_kmers)))
        xs_ = tmp
        ys_ = array([large_kmers.count(i) for i in tmp])
        xs = hstack((xs, xs_))
        ys = hstack((ys, ys_))
    plt.plot(xs, ys, color='blue', lw=5, alpha=.5)
    if x_ub is not None:
        plt.xlim(0, x_ub)
    plt.yscale('log')
    plt.xlabel('Abundance')
    plt.ylabel('Unique Kmers with Abundance X')
    plt.title('Kmer Counts')











def sequence_neighbor_graph(nodes, sequence, k, r, threshold, g = None):
    '''Produce a sequence neighbor graph.

    Parameters
    ----------
    nodes : Nodes object
        Nodes that will form the background to paint the sequence on. 
    sequence : str
        Sequence to trace through the graph.
    k : int
        Length of kmers to chop sequence into.
    r : int
        Radius of nodes around sequence to display. NOT USED NOW
    threshold : int
        The abundance that a node not on the sequence must have to be included.
    g : nx.DiGraph object
        networkx digraph.

    Expected number of nodes is 7n+2.
    n = (len(sequence) - k) + 1
    '''
    if g == None:
        g = nx.DiGraph()

    kmers = chop(sequence, k)

    for node in kmers:
        to_node = incoming_nodes(node) #nodes that can have edge to node
        from_node= outgoing_nodes(node) #nodes that can have edge from node

        g.add_node(node, {'abundance': nodes.nodes[node]['abundance'],
                          'is_seq': True})
        
        for n in to_node:
            if n not in g and nodes.nodes[n]['abundance'] > threshold:
                g.add_node(n, {'abundance': nodes.nodes[n]['abundance'],
                               'is_seq': False})
                
                ew = count_edge_weight(sum(nodes.nodes[n]['reads'].values(),[]),
                                       sum(nodes.nodes[node]['reads'].values(),
                                           []))
                g.add_edge(n, node, {'weight': ew})
            else:
                pass
        for n in from_node:
            if n not in g and nodes.nodes[n]['abundance'] > threshold:
                g.add_node(n, {'abundance': nodes.nodes[n]['abundance'],
                               'is_seq': False})
                ew = count_edge_weight(sum(nodes.nodes[n]['reads'].values(),[]),
                                       sum(nodes.nodes[node]['reads'].values(),
                                           []))
                g.add_edge(node, n, {'weight': ew})
            else:
                pass
    return g