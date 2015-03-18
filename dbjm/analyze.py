#!/usr/bin/env python

import matplotlib.pyplot as plt
import networkx as nx
from dbjm.dbj import outgoing_nodes, incoming_nodes, chop, count_edge_weight

def node_abundance_histogram(nodes):
    '''Plot a node abundance histogram.'''
    ab = [v['abundance'] for v in nodes.nodes.values()]
    counts, bins, patches = plt.hist(ab, bins=250, log=True)
    plt.xlabel('Unique Kmers')
    plt.ylabel('Count')
    plt.title('Kmer Count Histogram')
    return counts, bins, patches

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
    g : nx.DiGraph object
        networkx digraph.

    Expected number of nodes is 2*7 + (n-2)*6 +n.
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