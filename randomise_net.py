from multiprocessing import Pool
import networkx as nx
import random
import time
import pandas as pd
import numpy as np
import os
import itertools
# from matplotlib.pyplt import plt

# functions

def randomizeNetwork(G,rand_seed=100,no_swaps=0.2,graph_out = "random_test.tab"):
    """To randomize a directed graph (edge list) by using Networkx Directed edge swap."""
    graph_in = nx.read_edgelist(G,create_using=nx.DiGraph)
    nx.directed_edge_swap(graph_in,nswap=graph_in.number_of_edges()*no_swaps,max_tries=100000,seed = rand_seed)
    nx.write_edgelist(graph_in,graph_out,data=False, delimiter = '\t')
    return graph_in
    
def chunks(l, n):
    """Divide a list of nodes `l` in `n` chunks"""
    l_c = iter(l)
    while 1:
        x = tuple(itertools.islice(l_c, n))
        if not x:
            return
        yield x


def betweenness_centrality_parallel(G, processes=None):
    """Parallel betweenness centrality  function"""
    p = Pool(processes=processes)
    node_divisor = len(p._pool) * 12
    node_chunks = list(chunks(G.nodes(), G.order() // node_divisor))
    num_chunks = len(node_chunks)
    bt_sc = p.starmap(
        nx.betweenness_centrality_subset,
        zip(
            [G] * num_chunks,
            node_chunks,
            [list(G)] * num_chunks,
            [True] * num_chunks,
            [None] * num_chunks,
        ),
    )

     #Reduce the partial solutions
    bt_c = bt_sc[0]
    for bt in bt_sc[1:]:
        for n in bt:
            bt_c[n] += bt[n]
    return bt_c

