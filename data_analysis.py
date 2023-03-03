import os
import random
import networkx as nx
import time
import pandas as pd
import numpy as np
import randomise_net as rn



rand_ints = open("rand_ints_0.1_test.txt","a",encoding='utf-8',newline="\n")
# generate random networks
for i in range(10000):
    rand_int = random.randint(0,1000000)
    gout = "rand_" + str(i) + "_rewired.tab"
    g0 = rn.randomizeNetwork("string_spn_largecc.tab",rand_seed=rand_int,graph_out=gout,no_swaps=0.1)
    rand_ints.write(str(rand_int)+"\n")
    
rand_ints.close()

 # add all rand nets names into a txt file or call a python list
data2 = [i for i in os.listdir() if i.endswith("_rewired.tab") is True]
print(data2)

datafile = open("completed.txt",'r',encoding='utf-8',newline="\n")
data_done = [i.rstrip() for i in datafile.readlines()]
# create empty dataframe
header = ["id"]
original_graph = nx.read_edgelist("network_test.tab", create_using=nx.DiGraph)
genes = list(original_graph.nodes(data=False))
for g in genes:
    header.append(g)
df = pd.DataFrame(columns=[header])
df.fillna(0)
print(df.shape)

# start analysis of rand networks
start = time.time() 

for netw in data2:
    if netw not in data_done:
        gr = nx.read_edgelist(netw,create_using=nx.DiGraph)
        print(f"starting betweenness centrality measure for {netw}")
        bc = rn.betweenness_centrality_parallel(gr)
        print(f"Betweenness centrality measured for {netw}")
        mes = [netw]
        for gene in list(bc.keys()):
            cent = bc[gene]
            mes.append(cent)
        df.loc[data2.index(netw)] = mes
      
        print(df.head())

# write dataframe to a file
        df.to_csv("b2_0.1_centrality.csv")        

print(f"\t\tTime: {(time.time() - start):.4F} seconds") 
