# to read a graphml file and create an edgelist file
import networkx as nx

G = nx.read_graphml("string_spn_largecc.graphml",edge_key_type=int)

p = nx.to_pandas_edgelist(G)

sources = [i.split("(")[0].rstrip() for i in p["shared name"]]

targets = [ i.split(")")[-1].rstrip() for i in p["shared name"]]

fout = open("spn_largecc_edgelist.tab",'a',newline='\n',encoding='utf-8')
for i in range(len(sources)):
    edge = "\t".join([sources[i],targets[i]])
    print(edge)
    fout.write(edge+"\n")



fout.close()
