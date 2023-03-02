import networkx as nx

import network_creation as nm


res = nm.get_BLAST_hits_dict(blast_res1="170187_blast.csv",blast_res2="171101_blast.csv")

edges = nm.renameNodes(file1="ncbi_3_blast.tab", hits_dict=res)

# with open("spn_3_renamed.tab",'a') as file_new:
#     for e in edges:
#         file_new.write("\t".join(e)+"\n")

print(len(edges))

