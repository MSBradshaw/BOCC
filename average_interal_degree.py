from cdlib import evaluation
from cdlib import NodeClustering
from BOCC.BOCC import load_clusters
import networkx as nx
import sys

"""
calculate the average embeddedness for each com 
1. com file
2. graph edgelist
"""

coms = load_clusters(sys.argv[1])
G = nx.read_edgelist(sys.argv[2])

print('com_id\tsize\taverage_internal_degree')
for com in coms:
    communities = NodeClustering([com.members], graph=G, method_name="nada_nada_limonada")
    mod = evaluation.average_internal_degree(G, communities)
    print(str(com.name) + '\t' + str(len(com.members)) + '\t' +str(mod.score))
