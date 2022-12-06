import os
import sys
import BOCC
import networkx
import obonet
import requests
import json


"""
Provide a list of cluster files to be combined and output with HPO names rather than IDs
1. Output file name
2 - INF. the cluster files to combine and process
"""

url = 'http://purl.obolibrary.org/obo/hp.obo'
G = obonet.read_obo(url)
G_dict = G.nodes(data=True)
print(G.nodes(data=True)['HP:0001250']['name'])

output_file_name = sys.argv[1]

with open(output_file_name,'w') as outfile:
	for i,f in enumerate(sys.argv):
		if i < 2:
			continue
		# load the clusters
		coms = BOCC.load_clusters(f)
		for com in coms:
			num_HPOs = len([x for x in com.members if 'HP:' in x])
			if num_HPOs == 0 or num_HPOs == len(com.members) or len(com.members) < 3:
				continue
			members = []
			genes = []
			hpos = []
			for mem in com.members:
				if 'HP:' in mem and mem in G_dict:
					members.append(G_dict[mem]['name'])
					hpos.append(G_dict[mem]['name'])
				else:
					members.append(mem)
					genes.append(mem)
			outfile.write('{file}\t{id}\t{genes}\t{hpos}\n'.format(file = f,id=com.name,genes=','.join(genes),hpos=','.join(hpos)))
				

