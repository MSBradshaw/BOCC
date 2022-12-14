import networkx as nx

G21 = nx.read_edgelist('Edgelists/String_HPO_2021.phenotypic_branch.edgelist.txt')

# remove all genes from the graph
nodes = list(G21.nodes)
print(len(G21.nodes))
for n in nodes:
    if 'HP:' in n:
        continue
    G21.remove_node(n)
print(len(G21.nodes))

with open('Resources/hpo_node_depth.tsv','w') as outfile:
    outfile.write('node\tdepth\n')
    for n in G21.nodes:
        distance = len(nx.shortest_path(G21, source='HP:0000118', target=n))
        outfile.write('{}\t{}\n'.format(n,str(distance))) 
