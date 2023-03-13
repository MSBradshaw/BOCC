import networkx as nx

G22 = nx.read_edgelist('Edgelists/String_HPO_2022.phenotypic_branch.edgelist.txt')

# remove all genes from the graph
nodes = list(G22.nodes)
print(len(G22.nodes))
for n in nodes:
    if 'HP:' in n:
        continue
    G22.remove_node(n)
print(len(G22.nodes))

with open('Resources/hpo_node_depth.2022.tsv','w') as outfile:
    outfile.write('node\tdepth\n')
    for n in G22.nodes:
        distance = len(nx.shortest_path(G22, source='HP:0000118', target=n))
        outfile.write('{}\t{}\n'.format(n,str(distance))) 
