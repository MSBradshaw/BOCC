import matplotlib.pyplot as plt
import networkx as nx
import os
import seaborn as sns


plot_data = {'year': [], 'type': [], 'count': []}

for f in os.listdir('Edgelists/'):
    if 'phenotypic_branch.edgelist.txt' not in f: continue
    y = f.split('.')[0].replace('String_HPO_','')
    G = nx.read_edgelist('Edgelists/' + f)
    # break
    # get number of HPOs and Genes
    hpo_count = 0
    gene_count = 0
    for n in G.nodes:
        if 'HP:' in n:
            hpo_count += 1
        else:
            gene_count += 1

    hpo_edge_count = 0
    string_edge_count = 0
    g2p_edge_count = 0
    for e in G.edges:
        hp_sum = sum('HP:' in x for x in e)
        if hp_sum == 0:
            string_edge_count += 1
        elif hp_sum == 1:
            g2p_edge_count += 1
        elif hp_sum == 2:
            hpo_edge_count += 1
    for k,v in zip(['# HPO Nodes','# String Nodes','# HPO Edges','# G to P Edges','# String Edges'],
                   [hpo_count, gene_count, hpo_edge_count, g2p_edge_count, string_edge_count]):
        plot_data['year'].append(y)
        plot_data['type'].append(k)
        plot_data['count'].append(v)
