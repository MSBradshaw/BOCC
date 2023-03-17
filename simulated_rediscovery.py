import matplotlib.pyplot as plt
from BOCC.BOCC import load_clusters
import pandas as pd
import os
import numpy as np
import seaborn as sns
import networkx as nx
import random
import sys

np.random.seed(int(sys.argv[1]))

# load the 2019 communities
icoms = load_clusters('SubComs/{}/paris.{}.{}.coms.txt'.format(str(2019), 'infomap', str(2019)))
wcoms = load_clusters('SubComs/{}/paris.{}.{}.coms.txt'.format(str(2019), 'walktrap', str(2019)))
gcoms = load_clusters('SubComs/{}/paris.{}.{}.coms.txt'.format(str(2019), 'greedy', str(2019)))
ccoms = load_clusters('SubComs/{}/paris.{}.{}.coms.txt'.format(str(2019), 'cesna', str(2019)))
all_coms = {'paris.infomap': icoms, 'paris.walktrap': wcoms, 'paris.greedy': gcoms, 'paris.cesna': ccoms, }

years = [2019, 2020]
# load the the G2P Edges for 2019 and 2020
g2p_by_year = {}
for year in years:
    template = 'String_HPO_{}.phenotypic_branch.edgelist.txt'
    tmp_edges = []
    for line in open('Edgelists/' + template.format(str(year))):
        edge = line.strip().split('\t')
        edge.sort()
        # keep only the HPOs
        if sum('HP:' in x for x in edge) == 1:
            tmp_edges.append(str(edge))
    tmp_edges_set = set(tmp_edges)
    # un str the edge list
    g2p_by_year[year] = tmp_edges_set

added_in_2020 = [x.replace('[', '').replace(']', '').replace(' ', '').replace("'", '').split(',') for x in
                 g2p_by_year[2020] if x not in g2p_by_year[2019]]
added_in_2020 = set([(x[0], x[1]) for x in added_in_2020])

nodes2020 = []
for x in added_in_2020:
    nodes2020.append(x[0])
    nodes2020.append(x[1])

hpo2020 = [x for x in nodes2020 if 'HP:' in x]
genes2020 = [x for x in nodes2020 if 'HP:' not in x]

G = nx.read_edgelist('Edgelists/String_HPO_2019.phenotypic_branch.edgelist.txt')
HPOs = [n for n in G.nodes if 'HP:' in n]
Genes = [n for n in G.nodes if 'HP:' not in n]

for i in range(100):
    # choose X random g2p edges that are not edges in added_in_2020
    random_edges = set()
    edge_set = set(G.edges)
    num_added = 0
    stop_len = len(added_in_2020)
    while len(random_edges) < len(added_in_2020):
        hs = np.random.choice(HPOs, len(hpo2020))
        gs = np.random.choice(Genes, len(genes2020))
        for i in range(len(hs)):
            h = hs[i]
            g = gs[i]
            if (h, g) in edge_set or (g, h) in edge_set:
                # print('In G')
                continue
            # elif (h, g) in added_in_2020 or (g, h) in added_in_2020:
            #     print('In real list')
            #     continue
            elif (h, g) in random_edges or (g, h) in random_edges:
                # print('Already added')
                continue
            else:
                random_edges.add((h, g))
                num_added += 1
                if num_added == stop_len:
                    break

    list_random_edges = list(random_edges)
    # Score 2019 with this data
    results = {}
    results_counts = {}
    for key in all_coms.keys():
        coms = all_coms[key]
        results[key] = []
        results_counts[key] = []
        for c in coms:
            new_ones = c.get_new_edges(list_random_edges)
            results[key].append(new_ones)
            results_counts[key].append(len(new_ones))

    real_results = {}
    real_results_counts = {}
    for key in all_coms.keys():
        coms = all_coms[key]
        real_results[key] = []
        real_results_counts[key] = []
        for c in coms:
            new_ones = c.get_new_edges(added_in_2020)
            real_results[key].append(new_ones)
            real_results_counts[key].append(len(new_ones))

    for key in results_counts.keys():
        print('\t'.join(['simulation', key]+[str(x) for x in results_counts[key]]))

