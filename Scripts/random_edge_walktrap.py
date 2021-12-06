from cdlib import algorithms
import networkx as nx
import argparse
import random

def get_args():
    parser = argparse.ArgumentParser()


    parser.add_argument('--edgelist',
                        dest='edgelist',
                        required=True,
                        type=int,
                        help='tab separated edge list')

    parser.add_argument('--percent_random_edges',
                        dest='percent_random_edges',
                        required=True,
                        help='int, percent of total edges to be added to the network as random edges')

    parser.add_argument('--output',
                        dest='output',
                        required=True,
                        help='name of file to save the community to')

    args = parser.parse_args()
    return args

args = get_args()

G = nx.read_edgelist(args.edgelist)
# get list of genes
# get list of HPOs
# choose x at random that are not connected
genes = [x for x in G.nodes if 'HP:' not in x]
hpos = [x for x in G.nodes if 'HP:' in x]
random_edges = []
per = args.percent_random_edges
per = 10
while len(random_edges) < (len(10)/100) * len(G.edges):
    g = random.choice(genes)
    h = random.choice(hpos)
    if (h,g) in G.edges or (g,h) in G.edges:
        continue
    else:
        G.add_edge(g,h)

coms = algorithms.walktrap(G)

with open(args.output, 'w') as outfile:
    for i, com in enumerate(coms.communities):
        outfile.write('\t'.join([str(i)] + com))
        outfile.write('\n')
