from cdlib import algorithms
import networkx as nx
import argparse


def get_args():
    parser = argparse.ArgumentParser()


    parser.add_argument('--edgelist',
                        dest='edgelist',
                        required=True,
                        help='tab separated edge list')

    parser.add_argument('--output',
                        dest='output',
                        required=True,
                        help='name of file to save the community to')

    args = parser.parse_args()
    return args

args = get_args()

G = nx.read_edgelist(args.edgelist)
coms = algorithms.walktrap(G)

with open(args.output, 'w') as outfile:
    for i, com in enumerate(coms.communities):
        outfile.write('\t'.join([str(i)] + com))
        outfile.write('\n')
