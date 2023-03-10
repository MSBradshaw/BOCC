import networkx as nx
import argparse
import random
import pandas as pd

"""
python snowball.py --edgelist edgelists/String_HPO_2015.phenotypic_branch.edgelist.txt
--output snowball.infomap.String_HPO_2015.phenotypic_branch.tsv
--coms Coms/infomap.String_HPO_2015.phenotypic_branch.coms.txt
--new_edges Data/new_jenkins_edges.tsv
--reps 100
"""

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

    parser.add_argument('--new_edges',
                        dest='new_edges',
                        required=True,
                        help='file with new edges, tab separated')

    parser.add_argument('--coms',
                        dest='coms',
                        required=True,
                        help='file of communities and there members, each row is a com, '
                             'first item is the com name, all others are members')

    parser.add_argument('--reps',
                        dest='reps',
                        required=True,
                        default=100,
                        type=int,
                        help='number of repitions, default = 100')

    args = parser.parse_args()
    return args

def score_com(G, com, edges_dict):
    """

    :param com: list of members of a community
    :param edges_dict: dictionary of edges keys are nodes, values are list of neighbors
    :return: number of edges from the dict found in the com
    """
    count = 0
    sub = nx.subgraph(G, com)
    for n1 in com:
        for n2 in com:
            e = [n1,n2]
            try:
                if e[1] in edges_dict[e[0]] or e[0] in edges_dict[e[1]]:
                    count += 1
            except KeyError:
                continue
    return count


args = get_args()

el = args.edgelist
output = args.output
coms = args.coms
new_edge_list = args.new_edges
reps = args.reps

# el = 'Edgelists/String_HPO_2015.phenotypic_branch.edgelist.txt'
# output = 'del.txt'
# coms = 'Coms/infomap.String_HPO_2015.phenotypic_branch.coms.txt'
# new_edge_list = 'Data/new_jenkins_edges.tsv'
# reps = 100

G = nx.read_edgelist(el)

new_edges = {}
new_edges_count = 0
for line in open(new_edge_list,'r'):
    row = line.strip().split('\t')
    # remove edges already in the graph
    if [row[0], row[1]] in G.edges or [row[1], row[0]] in G.edges:
        continue
    if row[0] not in new_edges:
        new_edges[row[0]] = []
    if row[1] not in new_edges:
        new_edges[row[1]] = []
    new_edges[row[1]].append(row[0])
    new_edges[row[0]].append(row[1])
    new_edges_count += 1
print(new_edges_count)

data = {'com_id':[], 'com_score':[], 'replicate_id':[], 'replicate_score':[], 'rep_and_com_size':[]}

# for each com
for i,line in enumerate(open(coms, 'r')):
    row = line.strip().split('\t')
    com = row[1:]
    # this is one of the proteins that is in the weird connected component of size 2, it causes errors if not ignored
    if '9606.ENSP00000411694' in com: continue
    com_score = score_com(G, com, new_edges)
    for i in range(reps):
        # choose a random node from the community
        start = random.sample(list(G.nodes), 1)[0]
        # snowball till we have a set of size len(com)
        snowball = set()
        snowball_genes = set()
        snowball_hpos = set()
        used = set()
        snowball.add(start)
        current = start
        finished = set()
        caught_error = False
        com_hpos = [x for x in com if 'HP:' in x]
        com_genes = [x for x in com if 'HP:' not in x]
        print('starting while')
        itter = 0
        while len(snowball_genes) + len(snowball_hpos) < len(com):
            print(itter)
            itter += 1
            not_used = [x for x in snowball if x not in finished]
            try:
                current = random.sample(not_used, 1)[0]
            except ValueError as er:
                print('Warning, error, com')
                print(er)
                print('Not used',str(not_used))
                print('Used',str(finished))
                print('Cluster',str(snowball))
                caught_error = True
                print(com)
                break
            neighs = list(G.neighbors(current))
            deficit = len(com) - len(snowball_genes) + len(snowball_hpos)
            choose_x = min(deficit, len(neighs))
            choices = random.sample(neighs, choose_x)
            if len(choices) == len(neighs):
                finished.add(current)
            for x in choices:
                if 'HP:' in x:
                    if len(snowball_hpos) < len(com_hpos):
                        snowball_hpos.add(x)
                    else:
                        finished.add(current)
                        print('HP else',len(snowball_hpos), len(com_hpos))
                else:
                    if len(snowball_genes) < len(com_genes):
                        snowball_genes.add(x)
                    else:
                        finished.add(current)
                        print('Gene else',len(snowball_genes), len(com_genes))
            snowball = set(list(snowball_hpos) + list(snowball_genes))
        snowball = list(snowball_hpos) + list(snowball_genes)
        snowball_score = score_com(G, list(snowball), new_edges)
        data['com_id'].append(row[0])
        data['com_score'].append(com_score)
        data['replicate_id'].append(i)
        data['replicate_score'].append(snowball_score)
        data['rep_and_com_size'].append(len(com))
        if not caught_error:
            print(len(snowball_hpos), len(com_hpos))
            print(len(snowball_genes), len(com_genes))
            assert(len(snowball_hpos) == len(com_hpos))


df = pd.DataFrame(data)

df.to_csv(output, index=False, sep='\t')

