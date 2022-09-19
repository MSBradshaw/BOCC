import networkx as nx
import matplotlib.pyplot as plt
import argparse


def get_args():
    parser = argparse.ArgumentParser()


    parser.add_argument('--edgelist',
                        dest='edgelist',
                        required=True,
                        help='tab separated edge list')

    parser.add_argument('--all_coms',
                        dest='all_coms',
                        required=True,
                        help='file with all communities listed, first column is algorithm type, second is community id')

    parser.add_argument('--com',
                        dest='com',
                        required=True,
                        type=int,
                        help='community id of interest')

    parser.add_argument('--output',
                        dest='output',
                        required=True,
                        help='community id of interest')

    args = parser.parse_args()
    return args

args = get_args()
el = args.edgelist
all_coms = args.all_coms
com_num = args.com
output = args.output

# el = 'Data/HPO_String_edgelist_june_22_2021.tsv'
# all_coms = 'all_june_22_coms.txt'
# algo = 'infomap'
# com_num = 5
# output = 'filename2.png'

G = nx.read_edgelist(el)
com = None

for line in open(all_coms,'r'):
    row = line.strip().split('\t')
    if int(row[0]) != com_num:
        continue
    com = row[1:]

# get all neighbors of the com
list(nx.neighbors(G,com[0]))

com_neighbors = set()
for c in com:
    for n in nx.neighbors(G,c):
        com_neighbors.add(n)

g = nx.subgraph(G,com_neighbors)

colors = []

for node in g:
    if node in com:
        colors.append('green')
    else:
        colors.append('blue')

com_hpos = [n for n in g if 'HP:' in n and n in com]
com_genes = [n for n in g if 'HP:' not in n and n in com]

other_hpos = [n for n in g if 'HP:' in n and n not in com]
other_genes = [n for n in g if 'HP:' not in n and not n in com]

g_labels = {n:n for n in g.nodes}
cg = nx.subgraph(G, com)


size = 100
figure = plt.gcf()
figure.set_size_inches(8, 8)
pos = nx.spring_layout(g)
nx.draw_networkx_nodes(g,pos,nodelist=other_hpos,node_color='blue',node_shape='s',node_size=size)
nx.draw_networkx_nodes(g,pos,nodelist=other_genes,node_color='blue',node_size=size)
nx.draw_networkx_nodes(g,pos,nodelist=com_hpos,node_color='green',node_shape='s',node_size=size)
nx.draw_networkx_nodes(g,pos,nodelist=com_genes,node_color='green',node_size=size)
nx.draw_networkx_edges(g, pos, width=1.0, alpha=0.5)
nx.draw_networkx_edges(g, pos,edgelist=cg.edges, width=2.0, alpha=1, edge_color='green')
nx.draw_networkx_labels(g,pos,g_labels,font_size=10,font_color='black')
plt.savefig(output, dpi=600)
plt.clf()
