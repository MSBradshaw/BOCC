import networkx as nx
import argparse
import infomap


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

    parser.add_argument('--nodenames',
                        dest='names',
                        required=True,
                        help='file with node names and number separated by a tab')

    args = parser.parse_args()
    return args


def find_communities(G):
    """
    Partition network with the Infomap algorithm.
    Annotates nodes with 'community' id.
    """

    im = infomap.Infomap("--two-level")

    print("Building Infomap network from a NetworkX graph...")

    im.add_networkx_graph(G)

    print("Find communities with Infomap...")
    im.run()

    print(f"Found {im.num_top_modules} modules with codelength: {im.codelength}")

    communities = im.get_modules()
    nx.set_node_attributes(G, communities, 'community')
    print(type(communities))
    print(communities)
    return communities

args = get_args()
# el = 'Data/String_HPO_2021.edgelist'
# out = 'String_HPO_2021.infomap.coms.july_9_2021.txt'
# names = 'Data/String_HPO_2021.nodenames.tsv'
el = args.edgelist
out = args.output
names = args.names
num2name = {}
for line in open(names,'r'):
    row = line.strip().split('\t')
    num2name[int(row[0])] = row[1]

#args.edgelist
G = nx.read_edgelist(el)
coms = find_communities(G)
coms_dict = {}
for key in coms.keys():
    if coms[key] in coms_dict:
        coms_dict[coms[key]].append(num2name[key])
    else:
        coms_dict[coms[key]] = [num2name[key]]
#args.output
with open(out, 'w') as outfile:
    for i, key in enumerate(coms_dict.keys()):
        com = coms_dict[key]
        outfile.write('\t'.join([str(i)] + com))
        outfile.write('\n')
