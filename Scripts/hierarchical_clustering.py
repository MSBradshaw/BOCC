import sknetwork
import scipy
from sknetwork.hierarchy import Paris, Ward, LouvainHierarchy, cut_straight, cut_balanced
import pandas as pd
import argparse
import networkx as nx
import random
import os


def get_args():
    parser = argparse.ArgumentParser()

    parser.add_argument('--edgelist',
                        dest='edgelist',
                        required=True,
                        help='tab separated edge list')

    parser.add_argument('--output',
                        dest='output',
                        required=True,
                        help='prefix to be use in output files')

    parser.add_argument('--algo',
                        dest='algo',
                        required=True,
                        help='algorithm to be used [paris, ward, louvain]')

    parser.add_argument('--coms',
                        dest='coms',
                        required=True,
                        help='file of communities and there members, each row is a com, '
                             'first item is the com name, all others are members')

    parser.add_argument('--com_size',
                        dest='com_size',
                        required=False,
                        type=int,
                        help='the max com size')

    args = parser.parse_args()
    return args

args = get_args()
if args.algo not in ['paris', 'ward', 'louvain']:
    print('Error: --algo must be one of the following [paris, ward, louvain]')
    quit()

# load the graph
el = args.edgelist
output = args.output
coms = args.coms
algo = args.algo
# el = 'Edgelists/String_HPO_2015.phenotypic_branch.edgelist.txt'
# output = 'DELTest/infomap.String_HPO_2015.phenotypic_branch'
# coms = 'Coms/infomap.String_HPO_2015.phenotypic_branch.coms.txt'

G = nx.read_edgelist(el)
for line in open(coms,'r'):
    row = line.strip().split('\t')
    id = row[0]
    # for if the process needs to be restarted halfway through. If a line has already been processed, skip it
    if os.path.exists(output + '_' + args.algo + '_'  + id + '.{}.subcoms.txt'.format(str(args.com_size))):
        print("skipping",output + '_' + args.algo + '_'  + id + '.{}.subcoms.txt'.format(str(args.com_size)))
        continue
    else:
        print('doing', output + '_' + args.algo + '_'  + id + '.{}.subcoms.txt'.format(str(args.com_size)) )
    com = row[1:]
    if len(com) < 3:
        print('Warning com with size < 3, skipping')
        continue
    print(id)
    print(com)
    g = nx.subgraph(G, com)
    tmp_el = '.hierarchical_clustering_temp_edge_list_' + str(int(random.random()  * 100000)) + '.txt'
    with open(tmp_el,'w') as outfile:
        for edge in g.edges():
            outfile.write('\t'.join(edge))
            outfile.write('\n')
    if len(g.edges()) == 0:
        continue
    if args.algo == 'paris':
        model = Paris()
    elif args.algo == 'ward':
        model = Ward()
    elif args.algo == 'louvain':
        model = LouvainHierarchy()
    try:
        adjacency = sknetwork.data.load_edge_list(tmp_el)
        dendrogram = model.fit_transform(adjacency['adjacency'])
    except:
        print('ArpackError try to resolve')
        if len(com) <= args.com_size:
            print('Resolved because com size is < ', str(args.com_size))
            # if there was an error, but the com is under the com max size, just write it out anyways
            with open(output + '_' + args.algo + '_'  + id + '.{}.subcoms.txt'.format(str(args.com_size)), 'w') as outfile:
                outfile.write('0\t' + '\t'.join(com) )
            os.remove(tmp_el)
            exit()
        else:
            raise scipy.sparse.linalg.eigen.arpack.arpack.ArpackError
    # write the communities to a coms file
    print()
    print(dendrogram.shape)
    print(len(com))
    num_coms = min(args.com_size, dendrogram.shape[0])
    if num_coms < 2:
        print('Warning: community is to small, skipping')
        continue
    com_ids = cut_balanced(dendrogram, num_coms)
    coms_df = pd.DataFrame({'node': list(adjacency['names']), 'com': list(com_ids)})

    with open(output + '_' + args.algo + '_'  + id + '.{}.subcoms.txt'.format(str(args.com_size)), 'w') as outfile:
        for com in coms_df['com'].unique():
            sub = coms_df[coms_df['com'] == com]
            outfile.write('\t'.join([str(com)] + list(sub['node'])) + '\n')
    # delete the temp edge list
    os.remove(tmp_el)
