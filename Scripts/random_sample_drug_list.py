import sys
import matplotlib.pyplot as plt
import numpy as np
import networkx as nx

# load an edgelist as a nx graph
def load_edgelist(filename):
    G = nx.read_edgelist(filename, delimiter='\t')
    return G

# take a list of edges and return just those that are not already in the graph
def get_new_edges(edges,G):
    new_edges = []
    for edge in edges:
        if not G.has_edge(edge[0],edge[1]):
            new_edges.append(edge)
    return new_edges

# function to read in edge list 
def read_edgelist(filename):
    edges = []
    with open(filename) as infile:
        for line in infile:
            edge = line.strip().split('\t')[0:2]
            edges.append(edge)
    return edges

def select_edges(edges,num,outfile):
    # select X edges at random from the list
    np.random.seed(0)
    edges_numbered = range(len(edges))
    selected_edges_indexes = np.random.choice(edges_numbered,num,replace=False)
    selected_edges = [edges[i] for i in selected_edges_indexes]
    # write the selected edges to a file
    with open(outfile,'w') as outfile:
        for edge in selected_edges:
            outfile.write('\t'.join(edge)+'\n')

# main function
def main():
    # read in the edgelist
    edges1 = read_edgelist(sys.argv[1])
    edges2 = read_edgelist(sys.argv[2])
    # load the graph
    G = load_edgelist(sys.argv[3])
    # filter edges
    edges1 = get_new_edges(edges1, G)
    select_edges(edges2,len(edges1),'Resources/randomly_sampled_edge.proportional.txt')
if __name__ == '__main__':
    main()

# run the script
# python Scripts/random_sample_drug_list.py g2p_Edgelists/String_HPO_2022.phenotypic_branch.g2p_edgelist.txt Resources/novel_drug_edges_inferred_by_CTD.tsv Edgelists/String_HPO_2021.phenotypic_branch.edgelist.txt
    