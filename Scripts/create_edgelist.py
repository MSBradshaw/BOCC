import networkx as nx
import obonet
import number_edge_list
import sys

"""
Arguments
1. hp.obo HPO network in obo format
2. genes_to_phenotype.txt the Jenkins data
3. string_edge_list.txt the string edgelist
4. Data/9606.protein.info.v11.0.txt protein name information from STRING, used to rename nodes
5. prefix to be used when saving the files
"""
hp_file = sys.argv[1]
jenkins_file = sys.argv[2]
string_file = sys.argv[3]
string_info_file = sys.argv[4]
prefix = sys.argv[5]
"""
----------------------Make the Phenotypic Abnormality Only Edge List----------------------
"""
# for the edges of HPO
H = obonet.read_obo(hp_file)

# make it a simple graph, not a multi or directed graph
h = nx.Graph()
for edge in H.edges:
    h.add_edge(edge[0], edge[1])

# keep only thing that are children of Phenotypic abnormality
root = 'HP:0000001'
children = list(h.neighbors(root))
pheno_ab = 'HP:0000118'


def get_kids_that_pass_through_node(g, target, required_node):
    keepers = set()
    for n in g.nodes():
        path = nx.shortest_path(g, target, n)
        if required_node in path:
            for x in path:
                keepers.add(x)
    return keepers


# get the nodes that are part of the Phenotypic abnormality tree
ks = get_kids_that_pass_through_node(h, root, pheno_ab)
bad_kiddos0 = get_kids_that_pass_through_node(h, root, children[0])
bad_kiddos2 = get_kids_that_pass_through_node(h, root, children[2])
bad_kiddos3 = get_kids_that_pass_through_node(h, root, children[3])
print([x for x in bad_kiddos0 if x in ks])
print([x for x in bad_kiddos2 if x in ks])
print([x for x in bad_kiddos3 if x in ks])

# make a sub graph of just the Phenotypic abnormality tree
G = nx.Graph(h.subgraph(ks))

# Rename to gene symbols
protein_mapping = {}
for line in open(string_info_file):
    row = line.strip().split()
    protein_mapping[row[0]] = row[1]

# for the edges of String
error_count = 0
fine_count = 0
not_founds = set()
founds = set()
for line in open(string_file):
    row = line.strip().split()
    n1 = row[0]
    n2 = row[1]
    try:
        n1 = protein_mapping[n1]
        founds.add(n1)
    except KeyError:
        not_founds.add(n1)

    try:
        n2 = protein_mapping[n2]
        founds.add(n2)
    except KeyError:
        not_founds.add(n2)

    G.add_edge(n1, n2)

print(len(not_founds))
print(len(founds))

# for the edges of Jenkins
for line in open(jenkins_file):
    row = line.strip().split()
    if len(row) < 4: continue
    # this will skip adding edges that include edges pruned from HPO
    if row[3] not in G.nodes: continue
    G.add_edge(row[1], row[3])

# write to a file wile excluding duplicate edges (regardless of direction) and self loops
edges = set()
with open(prefix + '.phenotypic_branch.edgelist.txt', 'w') as outfile:
    for edge in G.edges:
        row = [edge[0], edge[1]]
        row.sort()
        if row[0] == row[1]: continue
        if str(row) in edges:
            continue
        edges.add(str(row))
        outfile.write(row[0])
        outfile.write('\t')
        outfile.write(row[1])
        outfile.write('\n')

number_edge_list.num_edgelist(prefix + '.phenotypic_branch.edgelist.txt',
                              prefix + '.phenotypic_branch.numbered.edgelist.txt',
                              prefix + '.phenotypic_branch.nodenames.txt')

"""
----------------------Make the full HPO tree Edge List----------------------
"""
# for the edges of HPO
H = obonet.read_obo(hp_file)

# make it a simple graph, not a multi or directed graph
G = nx.Graph()
for edge in H.edges:
    G.add_edge(edge[0], edge[1])

# Rename to gene symbols
protein_mapping = {}
for line in open(string_info_file):
    row = line.strip().split()
    protein_mapping[row[0]] = row[1]

# for the edges of String
error_count = 0
fine_count = 0
not_founds = set()
founds = set()
for line in open(string_file):
    row = line.strip().split()
    n1 = row[0]
    n2 = row[1]
    try:
        n1 = protein_mapping[n1]
        founds.add(n1)
    except KeyError:
        not_founds.add(n1)

    try:
        n2 = protein_mapping[n2]
        founds.add(n2)
    except KeyError:
        not_founds.add(n2)

    G.add_edge(n1, n2)

print(len(not_founds))
print(len(founds))

# for the edges of Jenkins
for line in open(jenkins_file):
    row = line.strip().split()
    if len(row) < 4: continue
    G.add_edge(row[1], row[3])

# write to a file wile excluding duplicate edges (regardless of direction) and self loops
edges = set()
with open(prefix + '.all_hpo.edgelist.txt', 'w') as outfile:
    for edge in G.edges:
        row = [edge[0], edge[1]]
        row.sort()
        if row[0] == row[1]: continue
        if str(row) in edges:
            continue
        edges.add(str(row))
        outfile.write(row[0])
        outfile.write('\t')
        outfile.write(row[1])
        outfile.write('\n')

number_edge_list.num_edgelist(prefix + '.all_hpo.edgelist.txt',
                              prefix + '.all_hpo.numbered.edgelist.txt',
                              prefix + '.all_hpo.nodenames.txt')
