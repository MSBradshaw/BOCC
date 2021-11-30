import matplotlib.pyplot as plt
import networkx as nx
import os
import seaborn as sns
import pandas as pd

plot_data = {'year': [], 'types': [], 'count': []}
num_jenkin_edges = {'year': [], 'count': []}
inferred_count = {'year': [], 'count': []}
for f in os.listdir('Edgelists/'):
    if 'phenotypic_branch.edgelist.txt' not in f: continue
    y = f.split('.')[0].replace('String_HPO_', '')
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
    hpo_edges = set()
    string_edge_count = 0
    string_edges = set()
    g2p_edge_count = 0
    g2p_edges = set()
    for e in G.edges:
        hp_sum = sum('HP:' in x for x in e)
        if hp_sum == 0:
            string_edge_count += 1
            string_edges.add(e)
        elif hp_sum == 1:
            g2p_edge_count += 1
            g2p_edges.add(e)
        elif hp_sum == 2:
            hpo_edge_count += 1
            hpo_edges.add(e)
    hpo_edge_count = len(hpo_edges)
    g2p_edge_count = len(g2p_edges)
    string_edge_count = len(string_edges)
    for k, v in zip(['# HPO Nodes', '# String Nodes', '# HPO Edges', '# G to P Edges', '# String Edges'],
                    [hpo_count, gene_count, hpo_edge_count, g2p_edge_count, string_edge_count]):
        plot_data['year'].append(int(y))
        plot_data['types'].append(k)
        plot_data['count'].append(v)
    # get the actual counts with filtering
    num_jenkin_edges['year'].append(y)
    num_jenkin_edges['count'].append(g2p_edge_count)
    # get the inferred number of edges
    year = str(y)
    edges = set()
    for line in open('work/{}/inferred_genes_to_phenotype.txt'.format(year), 'r'):
        row = line.strip().split('\t')
        # this will skip adding edges that include edges pruned from HPO
        if row[3] not in G.nodes: continue
        e = (row[1], row[3])
        edges.add(str(e))
    inferred_count['count'].append(len(edges))
    inferred_count['year'].append(str(year))

plotting_df = pd.DataFrame(plot_data)
sns.barplot(data=plotting_df, y='count', x='year', hue='types')
# sns.catplot(x=list(plot_data['count']),y=list(plot_data['year']))
plt.yscale('log')
# plt.show()
plt.savefig('Figures/network_stats_across_years_barplot.png')
plt.show()

tmp = plotting_df[plotting_df['types'] != '# String Edges']
sns.barplot(data=tmp, y='count', x='year', hue='types')
# sns.catplot(x=list(plot_data['count']),y=list(plot_data['year']))
# plt.yscale('log')
# plt.show()
plt.savefig('Figures/network_stats_across_years_barplot_without_string.png')
plt.show()

"""
Plot inferred vs actual
"""



jenkins_df = pd.DataFrame(num_jenkin_edges)
inferred_df = pd.DataFrame(inferred_count)
joint_df = pd.merge(jenkins_df, inferred_df, on='year')

year_color_map = {'2015': 'red', '2016': 'orange', '2017': 'yellow', '2018': 'green', '2019': 'blue', '2020': 'purple',
                  '2021': 'black'}
joint_df['year'] = joint_df['year'].astype(str)
for year in year_color_map.keys():
    sub = joint_df[joint_df['year'] == year]
    plt.scatter(sub['count_x'], sub['count_y'], color=year_color_map[str(year)], label=str(year))
# plt.show()
plt.legend()
plt.plot([20000, 200000], [20000, 200000], '--')
plt.xlim((20000, 200000))
plt.ylim((20000, 200000))
plt.xlabel('Actual Num. g2p Edges')
plt.ylabel('Inferred Num. g2p Edges')
plt.tight_layout()
plt.savefig('Figures/actual_vs_inferred_num_edges.png')
plt.show()
