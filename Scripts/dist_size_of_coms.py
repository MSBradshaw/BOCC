import matplotlib.pyplot as plt
import BOCC
from BOCC.BOCC import load_clusters
import pandas as pd


base_dir = '/Users/michael/PycharmProjects/ClusterComparison/Coms/'
com_files = ['greedy.String_HPO_2015.phenotypic_branch.coms.txt',
             'infomap.String_HPO_2015.phenotypic_branch.coms.txt',
             'walktrap.String_HPO_2015.phenotypic_branch.coms.txt']

new_edges_file = '/Users/michael/PycharmProjects/ClusterComparison/Data/new_jenkins_edges.tsv'
new_edges = []
for line in open(new_edges_file, 'r'):
    edge = line.strip().split('\t')
    new_edges.append(edge)

plt_info = {'algo': [], 'size': [], 'new_edges': []}

fig, axes = plt.subplots(1, 3)
fig.set_size_inches(9, 3)
for i, f in enumerate(com_files):
    al = f.split('.')[0]
    f = base_dir + f
    coms = load_clusters(f)
    sizes = []
    num_new_edges = []
    for c in coms:
        plt_info['algo'].append(al)
        plt_info['size'].append(len(c.members))
        plt_info['new_edges'].append(c.get_num_new_edges(new_edges))
        sizes.append(plt_info['size'][-1])
    axes[i].hist(sizes, bins=100)
    axes[i].set_title(al)
    axes[i].set_yscale('log')
    axes[i].set_ylabel('Count')
    axes[i].set_xlabel('Community Size')
    # axes[i].set_xscale('log')
plt.tight_layout()
plt.savefig('Figures/dist_2015_walk_info_greedy.png')
plt.show()

df = pd.DataFrame(plt_info)
fig, axes = plt.subplots(1, 3)
fig.set_size_inches(9, 3)
for i, al in enumerate(df.algo.unique()):
    sub = df[df['algo'] == al]
    axes[i].scatter(sub['size'], sub['new_edges'])
    axes[i].set_title(al)
    axes[i].set_xlabel('Community Size')
    axes[i].set_ylabel('Num. Rediscovered Edges')
plt.tight_layout()
plt.savefig('Figures/size_vs_new_edges_2015_walk_info_greedy.png')
plt.show()

df = pd.DataFrame(plt_info)
fig, axes = plt.subplots(1, 3)
fig.set_size_inches(9, 3)
for i,al in enumerate(df.algo.unique()):
    sub = df[df['algo'] == al]
    axes[i].scatter(sub['size'], sub['new_edges'])
    axes[i].set_title(al)
    axes[i].set_xlabel('Community Size')
    axes[i].set_ylabel('Num. Rediscovered Edges')
    axes[i].set_yscale('log')
plt.tight_layout()
plt.savefig('Figures/size_vs_new_edges_log_2015_walk_info_greedy.png')
plt.show()


"""
Same stuff as above but with the heirarchical models and crosses
"""


base_dir = '/Users/michael/PycharmProjects/ClusterComparison/AggSubComsBySize/'
com_files = ['greedy_paris_size_100.txt','infomap_paris_size_100.txt','walktrap_paris_size_100.txt','paris_paris_size_100.coms.txt']

new_edges_file = '/Users/michael/PycharmProjects/ClusterComparison/Data/new_jenkins_edges.tsv'
new_edges = []
for line in open(new_edges_file, 'r'):
    edge = line.strip().split('\t')
    new_edges.append(edge)

plt_info = {'algo': [], 'size': [], 'new_edges': []}

fig, axes = plt.subplots(1, 4)
fig.set_size_inches(9, 3)
for i, f in enumerate(com_files):
    al = f.split('.')[0].split('_')[0]
    f = base_dir + f
    coms = load_clusters(f)
    sizes = []
    num_new_edges = []
    for c in coms:
        plt_info['algo'].append(al)
        plt_info['size'].append(len(c.members))
        plt_info['new_edges'].append(c.get_num_new_edges(new_edges))
        sizes.append(plt_info['size'][-1])
    axes[i].hist(sizes, bins=10)
    axes[i].set_title(al)
    # axes[i].set_yscale('log')
    axes[i].set_ylabel('Count')
    axes[i].set_xlabel('Community Size')
    # axes[i].set_xscale('log')
plt.tight_layout()
plt.savefig('Figures/dist_2015_paris_cross_walk_info_greedy.png')
plt.show()

df = pd.DataFrame(plt_info)
fig, axes = plt.subplots(1, 4)
fig.set_size_inches(9, 3)
for i,al in enumerate(df.algo.unique()):
    sub = df[df['algo'] == al]
    axes[i].scatter(sub['size'], sub['new_edges'])
    axes[i].set_title(al)
    axes[i].set_xlabel('Community Size')
    axes[i].set_ylabel('Num. Rediscovered Edges')
plt.tight_layout()
plt.savefig('Figures/size_vs_new_edges_2015_paris_crosses_walk_info_greedy.png')
plt.show()

"""
Same as above but with cesna from 2020!!!!
"""


base_dir = './'
com_files = ['paris.cesna.100.coms.txt', 'paris.walktrap.100.coms.txt', 'paris.greedy.100.coms.txt']

new_edges_file = '/Users/michael/PycharmProjects/ClusterComparison/Data/new_jenkins_edges.tsv'
new_edges = []
for line in open(new_edges_file, 'r'):
    edge = line.strip().split('\t')
    new_edges.append(edge)

plt_info = {'algo': [], 'size': [], 'new_edges': []}

fig, axes = plt.subplots(1, 3)
fig.set_size_inches(9, 3)
for i, f in enumerate(com_files):
    al = f.split('.')[1]
    f = base_dir + f
    coms = load_clusters(f)
    sizes = []
    num_new_edges = []
    for c in coms:
        plt_info['algo'].append(al)
        plt_info['size'].append(len(c.members))
        plt_info['new_edges'].append(c.get_num_new_edges(new_edges))
        sizes.append(plt_info['size'][-1])
    axes[i].hist(sizes, bins=10)
    axes[i].set_title(al)
    # axes[i].set_yscale('log')
    axes[i].set_ylabel('Count')
    axes[i].set_xlabel('Community Size')
    # axes[i].set_xscale('log')
plt.tight_layout()
plt.savefig('Figures/cesna2020.png')
plt.show()

df = pd.DataFrame(plt_info)
fig, axes = plt.subplots(1, 3)
fig.set_size_inches(9, 3)
for i,al in enumerate(df.algo.unique()):
    sub = df[df['algo'] == al]
    axes[i].scatter(sub['size'], sub['new_edges'])
    axes[i].set_title(al)
    axes[i].set_xlabel('Community Size')
    axes[i].set_ylabel('Num. Rediscovered Edges')
plt.tight_layout()
plt.savefig('Figures/cesna2020_scatter.png')
plt.show()