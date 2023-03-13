import matplotlib.pyplot as plt
import sys

# read in the edgelist
edges = []
with open(sys.argv[1]) as infile:
    for line in infile:
        edge = line.strip().split('\t')
        edges.append(edge)

# count he number of times each gene is in the list
gene_counts = {}
for edge in edges:
    gene = edge[1]
    if gene not in gene_counts:
        gene_counts[gene] = 0
    gene_counts[gene] += 1

# plot a histogram of the gene counts
bs = plt.hist(gene_counts.values(),bins=100)
print(bs)
plt.xlabel('Number of phenotypes per gene')
plt.ylabel('Count')
plt.savefig(sys.argv[2])
# plt.show()
# 'Figures/drug_inferred_p_per_g.png'