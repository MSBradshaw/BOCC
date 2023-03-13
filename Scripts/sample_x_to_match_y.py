import sys
import matplotlib.pyplot as plt
import numpy as np

# function to read in edge list 
def read_edgelist(filename):
    edges = []
    with open(filename) as infile:
        for line in infile:
            edge = line.strip().split('\t')
            edges.append(edge)
    return edges

# function to count the number of times each gene is in the list
def count_genes(edges):
    gene_counts = {}
    for edge in edges:
        gene = edge[1]
        if gene not in gene_counts:
            gene_counts[gene] = 0
        gene_counts[gene] += 1
    return gene_counts

# function to plot a histogram of the gene counts
def get_gene_counts_hist_bins(gene_counts):
    bs = plt.hist(gene_counts.values(),bins=100)
    plt.clf()
    return bs

# function to sample a gene_count dict to match the bins returned from plot_gene_counts
def sample_to_match_bins(gene_counts,bs):
    # set the random seed
    np.random.seed(0)
    # get the bin edges
    bin_edges = bs[1]
    # get the bin counts
    bin_counts = bs[0]
    # get the bin centers

    # sample the gene counts to match the bin counts
    sampled_gene_counts = {}
    for i in range(len(bin_edges)-1):
        print(bin_edges[i],bin_edges[i+1])
        print(bin_counts[i])
        # get the number of genes in the bin
        num_genes = bin_counts[i]
        # get the number of phenotypes in the bin
        num_phenotypes = bin_counts[i]
        # get the genes in the bin
        genes_in_bin = [g for g in gene_counts.keys() if gene_counts[g] >= bin_edges[i] and gene_counts[g] < bin_edges[i+1]]
        print(len(genes_in_bin))
        # sample the genes in the bin to match the bin count
        # print(genes_in_bin)
        # print(num_genes)
        sampled_genes = np.random.choice(genes_in_bin,size=int(num_genes),replace=False)
        print(len(sampled_genes))
        print()
        # add the sampled genes to the sampled_gene_counts dict
        for g in sampled_genes:
            sampled_gene_counts[g] = gene_counts[g]
    return sampled_gene_counts

# function to remove top and right border form an axis
def remove_top_right_border(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# function to plot 3 histograms of the gene counts, each in their own axis
def plot_3_gene_counts(gene_counts1,gene_counts2,gene_counts3,figname):
    # get the number of phenotypes in each gene_counts dict
    num_phenotypes1 = sum(gene_counts1.values())
    num_phenotypes2 = sum(gene_counts2.values())
    num_phenotypes3 = sum(gene_counts3.values())
    # plot the histograms
    fig, axs = plt.subplots(1,3,figsize=(20,5))
    axs[0].hist(gene_counts1.values(),bins=100)
    axs[0].set_xlabel('Number of phenotypes per gene')
    axs[0].set_ylabel('Count')
    axs[0].set_title('A                                   n={}'.format(str(num_phenotypes1)),loc='left')
    axs[1].hist(gene_counts2.values(),bins=100)
    axs[1].set_xlabel('Number of phenotypes per gene')
    axs[1].set_ylabel('Count')
    axs[1].set_title('B                                   n={}'.format(str(num_phenotypes2)),loc='left')
    axs[2].hist(gene_counts3.values(),bins=100)
    axs[2].set_xlabel('Number of phenotypes per gene')
    axs[2].set_ylabel('Count')
    axs[2].set_title('C                                   n={}'.format(str(num_phenotypes3)),loc='left')
    remove_top_right_border(axs[0])
    remove_top_right_border(axs[1])
    remove_top_right_border(axs[2])
    plt.savefig(figname)
    plt.clf()
    # plt.show()
    # 'Figures/drug_inferred_p_per_g.png'

# function to selecte the edges related to genes from the gene_counts dict and write them to a file
def select_edges(edges,gene_counts,outfile):
    # get the genes from the gene_counts dict
    genes = gene_counts.keys()
    # open the output file
    with open(outfile,'w') as out:
        # loop over the edges
        for edge in edges:
            # get the gene
            gene = edge[1]
            # check if the gene is in the gene_counts dict
            if gene in genes:
                # write the edge to the output file
                out.write('\t'.join(edge)+'\n')


# main function
def main():
    # read in the edgelist
    edges1 = read_edgelist(sys.argv[1])
    edges2 = read_edgelist(sys.argv[2])
    # count the number of times each gene is in the list
    gene_counts1 = count_genes(edges1)
    gene_counts2 = count_genes(edges2)
    # plot a histogram of the gene counts
    bs1 = get_gene_counts_hist_bins(gene_counts1)
    bs2 = get_gene_counts_hist_bins(gene_counts2)
    # sample the gene counts to match the bins returned from plot_gene_counts
    sampled_gene_counts = sample_to_match_bins(gene_counts2,bs1)
    # plot 3 histograms of the gene counts, each in their own axis
    plot_3_gene_counts(gene_counts2,gene_counts1,sampled_gene_counts,sys.argv[3])
    select_edges(edges2,sampled_gene_counts,'Resources/new_drug_edges_sampled.txt')
if __name__ == '__main__':
    main()

# run the script
# python Scripts/sample_x_to_match_y.py g2p_Edgelists/String_HPO_2022.phenotypic_branch.g2p_edgelist.txt Resources/new_drug_edges.txt Figures/drug_inferred_p_per_g.png
    