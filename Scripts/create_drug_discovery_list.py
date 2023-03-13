import networkx as nx
import pandas as pd
import pickle
import os
import csv
import matplotlib.pyplot as plt

# create a list of gene -> hpo that pass through a drug

# load the HPO to MESH Mapping
m2h = {}
for line in open('Resources/hpo2mesh.txt','r'):
    row = line.strip().split('\t')
    k = None
    v = None
    if 'HP:' in row[0]:
        k=row[0]
        v=row[1]
    else:
        k=row[1]
        v=row[0]
    if v not in m2h:
        m2h[v] = []
    m2h[v].append(k) 
print(len(m2h)) 

# find chemicals that causes diseases
drug_causes_disease = {}
# for line in open('Resources/CTD_exposure_studies.csv','r'):
with open('Resources/CTD_exposure_studies.csv') as csvfile:
    reader = csv.reader(csvfile,  quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True)
    for row in reader:
        if row[0][0] == '#':
            continue
        # print(row)
        drugs = row[2]
        drugs = [ x.split('^')[1] for x in drugs.split('|')]
        diseases = row[7]
        if diseases == '':
            continue
        # print('"',diseases,'"')
        diseases = [ 'MESH:' + x.split('^')[1] for x in diseases.split('|')]
        for drug in drugs:
            if drug not in drug_causes_disease:
                drug_causes_disease[drug] = set()
            for dis in diseases:
                drug_causes_disease[drug].add(dis)

# create dict of drug -> gene
d2g = {}
for line in open('Resources/CTD_chem_gene_ixns.csv','r'):
    if line[0] == '#':
        continue
    row = line.strip().split(',')
    d = row[1]
    g = row[3]
    if d not in d2g:
        d2g[d] = set()
    d2g[d].add(g)

print('num d2g:',len(d2g))

# create a dictionary that maps drug to disease

d2d = {}
d2d_skips = 0
skipped_dis = set()
non_skipped_dis = set()
evidences = set()
causes_and_treats = list()
inferences = list()
inferences_by_drug = {}
# using the csv library open this file (Resources/CTD_chemicals_diseases.csv)
with open('Resources/CTD_chemicals_diseases.csv') as csvfile:
    reader = csv.reader(csvfile,  quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL, skipinitialspace=True)
    for row in reader:
        if row[0][0] == '#':
            continue
        if 'therapeutic' not in row:
            d2d_skips += 1
            skipped_dis.add(dis)
            continue
        # print(row)
        drug = row[1]
        dis = row[4]
        # if row[7] is not a float, don't store it
        if drug not in d2g:
            continue
        if row[7] != '':
            inferences.append(float(row[7]))
            if drug not in inferences_by_drug:
                inferences_by_drug[drug] = []
            inferences_by_drug[drug].append(float(row[7]))
        if 'therapeutic' in row and row[7] != '':
            print(row)
        evidences.add(row[5])
        # check if the drug also causes this disease
        if drug in drug_causes_disease and dis in drug_causes_disease[drug]:
            causes_and_treats.append([drug,dis])
            continue
        else:
            non_skipped_dis.add(dis)
        if drug not in d2d:
            d2d[drug] = set()
        d2d[drug].add(dis)


plt.hist(inferences, bins=100)
plt.xlabel('drug-disease inference score')
plt.savefig('Figures/CTD_inferences.png')
plt.clf()

plt.hist([x for x in inferences if x < 500], bins=100)
# plt.xscale('log')
plt.xlabel('drug-disease inference score')
plt.savefig('Figures/CTD_inferences.sub500.png')
plt.clf()
drug_inferences_vs_num_genes = {'drug':[],'inference':[],'num_genes':[]}
for drug in inferences_by_drug:
    if drug not in d2g:
        continue
    drug_inferences_vs_num_genes['drug'].append(drug)
    drug_inferences_vs_num_genes['inference'].append(sum(inferences_by_drug[drug])/len(inferences_by_drug[drug]))
    drug_inferences_vs_num_genes['num_genes'].append(len(d2g[drug]))

drug_inferences_vs_num_genes_df = pd.DataFrame(drug_inferences_vs_num_genes)
# sort drug_inferences_vs_num_genes_df by num_genes in descending order
drug_inferences_vs_num_genes_df = drug_inferences_vs_num_genes_df.sort_values(by='num_genes',ascending=False)
print(drug_inferences_vs_num_genes_df)


# scatter plot of drug inferences vs number of genes
plt.scatter(drug_inferences_vs_num_genes['inference'],drug_inferences_vs_num_genes['num_genes'],s=1,alpha=0.3)
plt.xlabel('mean drug-disease inference score')
plt.ylabel('number of genes')
plt.savefig('Figures/CTD_inferences_vs_num_genes.png')
plt.clf()

# plot of the number of drug-disease pairs included (y axis) if we threshold the inferences (x axis)
# inferences = sorted(inferences)
if len(inferences) > 0:
    inferences_xs = list(range(0,int(max(inferences)),int(max(inferences)/200)))
    num_pairs = []
    for i in inferences_xs:
        num_pairs.append(len([x for x in inferences if x > i]))
        print(i, num_pairs[-1])
    plt.plot(inferences_xs,num_pairs)
    plt.xlabel('drug-disease inference score (threshold)')
    plt.ylabel('number of drug-disease pairs')
    plt.yscale('log')
    # plot a horizontal line at 200,000
    plt.axhline(y=200000, color='r', linestyle='-')
    plt.savefig('Figures/CTD_inferences_vs_num_pairs.png')

print(evidences)
print('Number of drug-disease treat and cause pairs:', len(causes_and_treats))
print('num d2d:',len(d2d))
print('d2d_skips:',d2d_skips)
print('skipped_dis:',len(skipped_dis))
print('number of skipped diseases in mesh2hpo:',len(skipped_dis.intersection(set(m2h.keys()))))
print('number of non skipped diseases in mesh2hpo:',len(non_skipped_dis.intersection(set(m2h.keys()))))
# quit()

# create tios of disease, gene, drug
d2d2g = []
diseases = set()
for drug in d2d.keys():
    if drug not in d2g:
        continue
    for dis in d2d[drug]:
        for gene in d2g[drug]:
           d2d2g.append([dis,gene,drug])
           diseases.add(dis)
# print(diseases)           
print('num d2d2g:',len(d2d2g)) 
# print(m2h.keys())
# create trios of hpo, gene, drug
h2d2g = []
skips = 0
for row in d2d2g:
    dis = row[0]
    gene = row[1]
    drug = row[2]
    # print(dis,gene,drug)
    if dis not in m2h:
        skips += 1
        continue
    for hpo in m2h[dis]:
        h2d2g.append([hpo,gene,drug])
print('num h2d2g:',len(h2d2g)) 
print('skips:',skips)
       
# get unique gene disease pairs       
unique_edges = set([ str(row[:2]) for row in h2d2g])
print('g2p uniq pairs:',len(unique_edges))

with open('Resources/drug_edges.txt','w') as outfile:
    for r in d2d2g:
        outfile.write('\t'.join(r) + '\n') 

with open('Resources/drug_edges.unique.txt','w') as outfile:
    for r in unique_edges:
        outfile.write('\t'.join(r) + '\n') 

        
G22 = nx.read_edgelist('Edgelists/String_HPO_2022.phenotypic_branch.edgelist.txt')
know_drug_edges = []
new_drug_edges = []
seen = set()
for edge in h2d2g:
    edge_string = str(edge[:2])
    if edge_string in seen:
        continue
    else:
        seen.add(edge_string)
    if G22.has_edge(edge[0],edge[1]):
        know_drug_edges.append(edge)
    else:
        new_drug_edges.append(edge)

with open('Resources/new_drug_edges.txt','w') as outfile:
    for r in new_drug_edges:
        outfile.write('\t'.join(r) + '\n')

print('Number of known drug edges:',len(know_drug_edges))
print('Number of novel drug edges:',len(new_drug_edges))
print('Number of unique g2p edges:',len(seen))
# count the number of time each gene is in the list new_drug_edges
gene_counts = {}
for edge in new_drug_edges:
    gene = edge[1]
    if gene not in gene_counts:
        gene_counts[gene] = 0
    gene_counts[gene] += 1

# plot a histogram of the gene counts
plt.hist(gene_counts.values(),bins=100)
plt.xlabel('Number of phenotypes per gene')
plt.ylabel('Count')
plt.savefig('Figures/drug_inferred_p_per_g.png')
# plt.show()


# filter based on what tiff told me