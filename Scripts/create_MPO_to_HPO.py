"""
The purpose of this script is to map MPO terms to HPO terms via OMIM
Most of the work for this has been cut out for me by this article
https://academic.oup.com/nar/article/38/suppl_2/W165/1116251#supplementary-data
In there supplementatl material there is a file that maps MPO to OMIM.
I have converted and saved one of those supplemental files to csv via copy and paste... (I know not great but xls documents are a crime against science)
HPO publishes an annotation file mapping HPO to OMIM.
This produces a tsv with HPO -> Gene connections that are dirived from MPO
"""

import pandas as pd
import networkx as nx
import obonet

url = 'http://purl.obolibrary.org/obo/hp.obo'
G = obonet.read_obo(url)
G_dict = G.nodes(data=True)
print(G.nodes(data=True)['HP:0001250']['name'])

sup_data = pd.read_csv('Resources/PhenoHM_supplementary_file006.csv')
print(sup_data)
print(sup_data.columns)
"""
columns:
['MPID', 'MP Term', 'Term Submitted to OMIM', 'OMIM ID',
       'OMIM Description', 'Human Gene ID', 'Human Gene Symbol',
       'Mouse Gene ID', 'Mouse Gene Symbol', 'OMIM Allelic Variant',
       'Human Gene- From OMIM Allelic Variant',
       'Mutation -From OMIM Allelic Variant', 'Score']
"""

# load HPO annotations
hp_to_omim = pd.read_csv('Resources/phenotype.hpoa',sep='\t',header=None,comment='#')
print(hp_to_omim)
# OMIM:400004      RETINITIS PIGMENTOSA, Y-LINKED  NaN  HP:0000510  OMIM:400004  TAS  NaN         NaN  NaN  NaN  P     HPO:skoehler[2017-07-13]
hp_to_omim.columns = ['OMIM ID','name','unknown1','HP ID', 'Redundant ID', 'Gene', 'unknown2','unknown3','unknown4','unknown5','unknown6','annotation_info']
print(hp_to_omim)
print(hp_to_omim.shape)
hp_to_omim = hp_to_omim[['OMIM' in x for x in hp_to_omim['OMIM ID']]]
print(hp_to_omim.shape)

# merge the two
hpo_to_merge = hp_to_omim[['HP ID','OMIM ID','Gene','name']]
mpo_to_merge = sup_data[['MPID','OMIM ID', 'Human Gene Symbol', 'Mouse Gene Symbol']]
mpo_to_merge['OMIM ID'] = ['OMIM:' + str(x) for x in mpo_to_merge['OMIM ID']]
print(hpo_to_merge['OMIM ID'])
print(mpo_to_merge['OMIM ID'])
hpo_mpo = hpo_to_merge.merge(mpo_to_merge,how='inner',on='OMIM ID')
print(hpo_mpo)
print(hpo_mpo.shape)

# how many of these HPO -> Human Gene edges are in the 2021 network?
G21 = nx.read_edgelist('Edgelists/String_HPO_2021.phenotypic_branch.edgelist.txt')
print(len(G21.nodes))
with open('Resources/hpo_to_gene_derived_by_mpo.edgelist.2021.txt','w') as outfile:
	for i,row in hpo_mpo.iterrows():
		# remove autosomal resesive and dominate inheritance terms
		if row['HP ID'] in ['HP:0003745', 'HP:0000007', 'HP:0000006']:
			continue
		if [row['HP ID'], row['Human Gene Symbol']] not in G21.edges and [row['Human Gene Symbol'], row['HP ID']] not in G21.edges:
			print(row)
			outfile.write('{hpo}\t{gene}\t{name}\n'.format(hpo=row['HP ID'],gene=row['Human Gene Symbol'], name=G.nodes(data=True)[row['HP ID']]['name']))
