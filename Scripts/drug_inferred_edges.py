import networkx as nx
import pandas as pd

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

# load edgelist as nx graph
G = nx.read_edgelist('Edgelists/String_HPO_2022.phenotypic_branch.edgelist.txt', delimiter='\t')
genes = set([x for x in G.nodes() if 'HP:' not in x])

# get a list of drugs that cause a disease
drug_causes_disease = set()
for line in open('Resources/CTD_chemicals_diseases.csv','r'):
    if line[0] == '#':
        continue
    row = line.strip().split(',')
    if 'therapeutic' in line:
        continue
    drug_name = row[0]
    dis_id = row[4]
    drug_causes_disease.add(drug_name+'\t'+dis_id)


new_g2ps = {'genes':[],'phenotypes':[],'disease':[],'inference_score':[],'drug':[]}
num_drugs_cause_disease = 0
num_disease_not_in_hpo = 0  
num_missing_inference_scores = 0
preexisting_edges = 0
# for each line in Resources/CTD_genes_diseases.tsv
for i,line in enumerate(open('Resources/CTD_genes_diseases.tsv','r')):
    if i % 500000 == 0:
        print(i)
    if line[0] == '#':
        continue
    row = line.strip().split('\t')
    gene = row[0]
    dis_id = row[3]
    drug_name = row[5]
    if row[6] == '':
        num_missing_inference_scores += 1
        continue
    inference_score = float(row[6])
    # if this drug causes this disease skip it
    if drug_name + '\t' + dis_id in drug_causes_disease:
        num_drugs_cause_disease += 1
        continue
    # check if the gene exists
    if gene not in genes:
        continue
    # check if the disease is in the HPO
    if dis_id not in m2h:
        num_disease_not_in_hpo += 1
        continue
    for hpo in m2h[dis_id]:
        if G.has_edge(gene, hpo):
            preexisting_edges += 1
            continue
        new_g2ps['genes'].append(gene)
        new_g2ps['phenotypes'].append(hpo)
        new_g2ps['disease'].append(dis_id)
        new_g2ps['inference_score'].append(inference_score)
        new_g2ps['drug'].append(drug_name)
    
print('Number of edges skipped because the drug also causes the disease', num_drugs_cause_disease)
print('Number of edges skipped because the disease is not in the HPO', num_disease_not_in_hpo)
print('Number of edges skipped because the inference score is missing', num_missing_inference_scores)
print('Number of edges skipped because they already exist in the graph', preexisting_edges)

df = pd.DataFrame(new_g2ps)
df.to_csv('Resources/novel_drug_edges_inferred_by_CTD.tsv', sep='\t', index=False)
print(df.shape)