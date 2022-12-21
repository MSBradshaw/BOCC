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
    m2h[v] = k 
print(len(m2h)) 

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
print(len(d2g))

# create a dictionary that maps drug to disease
d2d = {}
d2d2g = []
for line in open('Resources/CTD_chemicals_diseases.csv','r'):
    if line[0] == '#':
        continue
    row = line.strip().split(',')
    drug = row[1]
    dis = row[4]
    g = row[6]
    if drug not in d2d:
        d2d[drug] = set()
    d2d[drug].add(dis)
    if drug != '' and dis != '' and g != '':
        if dis in m2h:
            d2d2g.append([m2h[dis],g,drug]) 
print(len(d2d))
print(len(d2d2g))
with open('Resources/drug_edges.txt','w') as outfile:
    for r in d2d2g:
        outfile.write('\t'.join(r) + '\n') 
