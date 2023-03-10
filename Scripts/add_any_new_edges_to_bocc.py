import sys
import pandas as pd

"""
1. BOCC results
2. community file
3. new edges list
4. output, where to save an updated version of BOCC results
"""

bocc_res = pd.read_csv(sys.argv[1],sep='\t')

if '2023' not in sys.argv[3]:
    new_edges = [line.strip().split('\t') for line in open(sys.argv[3],'r')]
else:
    new_edges = None

num_new = []
com_map = {}
for line in open(sys.argv[2],'r'):
    row = line.strip().split('\t')
    com = row[1:]
    com_id = row[0]
    com_map[com_id] = com

for com_id in bocc_res['cluster_id']:
    if new_edges is None:
        num_new.append('None: no edges for 2022')
        continue
    com = com_map[str(com_id)]
    tmp_num_new = 0
    for node in com:
        tmp_num_new += sum([node in e for e in new_edges])
    num_new.append(tmp_num_new)
print(bocc_res.shape)
print(len(num_new))
bocc_res['num_new_edges_on_any_node'] = num_new

bocc_res.to_csv(sys.argv[4],sep='\t',index=False)


