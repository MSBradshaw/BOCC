import sys
import pandas as pd

"""
1. bocc res
2. snowballing res
3. output location
"""

bocc = pd.read_csv(sys.argv[1], sep='\t')
snow = pd.read_csv(sys.argv[2], sep='\t')
# if there is no next year, the snowball data will not have anything and contains the message
# "Nada, next year does not exist", in this case, create an column for the p values that is full of nones and exit the script
for x in snow.columns:
    if 'Nada' in x:
        bocc['snowballing_pvalue'] = ['none'] * bocc.shape[0]
        bocc.to_csv(sys.argv[3],sep='\t',index=False)
        quit(0)

print(bocc)

pvals = []
# com_id  com_score  replicate_id  replicate_score  rep_and_com_size
for idx, dat in bocc.iterrows():
    sub = snow[snow['com_id'] == dat.cluster_id]
    p = None
    com_score = list(sub['com_score'])[0]
    p = 1 - (sum([com_score > x for x in sub['replicate_score']]) / sub.shape[0])
    pvals.append(p)
print(len(pvals))
print(bocc.shape)
bocc['snowballing_pvalue'] = pvals
print(bocc.shape)
bocc.to_csv(sys.argv[3],sep='\t',index=False)

