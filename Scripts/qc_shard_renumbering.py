import sys

"""
1. coms file
2. boc results
Print out if there is something unexpected about this results file
"""

# load the coms into a dictionary
def ss(l):
    return l.strip().split('\t')
coms = {ss(l)[0]:ss(l)  for l in open(sys.argv[1],'r')}
res = {ss(l)[0]:ss(l)  for l in open(sys.argv[2],'r')}
for key in res.keys():
    if key == 'cluster_id':
        continue
    # check the gene ratios match
    ratio = sum(['HP:' not in n for n in coms[key]]) / len(coms[key])
    if float(res[key][2]) != ratio:
        print('Unexcpected Gene Ratio Difference',' '.join([sys.argv[1],sys.argv[2],str(ratio),res[key][2]]))

