import sys

def ss(l):
    return l.strip().split('\t')

"""
Take a BOCC results file and remove all extra headers
1. Input file
2. output file
"""

with open(sys.argv[2],'w') as out:
    first = True
    for line in open(sys.argv[1],'r'):
        if line == '\n':
            continue
        row = ss(line)
        if 'cluster_id' in row[0] and first:
            out.write(line)
        elif 'cluster_id' in row[0] and not first:
            continue
        else:
            out.write(line)
        first = False
