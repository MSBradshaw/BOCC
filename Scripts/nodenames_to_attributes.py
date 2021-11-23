import sys

"""
take a node names file as input via standard input
output standard output node number and 1 or 0 if the node in HPO for not
"""

for line in sys.stdin:
    row = line.strip().split('\t')
    if 'HP:' not in row[1]:
        print(row[0] + '\t1')
    else:
        print(row[0] + '\t0')
