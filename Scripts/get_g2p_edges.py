import sys

"""
Input: a tab seporated edge list as standard input
Outpit: the edges that have only 1 HPO term (edges that connect HPO and Genes)
"""

for line in sys.stdin:
    row = line.strip().split('\t')
    if sum(['HP:' in x for x in row]) == 1:
        print(line.strip())
