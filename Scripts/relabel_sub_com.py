import sys

"""
Take a subcome file, renubmer the communities so that each one has its own unique id
1. Subcom file
2. where to save the re labeled subcom file
"""

with open(sys.argv[2],'w') as out:
    for i,line in enumerate(open(sys.argv[1],'r')):
        row = line.strip().split('\t')
        row[0] = str(i)
        out.write('\t'.join(row) + '\n')
