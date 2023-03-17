import sys

ids = []
for line in open(sys.argv[1],'r'):
    ids.append(line.strip().split('\t')[0])

if len(ids) != len(set(ids)):
    print('Error in ',sys.argv[1])
