import sys

"""
1. subcom file
2. output subcom file
"""

subcoms = {}
com_id = 0
sub_com_rows = []
non_trival_counts = 0
# for line in the subcoms
for i,line in enumerate(open(sys.argv[1], 'r')):
    row = line.strip().split('\t')
    row[0] = str(com_id)
    sub_com_rows.append(row)
    subcoms[com_id] = row[1:]
    com_id += 1

# write the subcoms
with open(sys.argv[2],'w') as outfile:
    for i in range(len(sub_com_rows)):
        outfile.write('\t'.join(sub_com_rows[i]) + '\n')
















