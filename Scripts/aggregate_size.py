import os
import sys

input_dir = sys.argv[1]
pattern = sys.argv[2]

if input_dir[-1] != '/':
    input_dir += '/'

com_id = 0
for f in os.listdir(input_dir):
    if pattern not in f: continue
    for line in open(input_dir + f):
        row = line.strip().split('\t')
        row[0] = str(com_id)
        print('\t'.join(row))
        com_id += 1
