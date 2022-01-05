import os
import sys

input_dir = sys.argv[1]
pattern = sys.argv[2]

if input_dir[-1] != '/':
    input_dir += '/'

for f in os.listdir(input_dir):
    if pattern not in f: continue
    for line in open(input_dir + f):
        print(line.strip())
