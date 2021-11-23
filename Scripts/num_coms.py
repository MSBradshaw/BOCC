import sys

for i,line in enumerate(open(sys.argv[1])):
    print(str(i) + '\t' + line.strip())
