import sys
import argparse


def get_args():
    parser = argparse.ArgumentParser()


    parser.add_argument('--cesna_res',
                        dest='cesna',
                        required=True,
                        help='tab separated community file output from cesna')

    parser.add_argument('--output',
                        dest='output',
                        required=True,
                        help='name of file to save the communities to')

    parser.add_argument('--node_names',
                        dest='node_name',
                        required=True,
                        help='node name file use in cesna, tab seporated, 2  columns first is the node name, second is the node number')


    args = parser.parse_args()
    return args

args = get_args()

num2node = {}
# built the name dict
for line in open(args.node_name,'r'):
    row = line.strip().split('\t')
    num2node[int(row[0])] = row[1]

with open(args.output,'w') as outfile:
    for i,line in enumerate(open(args.cesna,'r')):
        line = '\t'.join([num2node[int(x)] for x in line.strip().split('\t')]) + '\n'
        outfile.write(str(i) + '\t' + line)




