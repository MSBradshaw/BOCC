import sys
"""
python number_edge_list.py HPO_String_edgelist.tsv HPO_String_edgelist.numbered.tsv HPO_String_edgelist.nodenames.tsv
"""
def num_edgelist(input,numbered_out,node_names_out):
    name2num={}
    n=0
    subjects = []
    predicates = []
    weights = []
    with open(input,'r') as infile:
        for line in infile:
            row = line.strip().split('\t')
            if row[0] not in name2num:
                name2num[row[0]] = n
                n += 1
            if row[1] not in name2num:
                name2num[row[1]] = n
                n += 1
            subjects.append(name2num[row[0]])
            predicates.append(name2num[row[1]])
            if len(row) == 3:
                weights.append(row[2])

    with open(numbered_out,'w') as outfile:
        for i in range(len(subjects)):
            if len(weights) > 0:
                outfile.write(str(subjects[i]) + '\t' + str(predicates[i]) + '\t' + str(weights[i]) + '\n')
            else:
                outfile.write(str(subjects[i]) + '\t' + str(predicates[i]) + '\n')


    with open(node_names_out,'w') as outfile:
        for key in name2num.keys():
            outfile.write(str(name2num[key]) + '\t' + key  + '\n')


if __name__ == "__main__":
    num_edgelist(sys.argv[1], sys.argv[2], sys.argv[3])
