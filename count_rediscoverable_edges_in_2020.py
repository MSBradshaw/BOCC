
el_2020 = []
el_2019 = []

for line in open('Edgelists/String_HPO_2020.phenotypic_branch.edgelist.txt','r'):
    row = line.strip().split('\t')
    if sum(['HP:' in x for x in row]) == 1:
        row.sort()
        el_2020.append(str(row))
el_2020 = set(el_2020)


for line in open('Edgelists/String_HPO_2019.phenotypic_branch.edgelist.txt','r'):
    row = line.strip().split('\t')
    if sum(['HP:' in x for x in row]) == 1:
        row.sort()
        el_2019.append(str(row))
el_2019 = set(el_2019)

print(len(el_2020))
print(len(el_2019))
print(len([x for x in el_2020 if x not in el_2019]))
