import pandas as pd
import datetime
import argparse

"""
The purpose of this file to to infer the 2018 genes_to_phenotype.txt file from the current ones by using the date in
the Biocuration of the .hpoa file. This was what was recommended to me by Peter Robinson, a PI on the HPO project.
"""

anno_dict = {'reference':[],"HPO-ID":[],'curators':[]}

for line in open('Resources/phenotype_annotation.tab','r'):
    if line[0] == '#': continue
    row = line.strip().split('\t')
    anno_dict['reference'].append(row[5])
    anno_dict['HPO-ID'].append(row[4])
    anno_dict['curators'].append(row[12])
    if 'HP' not in anno_dict['HPO-ID'][-1]:
        print(line)
        print(row)
        print()

anno = pd.DataFrame(anno_dict)
#print(anno)

anno['date'] = [x.split('[')[1].split(']')[0] if '[' in str(x) else '2999-01-01' for x in anno['curators']]

dates = [ int(datetime.datetime.strptime(x, '%Y-%m-%d').strftime("%s")) for x in anno['date']]

cut_off_time = int(datetime.datetime.strptime('{}-01-01'.format('2022'), '%Y-%m-%d').strftime("%s"))
#print(dates)
for y in [str(x) for x in range(2015,2023)]:
    cut_off_time = int(datetime.datetime.strptime('{}-01-01'.format(y), '%Y-%m-%d').strftime("%s"))
    print(y, sum([ x < cut_off_time for x in dates]))

ds = list(set(dates))
ds.sort()
print(ds[0:5])
print(ds[-6:-1])
