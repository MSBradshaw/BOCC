import pandas as pd
import datetime
import argparse

def get_args():
    parser = argparse.ArgumentParser()


    parser.add_argument('--year',
                        dest='year',
                        required=True,
                        help='year to include data up till january first of. So --year 2019 will gather all data from'
                             'before January 1 2019')

    parser.add_argument('--output',
                        dest='output',
                        required=True,
                        help='name of file to save the inferred g2p info at')

    args = parser.parse_args()
    return args

args = get_args()
"""
The purpose of this file to to infer the 2018 genes_to_phenotype.txt file from the current ones by using the date in
the Biocuration of the .hpoa file. This was what was recommended to me by Peter Robinson, a PI on the HPO project.
"""

# HPO versions hpo-web@1.7.13 - hpo-obo@2021-10-10
g2p = pd.read_csv('https://ci.monarchinitiative.org/view/hpo/job/hpo.annotations/lastSuccessfulBuild/artifact/rare-diseases/util/annotation/genes_to_phenotype.txt', sep='\t', comment='#', header=None)
#anno = pd.read_csv('http://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa', sep='\t', comment='#', header=None)
#anno = pd.read_csv('Resources/phenotype_annotation.tab', sep='\t', comment='#', header=None, dtype=str)
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
#anno.columns = ['DatabaseID', 'DiseaseName', 'Qualifier', 'HPO_ID', 'Reference', 'Evidence', 'Onset', 'Frequency',
#                'Sex', 'Modifier', 'Aspect', 'Biocuration']

#anno.columns = ["#disease-db", "disease-identifier", "disease-name", "negation", "HPO-ID", "reference", "evidence-code", "onset", "frequencyHPO", "modifier", "sub-ontology", "alt-names", "curators", "frequencyRaw", "sex"]
#disease-db	disease-identifier	disease-name	negation	HPO-ID	reference	evidence-code	onset	frequencyHPO	modifier	sub-ontology	alt-names	curators	frequencyRaw	sex
#print(anno['curators'])
# extract the date, if it doesn't have one make it a really big date that will be excluded
anno['date'] = [x.split('[')[1].split(']')[0] if '[' in str(x) else '2999-01-01' for x in anno['curators']]
mapping = {}
for i in range(anno.shape[0]):
    hp = str(anno.iloc[i, :]['HPO-ID'])
    if 'HP:' not in hp:
        print(anno.iloc[i, :])
        print()
    db_id = anno.iloc[i, :]['reference']
    date = anno.iloc[i, :]['date']
    date = int(datetime.datetime.strptime(date, '%Y-%m-%d').strftime("%s"))
    joint_id = hp + str(db_id)
    if joint_id not in mapping:
        mapping[joint_id] = date
    else:
        # is an entry is documented more than once, choose the earliest date
        if date < mapping[joint_id]:
            mapping[joint_id] = date

cut_off_time = int(datetime.datetime.strptime('{}-01-01'.format(args.year), '%Y-%m-%d').strftime("%s"))

# if the pair is in the annotation mapping, add the time stamp, otherwise add an arbitrarily HUGE date
g2p['date'] = [ mapping[g2p.iloc[i, 2] + g2p.iloc[i, 8]] if g2p.iloc[i, 2] + g2p.iloc[i, 8] in mapping else 9000000000 for i in range(g2p.shape[0])]

sub = g2p[g2p['date'] < cut_off_time]

sub.to_csv(args.output, sep='\t', index=False, header=False)
