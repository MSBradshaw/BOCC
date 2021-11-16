import pandas as pd
import datetime

"""
The purpose of this file to to infer the 2018 genes_to_phenotype.txt file from the current ones by using the date in
the Biocuration of the .hpoa file. This was what was recommended to me by Peter Robinson, a PI on the HPO project.
"""

# HPO versions hpo-web@1.7.13 - hpo-obo@2021-10-10
g2p = pd.read_csv('https://ci.monarchinitiative.org/view/hpo/job/hpo.annotations/lastSuccessfulBuild/artifact/rare-diseases/util/annotation/genes_to_phenotype.txt', sep='\t', comment='#', header=None)
anno = pd.read_csv('http://purl.obolibrary.org/obo/hp/hpoa/phenotype.hpoa', sep='\t', comment='#', header=None)

anno.columns = ['DatabaseID', 'DiseaseName', 'Qualifier', 'HPO_ID', 'Reference', 'Evidence', 'Onset', 'Frequency',
                'Sex', 'Modifier', 'Aspect', 'Biocuration']

anno['date'] = [x.split('[')[1].split(']')[0] for x in anno['Biocuration']]

mapping = {}
for i in range(anno.shape[0]):
    hp = anno.iloc[i, :]['HPO_ID']
    db_id = anno.iloc[i, :]['DatabaseID']
    date = anno.iloc[i, :]['date']
    date = int(datetime.datetime.strptime(date, '%Y-%m-%d').strftime("%s"))
    if hp + db_id not in mapping:
        mapping[hp + db_id] = date
    else:
        # is an entry is documented more than once, choose the earliest date
        if date < mapping[hp + db_id]:
            mapping[hp + db_id] = date

cut_off_time = int(datetime.datetime.strptime('2019-01-01', '%Y-%m-%d').strftime("%s"))

# if the pair is in the annotation mapping, add the time stamp, otherwise add an arbitrarily HUGE date
g2p['date'] = [ mapping[g2p.iloc[i, 2] + g2p.iloc[i, 8]] if g2p.iloc[i, 2] + g2p.iloc[i, 8] in mapping else 9000000000 for i in range(g2p.shape[0])]

sub = g2p[g2p['date'] < cut_off_time]

sub.to_csv('genes_to_phenotype.txt', sep='\t', index=False, header=False)