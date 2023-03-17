import os
import pandas as pd

for f in os.listdir('AverageInteralDegreeCheckResults/'):
    df = pd.read_csv('AverageInteralDegreeCheckResults/' + f, sep='\t')
    total = sum(df['average_internal_degree'] == 0)
    if total > 0:
        print(f,str(total))
        sub = df[df['average_internal_degree'] == 0]
        subsub = sub[sub['size'] > 1]
        print(subsub.shape[0])
        if subsub.shape[0] > 0:
            print(subsub)
