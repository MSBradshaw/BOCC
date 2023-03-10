import sys
import pandas as pd

df = pd.read_csv(sys.argv[1],sep='\t')
# incse there are extra headers floating around
df = df[df['cluster_id'] != 'cluster_id']
if len(set(df['cluster_id'])) != df.shape[0]:
    print(sys.argv)
