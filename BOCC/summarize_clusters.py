from BOCC import load_clusters, summarize_clusters
import argparse
import networkx as nx

def get_args():
    parser = argparse.ArgumentParser()


    parser.add_argument('--coms',
                        dest='coms',
                        required=True,
                        help='tsv file where each row in a community, first area is com name followed by a tab ' +
                             'seporated list of community members (not all rows will be the same length)')

    parser.add_argument('--alpha',
                        dest='alpha',
                        type=float,
                        default=0.05,
                        help='threshold for significance in GO enrichment')

    parser.add_argument('--mg2',
                        dest='mg2',
                        type=str,
                        help='path to mygene2 family list file')

    parser.add_argument('--graph',
                        dest='graph',
                        type=str,
                        help='path to HPO-String edgelist to use')

    parser.add_argument('--out',
                        dest='out',
                        type=str,
                        required=True,
                        help='output file, should be a .tsv)')

    parser.add_argument('--ignore_api_requests',
                        dest='ignore_api_requests',
                        type=str,
                        default='False',
                        help='True or False, should api request stats be ignored?')

    parser.add_argument('--gnomad_file',
                        dest='gnomad_file',
                        type=str,
                        required=True,
                        help='gnomad all pLoF tsv file')
    # Data/gnomad.v2.1.1.all_lofs.txt

    args = parser.parse_args()
    return args

args = get_args()
coms = load_clusters(args.coms)
G = nx.read_edgelist(args.graph)
ignore_api_requests = args.ignore_api_requests == 'True'
print('First', str(ignore_api_requests))

print('args alpha', str(args.alpha))
df = summarize_clusters(clusters=coms, G=G, mygene2_file=args.mg2, gnomad_file=args.gnomad_file, alpha=args.alpha, ignore_api_requests=ignore_api_requests)
df.to_csv(args.out, sep='\t', index=False)
