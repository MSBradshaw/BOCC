from BOCC import load_clusters, summarize_clusters
import argparse
import networkx as nx
import pandas as pd

res_df_template = {'cluster_id': [], 'cluster_size': [], 'gene_ratio': [], 'HPO_ratio': [],
    'num_sig_go_enrichment_terms': [],
    'sig_go_enrichment_p_vals': [], 'sig_go_enrichment_fdr_corrected_p_vals': [],
    'sig_go_enrichment_terms': [], 'go_sig_threshold': [], 'max_norm_cell_type_specificity': [],
    'max_norm_cell_type_comma_sep_string': [], 'num_of_diseases': [], 'max_norm_disease_specificity': [],
    'max_norm_disease_comma_sep_string': [], 'mg2_pairs_count': [], 'mg2_not_pairs_count': [],
    'mg2_portion_families_recovered': [], 'avg_embeddedness': [], 'avg_internal_degree': [],
    'conductance': [], 'cut_ratio': [], 'normalized_cut': [], 'expansion': [],
    'triangle_participation_ratio': [], 'surprise': [], 'significance': [],
    'newman_girvan_modularity': [], 'internal_edge_density': [], 'edges_inside': [],
    'hub_dominance': [], 'max_plof': [], 'mean_plof': [], 'median_plof': [], 'std_plof': [], 'sum_plof': []}

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

def is_trivial(_com):
    _g = _com.get_genes()
    if len(_com.members) < 3:
        return True
    elif len(_g) == 0 or len(_g) == len(_com.members):
        return True
    else:
        return False

args = get_args()
coms = load_clusters(args.coms)
G = nx.read_edgelist(args.graph)
ignore_api_requests = args.ignore_api_requests == 'True'
print('First', str(ignore_api_requests))

print('args alpha', str(args.alpha))
df = summarize_clusters(clusters=coms, G=G, mygene2_file=args.mg2, gnomad_file=args.gnomad_file, alpha=args.alpha, ignore_api_requests=ignore_api_requests)
print('Coms:')
for com in coms:
    print(com.members)
print('--End coms--')

if df is None:
    print('Warning: df is None. Are all clusters trivial?')
    # check if all coms are trivial (< 3 members, not heterogenious)
    assert sum([is_trivial(com) for com in coms]) == len(coms)
    # if they are all indeed trivial, set df equal to a blank template df
    df = pd.DataFrame(res_df_template)
df.to_csv(args.out, sep='\t', index=False)
