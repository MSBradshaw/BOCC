import requests
import json
import pandas as pd
import typing
import matplotlib.pyplot as plt
import warnings
from statsmodels.stats.multitest import fdrcorrection
from cdlib import evaluation
import networkx as nx
from cdlib import NodeClustering
from statistics import median, stdev
from statistics import StatisticsError
import time

GO_ENRICHMENT_CACHE = {}

class BOCC:
    """
    Biological Ontology Cluster Comparison (BOCC) (Pronounced like Bach the musician)
    """

    def __init__(self):
        self.members = []
        self.types = []
        self.name = None
        self.genes = None
        self.go_results = None
        self.gene_plof_counts = None
        self.gene_plof_counts_list = None

    def reset(self) -> None:
        """
        Resets certain object values that are dependent on the content of the self.members list
        Anytime self.members is changed these values should be reset
        :return: None
        """
        self.genes = None
        self.go_results = None

    def add_members(self, mems: typing.List[str], types=None) -> None:
        """

        :param mems: List of cluster members to be appended
        :param types: List of node types associated with each of the items mems
        :return: None
        """
        self.reset()
        # if no types are given, create a list of equal length with mems to add
        if types is None:
            types = [None] * len(mems)
        self.members += mems
        self.types += types

    def get_genes(self) -> typing.List[str]:
        """
        Sets self.genes equal to a list of all the items in self.members that are listed as 'gene' in self.types
        If self.types contains only None then anything without the 'HP:' prefix is considered a gene
        :return: a list the genes in the self.members
        """
        if self.genes is not None:
            return self.genes
        # if there are no types listed assuming that anything without the HP: (denoting an HPO term) prefix is a gene
        elif self.types is None or sum(0 if x is None else 1 for x in self.types) == 0:
            self.genes = [x for x in self.members if 'HP:' not in x]
        else:
            self.genes = [self.members[i] for i in range(len(self.types)) if self.types[i] == 'gene']
        return self.genes

    def go_enrichment(self) -> pd.DataFrame:
        """
        Uses send the content of self.genes to the panther API for over representation analysis
        :return: pandas DataFrame of the panther API results, not all returned rows are significant
        """
        if self.go_results is not None:
            print('Returning Cached results')
            return self.go_results
        self.get_genes()

        gs = list(self.genes)
        gs.sort()
        # check for cached gene sets, this will speed up cesna (probably) 
        
        g_string = ','.join(self.genes)
        url = 'http://pantherdb.org/services/oai/pantherdb/enrich/overrep?geneInputList=' + g_string + \
              '&organism=9606&annotDataSet=GO%3A0008150&enrichmentTestType=FISHER&correction=FDR'

        if len(self.genes) > 5000:
            print('Too many genes for API request')
            self.go_results = pd.DataFrame()
            return pd.DataFrame()
        if len(self.genes) == 0:
            print('No genes for API request')                
            self.go_results = pd.DataFrame()
            return pd.DataFrame()
        print('Sending Request')
        # use 3 attempts to get the request before erroring out

        for i in range(3):
            print('request attempt ',str(i))
            try:
                resp = requests.get(url, timeout=60)
                break
            except requests.exceptions.RequestException as e:
                if i == 2:
                    print('Final attempt error')
                    raise(e)
                print('Time out error, sleeping for 3 minutes then re-attempting')
                time.sleep(180)

        print('Got request')
        resp_obj = json.loads(resp.content)

        results = {'number_in_list': [],
                   'fold_enrichment': [],
                   'fdr': [],
                   'expected': [],
                   'number_in_reference': [],
                   'pValue': [],
                   'id': [],
                   'label': [],
                   'plus_minus': []}
        try:
            for i in range(len(resp_obj['results']['result'])):
                for key in resp_obj['results']['result'][i].keys():
                    if key != 'term':
                        results[key].append(resp_obj['results']['result'][i][key])
                    else:
                        try:
                            results['id'].append(resp_obj['results']['result'][i]['term']['id'])
                        except KeyError:
                            results['id'].append('None')
                        try:
                            results['label'].append(resp_obj['results']['result'][i]['term']['label'])
                        except KeyError:
                            results['label'].append('None')

        except KeyError:
            # there are no results, return empty dataframe
            self.go_results = pd.DataFrame()
            return pd.DataFrame()

        self.go_results = pd.DataFrame(results)
        return self.go_results

    def get_gene_tissue_specificities(self, all_genes_dict: typing.Dict = None) -> \
            typing.Tuple[typing.Dict, typing.Dict, typing.Dict]:
        """
        Get the tissue specificity information for each gene in the cluster. This is the information taken from the
        Human Protein Atlas (https://www.proteinatlas.org/).
        The first dictionary contains 'RNA tissue specificity' with is a verbal measurement of how strong the
        specificity is.
        The second dictionary contains the info from 'RNA single cell type specificity' type of enrichment that occurs
        to the group.
        The third dictionary contains the info from 'RNA single cell type specific NX' with the cell-type / tissue
        the gene is specific for.
        :param all_genes_dict: dictionary of human protein atlas information. This will be auto loaded is not given
        directly, but giving the function a preloaded version makes this process quicker, especially if you are
        doing thousands of communities.
        :return: tuple of two dictionaries
        """

        self.get_genes()
        gene_ts = {}
        gene_sc_type = {}
        gene_sc_type_info = {}
        if all_genes_dict is None:
            with open('BOCC/all_genes_info.json', 'r') as f:
                all_genes_dict = json.load(f)
        for g in self.genes:
            try:
                gene_ts[g] = all_genes_dict[g]['RNA tissue specificity']
                gene_sc_type[g] = all_genes_dict[g]['RNA single cell type specificity']
                gene_sc_type_info[g] = all_genes_dict[g]['RNA single cell type specific NX']
            except KeyError:
                warnings.warn('Warning, gene not found in HPA ' + g)
                gene_ts[g] = None
                gene_sc_type[g] = None
                gene_sc_type_info[g] = None

        return gene_ts, gene_sc_type, gene_sc_type_info

    def get_diseases_associated_with_genes(self, all_genes_dict: typing.Dict = None) -> typing.Dict:
        """
        Get the disease associated with each gene in the cluster. This is the information taken from the
        Human Protein Atlas (https://www.proteinatlas.org/).
        The first dictionary contains 'RNA tissue specificity' with is a verbal measurement of how strong the
        specificity is.
        :param all_genes_dict: dictionary of human protein atlas information. This will be auto loaded is not given
        directly, but giving the function a preloaded version makes this process quicker, especially if you are
        doing thousands of communities.
        :return: tuple of two dictionaries
        """

        self.get_genes()

        gene_diseases = {}
        if all_genes_dict is None:
            with open('BOCC/all_genes_info.json', 'r') as f:
                all_genes_dict = json.load(f)
        for g in self.genes:
            try:
                gene_diseases[g] = all_genes_dict[g]['Disease involvement']
            except KeyError:
                # control for when there are no disease because the gene is not found in HPA
                gene_diseases[g] = None
        return gene_diseases

    def get_disease_counts(self, all_genes_dict: typing.Dict = None) -> typing.Dict:
        """
        Count the number of times each diseases is associated with a gene in the cluster
        :param all_genes_dict: dictionary of human protein atlas information. This will be auto loaded is not given
        directly, but giving the function a preloaded version makes this process quicker, especially if you are
        doing thousands of communities.
        :return: dictionary keys are disease and values are occurrence counts (number of times it appeared in the community)
        """
        gds = self.get_diseases_associated_with_genes(all_genes_dict)
        disease_counts = {}
        for key in gds.keys():
            diseases = gds[key]
            # control for when there are no disease because the gene is not found in HPA
            if diseases is None:
                continue
            for d in diseases:
                if d in disease_counts:
                    disease_counts[d] += 1
                else:
                    disease_counts[d] = 1

        return disease_counts

    def summarize_disease_associations(self, all_genes_dict: typing.Dict = None):
        """
        Count the number of times each diseases is associated with a gene in the cluster
        :param all_genes_dict: dictionary of human protein atlas information. This will be auto loaded is not given
        directly, but giving the function a preloaded version makes this process quicker, especially if you are
        doing thousands of communities.
        :return: the most common disease association
        """
        disease_counts = self.get_disease_counts(all_genes_dict)
        return get_max_in_dict(disease_counts)

    def get_number_of_diseases(self, all_genes_dict: typing.Dict = None):
        """
        Count the number of times each diseases is associated with a gene in the cluster
        :param all_genes_dict: dictionary of human protein atlas information. This will be auto loaded is not given
        directly, but giving the function a preloaded version makes this process quicker, especially if you are
        doing thousands of communities.
        :return: the most common disease association
        """
        disease_counts = self.get_disease_counts(all_genes_dict)
        return len(disease_counts)

    def get_cell_type_counts(self, all_genes_dict: typing.Dict = None) -> typing.Dict:
        """
        count up the number of times each cell type is listed as being specific for any of the genes in the community
        :param all_genes_dict: dictionary of human protein atlas information. This will be auto loaded is not given
        directly, but giving the function a preloaded version makes this process quicker, especially if you are
        doing thousands of communities.
        :return: dictionary keys are disease and values are occurrence counts (number of times it appeared in the community)
        """
        dicts = self.get_gene_tissue_specificities(all_genes_dict)
        # the third dict keys are the thing of interest here
        cell_type_counts = {}
        for gene in dicts[2].keys():
            if dicts[2][gene] is None:
                continue
            for cell_type in dicts[2][gene].keys():
                if cell_type in cell_type_counts:
                    cell_type_counts[cell_type] += 1
                else:
                    cell_type_counts[cell_type] = 1

        return cell_type_counts

    def summarize_cell_type_specificity(self, all_genes_dict: typing.Dict = None):
        """
        Get the most frequent cell type
        :param all_genes_dict: dictionary of human protein atlas information. This will be auto loaded is not given
        directly, but giving the function a preloaded version makes this process quicker, especially if you are
        doing thousands of communities.
        :return:
        """
        cell_type_counts = self.get_cell_type_counts(all_genes_dict)
        return get_max_in_dict(cell_type_counts)

    def get_summary_stats(self, mygene2_file: str, G: nx.Graph,  gnomad_file:str, alpha: float = 0.05,
                          ignore_api_requests=True) -> pd.DataFrame:
        """
        This function returns a pd.DataFrame with information summarizing the biological relevance of a cluster
        :return:
        """
        res = {'cluster_id': [], 'cluster_size': [], 'gene_ratio': [], 'HPO_ratio': [],
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
        # cluster id/name
        res['cluster_id'].append(self.name)
        # number size of cluster
        res['cluster_size'].append(len(self.members))
        # gene ratio
        res['gene_ratio'].append(len(self.get_genes()) / len(self.members))
        # HPO ratio
        res['HPO_ratio'].append((len(self.members) - len(self.get_genes())) / len(self.members))
        # go enrichment, # terms with p-value > p_tresh
        print('Ignore API', str(ignore_api_requests))
        if ignore_api_requests:
            go_df = pd.DataFrame()
        else:
            go_df = self.go_enrichment()

        if go_df.shape[0] == 0:
            print('No results from GO Enrichment to summarize')
            # return pd.DataFrame()

        if go_df.shape[0] == 0:
            res['num_sig_go_enrichment_terms'].append(-1)
            # go enrichment, comma separated string of terms with p-value > p_tresh
            res['sig_go_enrichment_p_vals'].append(-1)
            res['sig_go_enrichment_fdr_corrected_p_vals'].append(-1)
            res['sig_go_enrichment_terms'].append(-1)
        else:
            # FDR correction
            print('ALPHA!!!!!!!!!')
            print(alpha)
            print()
            print()
            print(list(go_df['pValue']))
            print()
            print()
            rejected_null, pvals = fdrcorrection(list(go_df['pValue']), alpha=alpha)
            # add FDR correction to the results df
            go_df['significant'] = rejected_null
            go_df['fdr_pVale'] = pvals

            go_df = go_df[go_df['significant'] == True]
            res['num_sig_go_enrichment_terms'].append(go_df.shape[0])
            # go enrichment, comma separated string of terms with p-value > p_tresh
            res['sig_go_enrichment_p_vals'].append(','.join([str(x) for x in go_df['pValue']]))
            res['sig_go_enrichment_fdr_corrected_p_vals'].append(','.join([str(x) for x in go_df['fdr_pVale']]))
            res['sig_go_enrichment_terms'].append(','.join([str(x) for x in go_df['id']]))
        res['go_sig_threshold'].append(alpha)
        # max normalized cell type specificity value
        cts = self.summarize_cell_type_specificity()
        try:
            res['max_norm_cell_type_specificity'].append(cts[1] / len(self.genes))
        except ZeroDivisionError:
            res['max_norm_cell_type_specificity'].append(-1)
        # comma separated list of cells with max normalized cell type specificity
        res['max_norm_cell_type_comma_sep_string'].append(','.join(cts[0]))
        # max normalized disease specificity value
        ds = self.summarize_disease_associations()
        try:
            res['max_norm_disease_specificity'].append(ds[1] / len(self.genes))
        except ZeroDivisionError:
            res['max_norm_disease_specificity'].append(-1)
        nd = self.get_number_of_diseases()
        res['num_of_diseases'].append(nd)
        # comma separated list of cells with max normalized disease specificity
        # the disease(s) with the most number of genes associated with them in a cluster
        # (if there is a tie all the highest count diseases are listed and separated by a comma)
        if len(ds[0]) == 0:
            res['max_norm_disease_comma_sep_string'].append(','.join(['No Associated Disease']))
        else:
            res['max_norm_disease_comma_sep_string'].append(','.join(ds[0]))
        if mygene2_file is not None and False:
            pc, npc, portion_recovered = self.mygene2_stats(mygene2_file)
        else:
            pc, npc, portion_recovered = 0, 0, 0
        res['mg2_pairs_count'].append(pc)
        res['mg2_not_pairs_count'].append(npc)
        res['mg2_portion_families_recovered'].append(portion_recovered)

        # use a bunch of the cdlib evaluation methods, originally intended for a group of communities,
        # reused here for individual communities
        communities = NodeClustering([self.members], graph=G, method_name="nada_nada_limonada")
        ave = evaluation.avg_embeddedness(G, communities)
        res['avg_embeddedness'].append(ave.score)
        mod = evaluation.average_internal_degree(G, communities)
        res['avg_internal_degree'].append(mod.score)
        cond = communities.conductance()
        res['conductance'].append(cond.score)
        cr = evaluation.cut_ratio(G, communities)
        res['cut_ratio'].append(cr.score)
        nc = evaluation.normalized_cut(G, communities)
        res['normalized_cut'].append(nc.score)
        ex = evaluation.expansion(G, communities)
        res['expansion'].append(ex.score)
        tpr = evaluation.triangle_participation_ratio(G, communities)
        res['triangle_participation_ratio'].append(tpr.score)
        spr = evaluation.surprise(G, communities)
        res['surprise'].append(spr.score)
        signif = evaluation.significance(G, communities)
        res['significance'].append(signif.score)
        ngm = evaluation.newman_girvan_modularity(G, communities)
        res['newman_girvan_modularity'].append(ngm.score)
        res['internal_edge_density'].append(communities.internal_edge_density().score)
        res['edges_inside'].append(communities.edges_inside().score)
        res['max_plof'].append(self.get_max_plof(gnomad_file))
        res['mean_plof'].append(self.get_mean_plof(gnomad_file))
        res['median_plof'].append(self.get_median_plof(gnomad_file))
        res['std_plof'].append(self.get_stdev_plof(gnomad_file))
        res['sum_plof'].append(self.get_sum_plof(gnomad_file))
        try:
            res['hub_dominance'].append(communities.hub_dominance().score)
        except ZeroDivisionError:
            warnings.warn('Waning ZeroDivisionError in hub_dominance, adding -1 instead of value')
            res['hub_dominance'].append(-1)
        return pd.DataFrame(res)

    def mygene2_stats(self, mygene2_file: str):
        """

        :param mygene2_file: path to tab separated file with columns: gene, HPO, familyID
        :return:
        """
        pairs_count = 0
        no_pairs_count = 0
        family_pairs_no_pairs = {}
        for line in open(mygene2_file, 'r'):
            row = line.strip().split()
            gene = row[0]
            hpo = row[1]
            family = row[2]
            if family not in family_pairs_no_pairs:
                family_pairs_no_pairs[family] = {'pairs': [], 'not_pairs': []}
            if gene in self.members and hpo in self.members:
                family_pairs_no_pairs[family]['pairs'].append((gene, hpo))
                pairs_count += 1
            else:
                family_pairs_no_pairs[family]['not_pairs'].append((gene, hpo))
                no_pairs_count += 1
        # total number of apirs
        # count number of family with > 50% pairs
        majority_pair_fams = []
        total_num_fams = 0
        for fam in family_pairs_no_pairs.keys():
            total_num_fams += 1
            if len(family_pairs_no_pairs[fam]['pairs']) > len(family_pairs_no_pairs[fam]['not_pairs']):
                majority_pair_fams.append(fam)
        return pairs_count, no_pairs_count, len(majority_pair_fams) / total_num_fams, family_pairs_no_pairs


    def get_num_new_edges(self, _new_edges: typing.List) -> int:
        _com = set(self.members)
        count = 0
        for _e in _new_edges:
            if _e[0] in _com and _e[1] in _com:
                count += 1
        return count

    def get_new_edges(self, _new_edges: typing.List) -> typing.List:
        """

        :param _new_edges:
        :return: list of strings of the new edges HPO listed first followed by the Gene
        """
        _com = set(self.members)
        count = 0
        _new_ones = list()
        for _e in _new_edges:
            if _e[0] in _com and _e[1] in _com:
                if 'HP:' in _e[0]:
                    _new_ones.append(str([_e[0], _e[1]]))
                else:
                    _new_ones.append(str([_e[1], _e[0]]))
        return _new_ones

    def collect_plof_stats_genes(self, gnomad_file: str) -> None:
        if self.gene_plof_counts is not None:
            return
        g_df = pd.read_csv(gnomad_file, sep='\t')
        self.gene_plof_counts = {}
        self.gene_plof_counts_list = []
        # for each gene, count the number of times it appears in the gnomad pLoF file
        # (the number of pLoF causing variants it has in that file)
        for g in self.get_genes():
            symbol_count = list(g_df['gene_symbols']).count(g)
            id_count = list(g_df['gene_ids']).count(g)
            # the gene is either listed as a symbol or an id, one will be zero one might not be
            self.gene_plof_counts[g] = max([symbol_count, id_count])
            self.gene_plof_counts_list.append(max([symbol_count, id_count]))

    def get_max_plof(self, gnomad_file: str):
        self.collect_plof_stats_genes(gnomad_file)
        return max(self.gene_plof_counts_list)

    def get_mean_plof(self, gnomad_file: str):
        self.collect_plof_stats_genes(gnomad_file)
        return sum(self.gene_plof_counts_list) / len(self.gene_plof_counts_list)

    def get_median_plof(self, gnomad_file: str):
        self.collect_plof_stats_genes(gnomad_file)
        return median(self.gene_plof_counts_list)

    def get_stdev_plof(self, gnomad_file: str):
        self.collect_plof_stats_genes(gnomad_file)
        try:
            return stdev(self.gene_plof_counts_list)
        except StatisticsError:
            print('StatisticsError only one value, cannot calc std')
            return -1

    def get_sum_plof(self, gnomad_file: str):
        self.collect_plof_stats_genes(gnomad_file)
        return sum(self.gene_plof_counts_list)


def get_max_in_dict(d):
    max_keys = []
    max_key_value = 0
    for key in d.keys():
        # TODO there are some disease that are not diseases like 'Disease mutation' that should be ignored
        diseases_to_ignore = ['Disease mutation']
        if key in diseases_to_ignore:
            continue

        if d[key] > max_key_value:
            max_key_value = d[key]
            max_keys = [key]
        elif d[key] == max_key_value:
            max_keys.append(key)

    return max_keys, max_key_value


def load_clusters(input_file: str) -> typing.List[BOCC]:
    """
    Load a community file into a list of BOCC objects
    :param input_file: path to file with each community listed on one tab separated line with the first item being the
    community name/ID
    :return: list of BOCC objects
    """
    clusters = []
    for line in open(input_file, 'r'):
        # check if line is commented out
        if line[0] == '#':
            continue
        row = line.strip().split('\t')
        clusters.append(BOCC())
        clusters[-1].add_members(row[1:])
        clusters[-1].name = row[0]
    return clusters


def plot_basic_com_stats(coms: typing.List[BOCC], output: str = None, logx: bool = False) -> typing.Dict[
    str, typing.List]:
    """
    Take a list of communities and plot the size vs the HPO:Gene ratio of each community
    :param coms: list of BOCC objects to plot
    :param output: where to save figure, if None fig is not created or saved
    :param logx: if the x scale should be log
    :return: dictionary with keys 'ratio' and 'size' all values are lists
    """
    results = {'ratio': [], 'size': []}
    for com in coms:
        gl = len(com.get_genes())
        tl = len(com.members)
        # using pseudo counts to avoid a dividing by zero issue
        hpo_gene_ratio = ((tl + 1) - (gl + 1)) / float((gl + 1))
        results['ratio'].append(hpo_gene_ratio)
        results['size'].append(tl)
    if output is not None:
        plt.scatter(results['size'], results['ratio'])
        plt.xlabel('Size')
        if logx:
            plt.xscale('log')
        plt.ylabel('Ratio (HPO:Gene)')
        plt.savefig(output)
    return results


def summarize_clusters(clusters: typing.List[BOCC], G: nx.Graph, mygene2_file: str, gnomad_file: str, alpha: float,
                       ignore_api_requests=True) -> pd.DataFrame:
    """
    :param ignore_api_requests:
    :param gnomad_file:
    :param alpha:
    :param mygene2_file:
    :param G:
    :param clusters: list of BOCC objects
    :param p_tresh: threshold to be used for go term significance
    :return: pd.DataFrame of biologically relevant information for each cluster
    """
    print('BOCC/BOCC.py summarize_clusters')
    df = None
    for c in clusters:
        print('Cluster ID',str(c.name))
        if len(c.members) < 3:
            continue
        if len(c.get_genes()) == len(c.members) or len(c.get_genes())  == 0:
            continue
        print('Prior', str(ignore_api_requests))
        d = c.get_summary_stats(mygene2_file=mygene2_file, G=G, alpha=alpha, gnomad_file=gnomad_file, ignore_api_requests=ignore_api_requests)
        if df is None:
            df = d
        else:
            df = pd.concat([df, d])
    return df
