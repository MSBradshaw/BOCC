import argparse
from pysam import TabixFile
import pysam
import typing
import networkx as nx
from BOCC.BOCC import load_clusters
from BOCC.BOCC import BOCC
import time
from multiprocessing import Pool
import os

# create a function that reads the VCF file that is bgzipped as a pytabix object
def read_vcf(vcf_path:str) -> pysam.VariantFile:
    return pysam.VariantFile(vcf_path)

# create a function to parse the command line arguments
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description='This is a simple command line interface that takes a VCF file, a genes file, a list of HPO terms,'
                                                    'a BOCC clusters file, an edgelist file, and an output directory.'
                                                    ' The code then finds the genes that are affected by the variants in the VCF file,'
                                                    ' finds the BOCC clusters that contain the affected genes and the HPO terms, and writes the matches  a file in the output directory.'
                                                    'The code also checks if the matches already exist in the edgelist file, and writes those matches to a separate file.')
    parser.add_argument('-v', '--vcf', required=True, help='Path to the VCF file')
    parser.add_argument('-g', '--genes', required=True, help='Path to a .bed file that lists the genes, the 4th column should contain the gene name')
    parser.add_argument('-p', '--hpos', required=True, help='Path to the file with HPO terms, one per line')
    parser.add_argument('-c', '--clusters', required=True, help='path to the file with the clusters file')
    parser.add_argument('-o', '--output', required=True, help='where to write the output dirrectory')
    parser.add_argument('-e', '--edgelist', required=True, help='Path to the edgelist file')
    # add a parameter for the number of processes to use
    parser.add_argument('-n', '--num_processes', required=False, default=1, help='Number of processes to use')
    return parser.parse_args()

# create a function that reads the genes file as a tabix object 
def read_genes(genes_path:str) -> TabixFile:
    return TabixFile(genes_path)
    
# create a function that find the intersecting genes for each record in the VCF file
def find_intersection(genes:TabixFile, vcf:pysam.VariantFile) -> typing.Set[str]:
    affected_genes = set()
    for record in vcf.fetch():
        for gene_string in genes.fetch(record.chrom, record.pos, record.pos+1):
            gene = gene_string.split('\t')
            affected_genes.add(gene[3])
    return affected_genes

# create a function that searchs for pairs of HPOs and affected genes in side the .members object of each BOCC cluster
def find_hpo_genes(clusters:typing.List[BOCC], affected_genes:set, hpos:list, G:nx.Graph) -> tuple[typing.List[str], typing.List[str]]:
    matches = []
    preexisting = []
    for com in clusters:
        for gene in affected_genes:
            for hpo in hpos:
                if hpo in com.members and gene in com.members:
                    # check if the hpo gene edge already exists in the graph
                    if G.has_edge(gene, hpo):
                        preexisting.append(f'{gene}\t{hpo}\t{com.name}')
                    else:
                        matches.append(f'{gene}\t{hpo}\t{com.name}')
    return matches, preexisting

# create a function to load the HPO file
def load_hpos(hpos_path:str) -> list:
    with open(hpos_path) as f:
        return [hpo.strip() for hpo in f]

# create a function to read the edgelist in as a networkx graph
def read_edgelist(edgelist_path:str) -> nx.Graph:
    return nx.read_edgelist(edgelist_path)

# create a function to write the matches and preexisting matches to seporate files in the output dirrectory
def write_matches(matches:typing.List[str], preexisting:typing.List[str], output_path:str) -> None:
    # if ouput path doesn end in a slash, add one
    if not output_path.endswith('/'):
        output_path += '/'
    # check if the output dirrectory exists, if not create it
    check_dir(output_path)
    # write the preexisting matches to a file
    with open(output_path + 'preexisting_matches.txt', 'w') as f:
        for match in preexisting:
            f.write(match + '\n')
    # write the preexisting matches to a file
    with open(output_path + 'novel_matches.txt', 'w') as f:
        for match in matches:
            f.write(match + '\n')

# create a another version of the find_hpo_genes function that search for each cluster indiviudally and is parallelized to use X cores
def find_hpo_genes_parallel(clusters:typing.List[BOCC], affected_genes:set, hpos:list, G:nx.Graph, cores:int) -> tuple[typing.List[str], typing.List[str]]:
    # create a pool of workers
    pool = Pool(cores)
    # create a list of arguments for each worker
    args = [([com], affected_genes, hpos, G) for com in clusters]
    # run the workers
    results = pool.starmap(find_hpo_genes, args)
    # close the pool
    pool.close()
    # create a list of matches and preexisting matches
    matches = []
    preexisting = []
    # add the results from each worker to the list of matches and preexisting matches
    for result in results:
        matches += result[0]
        preexisting += result[1]
    return matches, preexisting

# write a function to check if a dirrectory exists and if not create it
def check_dir(dir_path:str) -> None:
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)

def main():
    args = parse_args()
    vcf =   read_vcf(args.vcf)
    genes = read_genes(args.genes)
    affect_genes = find_intersection(genes, vcf)
    clusters = load_clusters(args.clusters)
    hpos = load_hpos(args.hpos)
    G = read_edgelist(args.edgelist)
    print('parallel starting')
    # measure the amount of time it takes to run the parallel version of the function
    start = time.time()
    matches2, preexisting2 = find_hpo_genes_parallel(clusters, affect_genes, hpos, G, args.num_processes)
    end = time.time()
    print(f'parallel took {end-start} seconds')
    print('parallel done')
    print('non-parallel starting')
    # measure the amount of time it takes to run the non-parallel version of the function
    start = time.time()
    matches1, preexisting1 = find_hpo_genes(clusters, affect_genes, hpos, G)
    end = time.time()
    print(f'non-parallel took {end-start} seconds')
    assert len(matches1) == len(matches2)
    assert len(preexisting1) == len(preexisting2)
    write_matches(matches2, preexisting2, args.output)

if __name__ == '__main__':
    main()
