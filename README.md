# BOCC

Biological Ontology Cluster Comparison (BOCC) - Leverages the knowledge contained with in resources like the STRING protein-ptotein-interaction networks and the Human Phenotype Ontology to hypothesis about latent gene-to-phenotype connections. The CLI and Web App below provides access to clusters of genes and phenotypes where currently unconnected genes are phenotype are much more likely to be connected in the future.  

## Visualization Web App

Use the [BOCC web app](https://ryanlayerlab.github.io/BOCC/) to explore and visualize the clusters in a fun interactive way!

##  Command Line Interface (CLI)
This CLI is useful for integration into automated pipeline and or for searching for large numbers of genes and phenoype pairs quickly.

### Set-up

Download the CLI script

`wget https://raw.githubusercontent.com/MSBradshaw/BOCC/cli/bocc_cli.py`

Download the edge list

`wget --no-check-certificate 'https://docs.google.com/uc?export=download&id=1BDrUa4Y55fFjx3j2tz7CzwD4qmgZS3fp' -O edgelist.2022.txt`

Download genes file (HG37)

```
wget --no-check-certificate 'https://docs.google.com/uc?export=download&id=17EHfiA3aiu2LouwhpNgX5VhC6ZaMcEpk' -O hg37.exons.bed.gz
wget --no-check-certificate 'https://docs.google.com/uc?export=download&id=1jxCuK7xKzz5r2cFlVkr2nFAfeW5c1EsG' -O hg37.exons.bed.gz.tbi
```

Download clusters

`wget --no-check-certificate 'https://docs.google.com/uc?export=download&id=1ZfjtKn3PxB8tTe0kg-OxH5BJlvcU30fz' -O all_clusters.2022.txt `

Download example data

```
wget https://raw.githubusercontent.com/MSBradshaw/BOCC/cli/Example/test.hpo.txt
wget https://github.com/MSBradshaw/BOCC/raw/cli/Example/test.vcf.gz
wget https://github.com/MSBradshaw/BOCC/raw/cli/Example/test.vcf.gz.tbi
```

Configure environment

`conda env create -f  env.yml`

### Usage
```
conda activate bocc

python bocc_cli.py -v test.vcf.gz \
-g hg37.exons.bed.gz \
-p test.hpo.txt \
-c all_clusters.2022.txt \
-o TestDir \
-e edgelist.2022.txt \
-n 1
```

### Output

BOCC will create the dirrectory specified with `-o` and place in it two files: 

`novel_matches.txt` - which contains a list of affected genes and HPO terms that co-occur in the same cluster and are currently not connected - suggestive of a talent gene-to-phenotype connection.

`preexisting_matches.txt` - which contains any gene-to-phenotype pairs from the input data that already have a known connection.

### Arguments
```
-v, --vcf: Path to the tabixed and BGZipped VCF file. Alternatively you an also provide a tabixed and BGZipped .bed file.
-g, --genes: Path to a .bed file that lists the genes, the 4th column should contain the gene name
-p, --hpos: Path to the file containing the HPO terms - one per line
-c, --clusters: Path to the BOCC cluster file
-o, --output: Path to the output dirrectory
-e, --edgelist: Path to the edgelist
-n, --num_processes: Number of cores to use (default 1). If searching for <= 50 genes and phenotypes the non-parallelized version is quicker.
--verbose: flag to enable verbose mode
```


## Reproduce Results

We are committed to open and reproducible science. To best facilitate the reproducibiliy of our results we have bundled all of our methods into a snakemake pipeline and accompanying conda enviroment. Fair warning just because it is reproducible does not mean it is quick. The graph we repretendly cluster is large and some of the algorithms are slow and can take a week or more to complete even with 64 cores running in parallel. Generating the results depends on the Panther API for enrichment analysis which can also be an extremely slow processes and is prone to blocking IP addresses if you are try making more than 20 requests at the same time, FYI.

### Conda Env

`conda env create -f  reproduction.env.yml`

### Reproduce Results

The figures and tables seen in the article can be reproduced with the jupyter notebooks in `Experiments`.

### Reproduce data & clusters

```
conda activate reproduce_bocc
snakemake -s run_everything.smk --cores 64
```

This pipeline will produce the clusters them selves and .tsv files containing the features used by XGBoost for usefulness predictions. The training, testing and use of XGBoost reproduction can be found [here](https://github.com/ConGibbs10/BOCCRank).
