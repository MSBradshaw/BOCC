YEARS = ['2019', '2020', '2021','2022']
YEARS_short = ['2019', '2020','2021']
ALGOS = ['greedy','walktrap','infomap','cesna']
HIER_ALGOS = ['paris','ward','louvain']
HIER_ALGOS = ['paris']
rule all:
	input:
		'work/start.start',
		'BOCC/all_genes_info.json',
		expand('work/{year}/hp.obo',year=YEARS),
		expand('work/{year}/genes_to_phenotype.txt',year=YEARS),
		expand('work/{year}/string_edge_list.txt',year=YEARS),
		'work/Resources/9606.protein.info.v11.5.txt',
		expand('Edgelists/String_HPO_{year}.phenotypic_branch.edgelist.txt',year=YEARS),
		expand('Edgelists/String_HPO_{year}.phenotypic_branch.numbered.edgelist.txt',year=YEARS),
		expand('Edgelists/String_HPO_{year}.phenotypic_branch.nodenames.txt',year=YEARS),
		expand('Edgelists/String_HPO_{year}.phenotypic_branch.nodeattributes.txt',year=YEARS),
		expand('Edgelists/String_HPO_{year}.all_hpo.edgelist.txt',year=YEARS),
		expand('Edgelists/String_HPO_{year}.all_hpo.numbered.edgelist.txt',year=YEARS),
		expand('Edgelists/String_HPO_{year}.all_hpo.nodenames.txt',year=YEARS),
		expand('Edgelists/String_HPO_{year}.all_hpo.nodeattributes.txt',year=YEARS),
		expand('work/{year}/inferred_genes_to_phenotype.txt',year=YEARS),
		expand('Clusters/{year}/greedy.{year}.coms.txt',year=YEARS),
		expand('Clusters/{year}/walktrap.{year}.coms.txt',year=YEARS),
		expand('g2p_Edgelists/String_HPO_{year}.phenotypic_branch.g2p_edgelist.txt',year=YEARS),
		expand('Clusters/{year}/cesna.{year}.coms.txt',year=YEARS),
		'Algorithms/cesna',
		expand('Clusters/{year}/infomap.{year}.coms.txt',year=YEARS),
		expand('SubComs/{year}/paris.{algo}.{year}.coms.txt',algo=ALGOS,year=YEARS),
		expand('SubComs/{year}/ward.{algo}.{year}.coms.txt',algo=ALGOS,year=YEARS),
		expand('SubComs/{year}/louvain.{algo}.{year}.coms.txt',algo=ALGOS,year=YEARS),
		#expand('ResultsBOCC/{year}/paris.{algo}.{year}.bocc_res.tsv',algo=ALGOS,year=YEARS),
		expand('ResultsBOCC/{year}/paris.{algo}.{year}.bocc_res.tsv',algo=ALGOS,year=['2022']),
		expand('DrugRediscovery/snowball.{hier}.{algo}.String_HPO_2022.phenotypic_branch.tsv',algo=ALGOS,hier=HIER_ALGOS),
		#expand('ResultsBOCC/{year}/{hier}.{algo}.{year}.bocc_res.tsv',algo=ALGOS,year=YEARS,hier=HIER_ALGOS),
		#expand("SnowballResults/snowball.{hier}.{algo}.String_HPO_{year}.phenotypic_branch.tsv",algo=ALGOS,year=YEARS_short,hier=HIER_ALGOS),
		#expand('CleanedSplit100BOCCResultsCombined/{hier}.{algo}.{year}.bocc_res.tsv',algo=ALGOS,year=YEARS,hier=HIER_ALGOS),
		expand('SnowballedCleanedSplit100BOCCResultsCombinedFixed/{year}/{hier}.{algo}.{year}.bocc_res.tsv',algo=ALGOS,year=YEARS,hier=HIER_ALGOS),
		expand('AnyNewEdgeCountSnowballedCleanedSplit100BOCCResultsCombinedFixed/{year}/{hier}.{algo}.{year}.bocc_res.tsv',algo=ALGOS,year=YEARS_short,hier=HIER_ALGOS),
		#expand('AverageInteralDegreeCheckResults/{algo}.{year}.coms.txt',year=YEARS,algo=ALGOS),
		#expand('MouseSnowballResults/snowball.{hier}.{algo}.String_HPO_2021.mouse.phenotypic_branch.tsv', hier=HIER_ALGOS, algo=ALGOS),
		expand('SnowballResultsFixed/snowball.paris.{algo}.String_HPO_{year}.phenotypic_branch.tsv',algo=ALGOS,year=YEARS_short),
		expand('DrugRediscoveryResults/snowball.results.{hier}.{algo}.String_HPO_2022.phenotypic_branch.tsv',algo=ALGOS,hier=HIER_ALGOS)

"""
--- Collect The Data ---
"""

rule download_hpo:
	input:
		'work/start.start'
	output:
		hp15='work/2015/hp.obo',
		hp16='work/2016/hp.obo',
		hp17='work/2017/hp.obo',
		hp18='work/2018/hp.obo',
		hp19='work/2019/hp.obo',
		hp20='work/2020/hp.obo',
		hp21='work/2021/hp.obo',
		hp22='work/2022/hp.obo'
	shell:
		"""
		mkdir -p work/2017
		wget --no-check-certificate https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/9a977e2b9d68eb964deef586c7645531733ae5f7/hp.obo
		mv hp.obo {output.hp17}

		mkdir -p work/2021
		wget --no-check-certificate https://github.com/obophenotype/human-phenotype-ontology/raw/5ebbeb47e051f849241e9f055ce1afc91f1c906c/hp.obo
		mv hp.obo {output.hp21}
		
		mkdir -p work/2015
		wget --no-check-certificate https://github.com/drseb/HPO-archive/raw/master/2014-2015/2015_week_46/hpo/artefacts/hp.obo.gz
		gunzip hp.obo.gz
		mv hp.obo {output.hp15}
		
		mkdir -p work/2016
		wget --no-check-certificate https://github.com/drseb/HPO-archive/raw/master/2016-2017/hp/1702/archive.zip 
		unzip archive.zip
		mv archive/hp/hp.obo {output.hp16}
		rm -rf archive*

		mkdir -p work/2018
		wget --no-check-certificate https://github.com/obophenotype/human-phenotype-ontology/raw/b745ba91f22d641629997aa4de512ac32a1988ea/hp.obo
		mv hp.obo {output.hp18}

		mkdir -p work/2019
		wget --no-check-certificate https://github.com/obophenotype/human-phenotype-ontology/raw/13f9d1cf2f1a33473996a72f13b87420c4a8ad95/hp.obo
		mv hp.obo {output.hp19}
		
		mkdir -p work/2020
		wget --no-check-certificate https://github.com/obophenotype/human-phenotype-ontology/raw/15b018a427667e33bd6e51c136414c17c952795d/hp.obo
		mv hp.obo {output.hp20}

		mkdir -p work/2022
		wget --no-check-certificate https://raw.githubusercontent.com/obophenotype/human-phenotype-ontology/65395ade3a3cb1bec2929242cc358484730172c3/hp.obo
                mv hp.obo {output.hp22}
		"""

rule download_jenkins:
	input:
		'work/start.start'
	output:
		j17='work/2017/genes_to_phenotype.txt',
		j15='work/2015/genes_to_phenotype.txt',
		j21='work/2021/genes_to_phenotype.txt',
		j16='work/2016/genes_to_phenotype.txt',
		j19='work/2019/genes_to_phenotype.txt',
		j20='work/2020/genes_to_phenotype.txt',
		j18='work/2018/genes_to_phenotype.txt',
		j22='work/2022/genes_to_phenotype.txt'
	log:
		'logs/download_jenkins.txt'
	shell:
		"""
		for x in $(seq 2015 2021);
		do
			mkdir -p work/$x
		done
		#This is the data from Feb 2017, we know this because the sql dump file in this archive is labeled MYHPO_02_2017.sql
		wget --no-check-certificate https://github.com/drseb/HPO-archive/raw/master/hpo.annotations.monthly/117/archive.zip
		unzip archive.zip 2>> {log}
		cp archive/annotation/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt {output.j17}  2>> {log}
		rm -rf archive*  2>> {log}

		# get the December 2015 data
		wget --no-check-certificate https://github.com/drseb/HPO-archive/raw/master/hpo.annotations.monthly/2015-12-01_00-00-05/archive/annotation/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt.gz 2>> {log}
		gunzip ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt.gz 2>> {log}
		mv ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt {output.j15} 2>> {log}

		# get the 2021 data
		# This URL is no longer accessible, this change was discovered by Michael on Dec 6 2022
		wget --no-check-certificate https://ci.monarchinitiative.org/view/hpo/job/hpo.annotations/11/artifact/rare-diseases/util/annotation/genes_to_phenotype.txt 2>> {log}
		cat genes_to_phenotype.txt | awk '{{print $1"\t"$2"\t"$4"\t"$3}}' > {output.j21} 2>> {log}
		rm genes_to_phenotype.txt 2>> {log}

		# this is techically from January 2017, but since the only other 2016 version I can find is from Jan 2016 which is right after and persumably very similar to the Dec 2015 data I am using this as the 2016 data.
		wget --no-check-certificate https://github.com/drseb/HPO-archive/raw/master/hpo.annotations.monthly/115/archive.zip 2>> {log}
		unzip archive.zip 2>> {log}
		cp archive/annotation/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt {output.j16} 2>> {log}
		rm -rf archive* 2>> {log}

		# 2020
		wget --no-check-certificate https://ci.monarchinitiative.org/view/hpo/job/hpo.annotations/4/artifact/rare-diseases/util/annotation/genes_to_phenotype.txt 2>> {log}
		cat  genes_to_phenotype.txt | awk '{{print $1"\t"$2"\t"$4"\t"$3}}' > {output.j20} 2>> {log}
		rm genes_to_phenotype.txt 2>> {log}
		
		# 2019, this is from the wayback machine, but is not wget able, I copy and pasted the info into a file to save
		# https://web.archive.org/web/20190902195916/http://compbio.charite.de:80/jenkins/job/hpo.annotations.monthly/lastSuccessfulBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt 2>> {log}
		cp Resources/jenkins2019.txt {output.j19} 2>> {log}

		# Infer the 2018 data from the current data (per recommendation of Peter Robinson (PI on HPO)
		wget https://ci.monarchinitiative.org/view/hpo/job/hpo.annotations/11/artifact/rare-diseases/misc/phenotype_annotation.tab 2>> {log}
		mv phenotype_annotation.tab Resources/ 2>> {log}
		python Scripts/infer_2018.py 2>> {log}
		cat genes_to_phenotype.txt | awk '{{print $1"\t"$2"\t"$4"\t"$3}}' > {output.j18} 2>> {log}
		rm genes_to_phenotype.txt 2>> {log}

		# Current date (2022 at the time of typing this)
		wget http://purl.obolibrary.org/obo/hp/hpoa/genes_to_phenotype.txt
		cat genes_to_phenotype.txt | awk '{{print $1"\t"$2"\t"$4"\t"$3}}' > {output.j22} 2>> {log}
                rm genes_to_phenotype.txt 2>> {log}
		"""

rule download_string:
	input:
		'work/start.start'
	output:
		s17='work/2017/string_edge_list.txt',
		s15='work/2015/string_edge_list.txt',
		s21='work/2021/string_edge_list.txt',
		s18='work/2018/string_edge_list.txt',
		s19='work/2019/string_edge_list.txt',
		s20='work/2020/string_edge_list.txt',
		s16='work/2016/string_edge_list.txt',
		s22='work/2022/string_edge_list.txt'
	shell:
		"""
		for x in $(seq 2015 2022);
		do
			mkdir -p work/$x
		done

		wget --no-check-certificate http://version10.string-db.org/download/protein.links.v10/9606.protein.links.v10.txt.gz
		gunzip 9606.protein.links.v10.txt.gz
		mv 9606.protein.links.v10.txt {output.s17}

		wget --no-check-certificate http://string91.embl.de/newstring_download/protein.links.v9.1/9606.protein.links.v9.1.txt.gz
		gunzip 9606.protein.links.v9.1.txt.gz
		mv 9606.protein.links.v9.1.txt {output.s15}

		wget --no-check-certificate https://stringdb-static.org/download/protein.links.v11.5/9606.protein.links.v11.5.txt.gz
		gunzip 9606.protein.links.v11.5.txt.gz
		mv 9606.protein.links.v11.5.txt {output.s21}
		# There has been no update since 2021
		cp {output.s21} {output.s22}

		# 2020 and 2019 are from different web pages 11b and 11 of string but route to the same list of protein links
		wget --no-check-certificate https://stringdb-static.org/download/protein.links.v11.0/9606.protein.links.v11.0.txt.gz
		gunzip 9606.protein.links.v11.0.txt.gz
		mv 9606.protein.links.v11.0.txt {output.s20}

		wget --no-check-certificate https://stringdb-static.org/download/protein.links.v11.0/9606.protein.links.v11.0.txt.gz
		gunzip 9606.protein.links.v11.0.txt.gz
		mv 9606.protein.links.v11.0.txt {output.s19}
		
		wget --no-check-certificate https://version-10-5.string-db.org/download/protein.links.v10.5/9606.protein.links.v10.5.txt.gz
		gunzip 9606.protein.links.v10.5.txt.gz
		mv 9606.protein.links.v10.5.txt {output.s18}

		wget --no-check-certificate http://version10.string-db.org/download/protein.links.v10/9606.protein.links.v10.txt.gz
		gunzip 9606.protein.links.v10.txt.gz
		mv 9606.protein.links.v10.txt {output.s17}

		wget --no-check-certificate http://version10.string-db.org/download/protein.links.v10/9606.protein.links.v10.txt.gz
		gunzip 9606.protein.links.v10.txt.gz
		mv 9606.protein.links.v10.txt {output.s16}
		"""

rule unzip_gene_info:
	input:
		'work/start.start'
	output:
		'BOCC/all_genes_info.json'
	shell:
		"""
		if [ -f gunzip BOCC/all_genes_info.json.gz ]
		then
   			gunzip BOCC/all_genes_info.json.gz
		else
   			echo "Doing nothing"
		fi
		"""

rule download_gnomad:
	input:
		'work/start.start'
	output:
		'work/Resources/gnomad.v2.1.1.all_lofs.txt'
	shell:
		"""
		wget https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-flagship-lof/v1.0/gnomad.v2.1.1.all_lofs.txt.bgz
		mv gnomad.v2.1.1.all_lofs.txt.bgz gnomad.v2.1.1.all_lofs.txt.gz
		bgzip -d gnomad.v2.1.1.all_lofs.txt.gz
		mv gnomad.v2.1.1.all_lofs.txt {output}	
		"""

rule collect_string_info:
	input:
		'work/start.start'
	output:
		'work/Resources/9606.protein.info.v11.5.txt'
	shell:
		"""
		mkdir -p work/Resources
		wget --no-check-certificate https://stringdb-static.org/download/protein.info.v11.5/9606.protein.info.v11.5.txt.gz
		gunzip 9606.protein.info.v11.5.txt.gz
		mv 9606.protein.info.v11.5.txt {output}
		"""

rule make_edgelists:
	input:
		string='work/{year}/string_edge_list.txt',
		jenkins='work/{year}/genes_to_phenotype.txt',
		hpo='work/{year}/hp.obo',
		string_info='work/Resources/9606.protein.info.v11.5.txt'
	output:
		'Edgelists/String_HPO_{year}.phenotypic_branch.edgelist.txt',
		'Edgelists/String_HPO_{year}.phenotypic_branch.numbered.edgelist.txt',
		'Edgelists/String_HPO_{year}.phenotypic_branch.nodenames.txt',
		'Edgelists/String_HPO_{year}.phenotypic_branch.nodeattributes.txt',
		'Edgelists/String_HPO_{year}.all_hpo.edgelist.txt',
		'Edgelists/String_HPO_{year}.all_hpo.numbered.edgelist.txt',
		'Edgelists/String_HPO_{year}.all_hpo.nodenames.txt',
		'Edgelists/String_HPO_{year}.all_hpo.nodeattributes.txt'
	shell:
		"""
		mkdir -p Edgelists
		python Scripts/create_edgelist.py {input.hpo} {input.jenkins} {input.string} {input.string_info} Edgelists/String_HPO_{wildcards.year}
		cat Edgelists/String_HPO_{wildcards.year}.phenotypic_branch.nodenames.txt  | python Scripts/nodenames_to_attributes.py > Edgelists/String_HPO_{wildcards.year}.phenotypic_branch.nodeattributes.txt
		cat Edgelists/String_HPO_{wildcards.year}.all_hpo.nodenames.txt  | python Scripts/nodenames_to_attributes.py > Edgelists/String_HPO_{wildcards.year}.all_hpo.nodeattributes.txt 
		"""	

"""
--- Clustering ---
"""


rule cluster_greedy:
	input:
		'Edgelists/String_HPO_{year}.phenotypic_branch.edgelist.txt'
	output:
		'Clusters/{year}/greedy.{year}.coms.txt'
	shell:
		"""
		python Algorithms/greedy.py --edgelist {input} --out {output}
		"""

rule cluster_walktrap:
	input:
		'Edgelists/String_HPO_{year}.phenotypic_branch.edgelist.txt'
	output:
		'Clusters/{year}/walktrap.{year}.coms.txt'
	shell:
		"""
		python Algorithms/walktrap.py --edgelist {input} --out {output}
		"""

rule cluster_infomap:
	input:
		el='Edgelists/String_HPO_{year}.phenotypic_branch.edgelist.txt',
		names='Edgelists/String_HPO_{year}.phenotypic_branch.nodenames.txt' 
	output:
		'Clusters/{year}/infomap.{year}.coms.txt'
	shell:
		"""
		python Algorithms/cluster_infomap.py --edgelist {input.el} --out {output} --nodenames {input.names}
		"""

rule install_cesna:
	input:
		'work/start.start',
	output:
		'Algorithms/cesna'
	shell:
		"""
		git clone git@github.com:snap-stanford/snap.git
		cd snap/examples/cesna
		make
		cp cesna ../../../Algorithms/
		cd ../../..
		rm -rf snap/
		"""

rule cluster_cesna:
	input:
		cesna='Algorithms/cesna',
		el='Edgelists/String_HPO_{year}.phenotypic_branch.numbered.edgelist.txt',
		el_names='Edgelists/String_HPO_{year}.phenotypic_branch.nodenames.txt',
		attr='Edgelists/String_HPO_{year}.phenotypic_branch.nodeattributes.txt'
	output:
		'Clusters/{year}/cesna.{year}.coms.txt'
	threads: 60
	shell:
		"""
		cp CachedClusters/cesna.{wildcards.year}.coms.txt {output}
		"""

# this is the actual censa code, leaving it here for later use
"""
		cd Clusters/{wildcards.year}	
		../../Algorithms/cesna -i:../../{input.el} -l:../../{input.el_names} -n:../../{input.attr} -a:../../{input.attr} -c:-1 -nt:60 -o:hpo_string_cesna -mc:10000 -xc:15000
		python ../../Scripts/num_coms.py hpo_string_cesnacmtyvv.txt > cesna.{wildcards.year}.coms.txt
		#python ../../Algorithms/number_cesna_results.py --cesna_res hpo_string_cesnacmtyvv.txt --output cesna.{wildcards.year}.coms.txt --node_names ../../{input.el_names}
		cd ../..
"""	

rule infer_all_jenkins_years:
	input:
		'work/start.start'
	output:
		ij15='work/2015/inferred_genes_to_phenotype.txt',
		ij16='work/2016/inferred_genes_to_phenotype.txt',
		ij17='work/2017/inferred_genes_to_phenotype.txt',
		ij18='work/2018/inferred_genes_to_phenotype.txt',
		ij19='work/2019/inferred_genes_to_phenotype.txt',
		ij20='work/2020/inferred_genes_to_phenotype.txt',
		ij21='work/2021/inferred_genes_to_phenotype.txt',
		ij22='work/2022/inferred_genes_to_phenotype.txt'
	shell:
		"""
		python Scripts/infer_any_year.py --year 2016 --output inferred_genes_to_phenotype.txt
		cat inferred_genes_to_phenotype.txt | awk '{{print $1"\t"$2"\t"$4"\t"$3}}' > {output.ij15}
		rm inferred_genes_to_phenotype.txt
		echo "15 done"
		
		python Scripts/infer_any_year.py --year 2017 --output inferred_genes_to_phenotype.txt
		cat inferred_genes_to_phenotype.txt | awk '{{print $1"\t"$2"\t"$4"\t"$3}}' > {output.ij16}
		rm inferred_genes_to_phenotype.txt
		echo "16 done"
		
		python Scripts/infer_any_year.py --year 2018 --output inferred_genes_to_phenotype.txt
		cat inferred_genes_to_phenotype.txt | awk '{{print $1"\t"$2"\t"$4"\t"$3}}' > {output.ij17}
		rm inferred_genes_to_phenotype.txt
		echo "17 done"
		
		python Scripts/infer_any_year.py --year 2019 --output inferred_genes_to_phenotype.txt
		cat inferred_genes_to_phenotype.txt | awk '{{print $1"\t"$2"\t"$4"\t"$3}}' > {output.ij18}
		rm inferred_genes_to_phenotype.txt
		echo "18 done"
		
		python Scripts/infer_any_year.py --year 2020 --output inferred_genes_to_phenotype.txt
		cat inferred_genes_to_phenotype.txt | awk '{{print $1"\t"$2"\t"$4"\t"$3}}' > {output.ij19}
		rm inferred_genes_to_phenotype.txt
		echo "19 done"
		
		python Scripts/infer_any_year.py --year 2021 --output inferred_genes_to_phenotype.txt
		cat inferred_genes_to_phenotype.txt | awk '{{print $1"\t"$2"\t"$4"\t"$3}}' > {output.ij20}
		rm inferred_genes_to_phenotype.txt
		echo "20 done"
		
		python Scripts/infer_any_year.py --year 2022 --output inferred_genes_to_phenotype.txt
		cat inferred_genes_to_phenotype.txt | awk '{{print $1"\t"$2"\t"$4"\t"$3}}' > {output.ij21}
		rm inferred_genes_to_phenotype.txt
		echo "21 done"

		python Scripts/infer_any_year.py --year 2023 --output inferred_genes_to_phenotype.txt
		cat inferred_genes_to_phenotype.txt | awk '{{print $1"\t"$2"\t"$4"\t"$3}}' > {output.ij22}
		echo "22 done"
		"""

rule paris_crosses:
	input:
		el='Edgelists/String_HPO_{year}.phenotypic_branch.edgelist.txt',
		com='Clusters/{year}/{algo}.{year}.coms.txt'
		# commented out to avoid needing to rerun cesna in the snakemake pipeline
		#com='CachedClusters/{algo}.{year}.coms.txt'
	output:
		'SubComs/{year}/paris.{algo}.{year}.coms.txt'
	shell:
		# this is the real code if you are reading this, put it back
		"""
		echo "paris " {input.el} {input.com}
		mkdir -p SubComs
		mkdir -p SubComs/{wildcards.year}
		python Scripts/hierarchical_clustering.py --edgelist {input.el}  --algo paris --coms {input.com} --output SubComs/{wildcards.year}/tmp.paris.{wildcards.algo}.{wildcards.year} --com_size 100
		python Scripts/aggregate_size.py SubComs/{wildcards.year} tmp.paris > {output}
		rm SubComs/{wildcards.year}/tmp.paris.*
		"""

"""
		cp RelabeledSubComs/{wildcards.year}/paris.{wildcards.algo}.{wildcards.year}.coms.txt {output}
"""
rule ward_crosses:
	input:
		el='Edgelists/String_HPO_{year}.phenotypic_branch.edgelist.txt',
		com='Clusters/{year}/{algo}.{year}.coms.txt'
		#com='CachedClusters/{algo}.{year}.coms.txt'
	output:
		'SubComs/{year}/ward.{algo}.{year}.coms.txt'
	shell:
		"""
                mkdir -p SubComs
                mkdir -p SubComs/{wildcards.year}
		echo "ward " {input.el} {input.com}
                python Scripts/hierarchical_clustering.py --edgelist {input.el}  --algo ward --coms {input.com} --output SubComs/{wildcards.year}/tmp.ward.{wildcards.algo}.{wildcards.year} --com_size 100
                python Scripts/aggregate_size.py SubComs/{wildcards.year} tmp.ward > {output}
                rm SubComs/{wildcards.year}/tmp.ward.*
		"""

"""
			cp RelabeledSubComs/{wildcards.year}/ward.{wildcards.algo}.{wildcards.year}.coms.txt {output}
"""
rule louvain_crosses:
	input:
		el='Edgelists/String_HPO_{year}.phenotypic_branch.edgelist.txt',
		com='Clusters/{year}/{algo}.{year}.coms.txt'
		#com='CachedClusters/{algo}.{year}.coms.txt'
	output:
		'SubComs/{year}/louvain.{algo}.{year}.coms.txt'
	shell:
		"""
		#cp RelabeledSubComs/{wildcards.year}/louvain.{wildcards.algo}.{wildcards.year}.coms.txt {output}
                mkdir -p SubComs
                mkdir -p SubComs/{wildcards.year}
                python Scripts/hierarchical_clustering.py --edgelist {input.el}  --algo louvain --coms {input.com} --output SubComs/{wildcards.year}/tmp.louvain.{wildcards.algo}.{wildcards.year} --com_size 100
                python Scripts/aggregate_size.py SubComs/{wildcards.year} tmp.louvain > {output}
                rm SubComs/{wildcards.year}/tmp.louvain.*
		"""
rule g2p_edge_lists:
	input:
		el='Edgelists/String_HPO_{year}.phenotypic_branch.edgelist.txt'
	output:
		'g2p_Edgelists/String_HPO_{year}.phenotypic_branch.g2p_edgelist.txt'
	shell:
		"""
		cat {input.el} | python Scripts/get_g2p_edges.py > {output}
		"""



rule BOCC_results:
	input:
		coms='SubComs/{year}/paris.{algo}.{year}.coms.txt',
		el='Edgelists/String_HPO_{year}.phenotypic_branch.edgelist.txt',
		plof='work/Resources/gnomad.v2.1.1.all_lofs.txt'
	output:
		'ResultsBOCC/{year}/paris.{algo}.{year}.bocc_res.tsv'
	shell:
		"""
		mkdir -p ResultsBOCC/{wildcards.year}
		python BOCC/summarize_clusters.py --coms {input.coms} --alpha 0.05 --mg2 unused_place_holder --graph {input.el} --gnomad_file {input.plof} --out {output}
		"""

#rule all_BOCC_results:
#	input:
#                coms='SubComs/{year}/{hier}.{algo}.{year}.coms.txt',
#                el='Edgelists/String_HPO_{year}.phenotypic_branch.edgelist.txt',
#                plof='work/Resources/gnomad.v2.1.1.all_lofs.txt'
#	output:
#                'ResultsBOCC/{year}/{hier}.{algo}.{year}.bocc_res.tsv'
#        shell:
#                """
#                python BOCC/summarize_clusters.py --coms {input.coms} --alpha 0.05 --mg2 unused_place_holder --graph {input.el} --gnomad_file {input.plof} --out {output}
#                """

rule snowball:
        input:
                el='Edgelists/String_HPO_{year}.phenotypic_branch.edgelist.txt',
                coms='SubComs/{year}/{hier}.{algo}.{year}.coms.txt',
		g2p=expand('g2p_Edgelists/String_HPO_{year}.phenotypic_branch.g2p_edgelist.txt',year=YEARS_short)
        output:
                "SnowballResults/snowball.{hier}.{algo}.String_HPO_{year}.phenotypic_branch.tsv"
        shell:
                """
		# add one to the year for get the edgelist for scoring
                y=$(({wildcards.year}+1))
		if [ $y -eq "2022" ]
		then
			echo "Nada, next year does not exist" > {output}
		else
			# CachedSnowballResults
			if test -f "CachedSnowballResults/snowball.{wildcards.hier}.{wildcards.algo}.String_HPO_{wildcards.year}.phenotypic_branch.tsv";
			then
				cp CachedSnowballResults/snowball.{wildcards.hier}.{wildcards.algo}.String_HPO_{wildcards.year}.phenotypic_branch.tsv {output}
			else
                		python Scripts/snowball.py --edgelist {input.el} --output {output} --coms {input.coms} --new_edges g2p_Edgelists/String_HPO_${{y}}.phenotypic_branch.g2p_edgelist.txt --reps 100
			fi
		fi
                """

rule snowball_mouse:
	input:
		el='Edgelists/String_HPO_2021.phenotypic_branch.edgelist.txt',
		coms='SubComs/2021/{hier}.{algo}.2021.coms.txt',
		mouse_g2p='Resources/hpo_to_gene_derived_by_mpo.edgelist.2021.unique.txt'
	output:
		"MouseSnowballResults/snowball.{hier}.{algo}.String_HPO_2021.mouse.phenotypic_branch.tsv"
	shell:
		"""
		mkdir -p MouseSnowballResults
		python Scripts/snowball.py --edgelist {input.el} --output {output} --coms {input.coms} --new_edges {input.mouse_g2p} --reps 100
		"""

rule make_drug_list:
	input:
	output:
		"Resources/new_drug_edges.txt"
	shell:
		"""
		python Scripts/create_drug_discovery_list.py
		"""

rule snowball_drugs:
	input:
		el='Edgelists/String_HPO_2021.phenotypic_branch.edgelist.txt',
		coms='SubComs/2022/{hier}.{algo}.2022.coms.txt',
		mouse_g2p='Resources/new_drug_edges_sampled.txt'
	output:
		"DrugRediscovery/snowball.{hier}.{algo}.String_HPO_2022.phenotypic_branch.tsv"
	shell:
		"""
		mkdir -p DrugRediscovery
		python Scripts/snowball_fixed.py --edgelist {input.el} --output {output} --coms {input.coms} --new_edges {input.mouse_g2p} --reps 100
		"""

rule snowball_drugs_results:
	input:
		'DrugRediscovery/snowball.{hier}.{algo}.String_HPO_2022.phenotypic_branch.tsv'
	output:
		'DrugRediscoveryResults/snowball.results.{hier}.{algo}.String_HPO_2022.phenotypic_branch.tsv'
	shell:
		"""
		python Scripts/add_snowballing_drug_results.py {input} > {output}
		"""

rule snowball_fixed:
        input:
                el='Edgelists/String_HPO_{year}.phenotypic_branch.edgelist.txt',
                coms='SubComs/{year}/paris.{algo}.{year}.coms.txt',
		g2p=expand('g2p_Edgelists/String_HPO_{year}.phenotypic_branch.g2p_edgelist.txt',year=YEARS_short)
        output:
                "SnowballResultsFixed/snowball.paris.{algo}.String_HPO_{year}.phenotypic_branch.tsv"
        shell:
                """
                mkdir -p SnowballResultsFixed
		y=$(({wildcards.year}+1))
                if [ $y -eq "2023" ]
		then
			echo "Nada, next year does not exist" > {output}
		else
                	# python Scripts/snowball_fixed.py --edgelist {input.el} --output {output} --coms {input.coms} --new_edges g2p_Edgelists/String_HPO_${{y}}.phenotypic_branch.g2p_edgelist.txt --reps 100
			python Scripts/snowball_g_p_specific.py --edgelist {input.el} --output {output} --coms {input.coms} --new_edges g2p_Edgelists/String_HPO_${{y}}.phenotypic_branch.g2p_edgelist.txt --reps 100
		fi
                """

rule create_bocc_results:
	input:
		expand('SubComs/{year}/{hier}.{algo}.{year}.coms.txt',algo=ALGOS,year=YEARS,hier=HIER_ALGOS)	
	output:
		expand('CleanedSplit100BOCCResultsCombined/{hier}.{algo}.{year}.bocc_res.tsv',algo=ALGOS,year=YEARS,hier=HIER_ALGOS)
	shell:
		"""
		touch .start.txt
		snakemake -s prepare_subcoms.snake --cores 10
		snakemake -s generate_bocc_results.snake --cores 20
		"""

rule add_snowballing_to_bocc_res:
	input:
		bocc='CleanedSplit100BOCCResultsCombined/{hier}.{algo}.{year}.bocc_res.tsv',
		snow="SnowballResultsFixed/snowball.{hier}.{algo}.String_HPO_{year}.phenotypic_branch.tsv"
	output:
		'SnowballedCleanedSplit100BOCCResultsCombinedFixed/{year}/{hier}.{algo}.{year}.bocc_res.tsv'
	shell:
		"""
		mkdir -p SnowballedCleanedSplit100BOCCResultsCombinedFixed/{wildcards.year}
		python Scripts/add_snowballing_to_bocc.py {input.bocc} {input.snow} {output}
		"""

rule add_any_new_edge_count_to_bocc_res:
	input:
		bocc='SnowballedCleanedSplit100BOCCResultsCombinedFixed/{year}/{hier}.{algo}.{year}.bocc_res.tsv',
		subcom='SubComs/{year}/{hier}.{algo}.{year}.coms.txt',
		g2p=expand('g2p_Edgelists/String_HPO_{year}.phenotypic_branch.g2p_edgelist.txt',year=YEARS_short)
	output:
		'AnyNewEdgeCountSnowballedCleanedSplit100BOCCResultsCombinedFixed/{year}/{hier}.{algo}.{year}.bocc_res.tsv'
	shell:
		"""
		mkdir -p AnyNewEdgeCountSnowballedCleanedSplit100BOCCResultsCombinedFixed/{wildcards.year}
		y=$(({wildcards.year}+1))
		python Scripts/add_any_new_edges_to_bocc.py {input.bocc} {input.subcom} g2p_Edgelists/String_HPO_${{y}}.phenotypic_branch.g2p_edgelist.txt {output}
		"""

rule check_avg_internal_degree_of_parent_clusters:
	input:
		el='Edgelists/String_HPO_{year}.phenotypic_branch.edgelist.txt',
		com='Clusters/{year}/{algo}.{year}.coms.txt'
	output:
		'AverageInteralDegreeCheckResults/{algo}.{year}.coms.txt'
	shell:
		"""
		mkdir -p AverageInteralDegreeCheckResults
		python Scripts/average_interal_degree.py {input.com} {input.el} > {output}
		"""



