"""
Using the sub com files split into 100 com chunks from the previous pipeline, generate the bocc results for each one
Only 20 jobs can ever run at a time (the API will freak out otherwise) and if it takes more than 1 hours something is wrong and the jobs shoul dbe killed and restarted.
"""
YEARS = ['2015','2016','2017','2018','2019','2020','2021']
YEARS = ['2022','2021','2020','2019','2017','2016','2015','2018']
ALGOS = ['cesna','greedy','walktrap','infomap']
HIERS = ['ward','paris','louvain']

import os

shards = []
for y in YEARS:
	for f in os.listdir('Split100SubComs/{y}/'.format(y=y)):
		if '.' == f[0]:
			continue
		shards.append(f.replace('.coms.txt',''))	

rule all:
        input:  
                expand('Split100BOCCResults/{s}.bocc_res.tsv',s=shards),
		'done.txt',
		'clean_bocc_res.done'
rule get_bocc_results:
	input:
		'Split100SubComs/{s}.coms.txt'
	output:
		'Split100BOCCResults/{s}.bocc_res.tsv'
	log:
		'Logs/Split100BOCCResults/{s}.log'
	resources:
		time='1:00:00'
	shell:
		"""
		IN="{wildcards.s}"
		IFS='.' read -r -a array <<< "$IN"
		year="${{array[2]}}"

		echo $year > {log}
		python BOCC/summarize_clusters.py --coms {input} --alpha 0.05 --mg2 unused_place_holder --graph Edgelists/String_HPO_${{year}}.phenotypic_branch.edgelist.txt --gnomad_file work/Resources/gnomad.v2.1.1.all_lofs.txt --out {output} >> {log}
		"""

# cat Split100BOCCResults/${{h}}.${{a}}.${{y}}.*.bocc_res.tsv > Split100BOCCResultsCombinedi/${{h}}.${{a}}.${{y}}.bocc_res.tsv

rule combine_shards:
	input:
		expand('Split100BOCCResults/{s}.bocc_res.tsv',s=shards)
	output:
		'done.txt'
	shell:
		"""
		touch Split100BOCCResultsCombined
		declare -a ys=("2015" "2016" "2017" "2018" "2019" "2020" "2021" "2022")
		declare -a as=("greedy" "walktrap" "cesna" "infomap")
		declare -a hs=("paris" "ward" "louvain")

		for y in "${{ys[@]}}"
		do
			for a in "${{as[@]}}"
			do
				for h in "${{hs[@]}}"
				do
					cat Split100BOCCResults/${{h}}.${{a}}.${{y}}.*.bocc_res.tsv > Split100BOCCResultsCombined/${{h}}.${{a}}.${{y}}.bocc_res.tsv
				done
			done
		done
		touch done.txt
		"""

rule clean_bocc_res:
	input:
		'done.txt'
	output:
		'clean_bocc_res.done'
	shell:
		"""
		mkdir -p CleanedSplit100BOCCResultsCombined
		declare -a ys=("2015" "2016" "2017" "2018" "2019" "2020" "2021" "2022")
                declare -a as=("greedy" "walktrap" "cesna" "infomap")
                declare -a hs=("paris" "ward" "louvain")

                for y in "${{ys[@]}}"
                do
                        for a in "${{as[@]}}"
                        do
                                for h in "${{hs[@]}}"
                                do
					python Scripts/remove_extra_headers.py Split100BOCCResultsCombined/${{h}}.${{a}}.${{y}}.bocc_res.tsv CleanedSplit100BOCCResultsCombined/${{h}}.${{a}}.${{y}}.bocc_res.tsv
                                done
                        done
                done
		touch clean_bocc_res.done
		"""
