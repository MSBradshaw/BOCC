"""
Take the subcoms generaged by the first pipeline, relabel them, then split them into small digestable chunks for use in the next pipeline
Check the QC report: Renumbering_QC_Reports/all_reports.txt should not have any lines in it. If it does it means some of the communities in a file have the same id, which will RUIN the BOCC results

Example usage:
snakemake -s prepare_subcoms.snake --cores 10
"""
YEARS = ['2015','2016','2017','2018','2019','2020','2021','2022']
ALGOS = ['cesna','greedy','walktrap','infomap']
HIERS = ['ward','paris','louvain']


rule all:
	input:
		'.start.txt',
		expand('RelabeledSubComs/{y}/{h}.{a}.{y}.coms.txt',a=ALGOS,y=YEARS,h=HIERS),
		expand('Renumbering_QC_Reports/{h}.{a}.{y}.report',a=ALGOS,y=YEARS,h=HIERS),
		'Renumbering_QC_Reports/all_reports.txt',
		expand('Split100SubComs/{y}/.{h}.{a}.{y}.completed',a=ALGOS,y=YEARS,h=HIERS)
	
	

rule renumber_subcoms:
	input:
		'SubComs/{y}/{h}.{a}.{y}.coms.txt'
	output:
		'RelabeledSubComs/{y}/{h}.{a}.{y}.coms.txt'		
	shell:
		"""
		python Scripts/relabel_sub_com.py {input} {output}
		"""

rule qc_renumbering:
	input:
		'RelabeledSubComs/{y}/{h}.{a}.{y}.coms.txt'
	output:
		'Renumbering_QC_Reports/{h}.{a}.{y}.report'
	shell:
		"""
		mkdir -p Renumbering_QC_Reports
		python Scripts/qc_renumbering.py {input} > {output}
		"""

rule combine_qc_renumbering_reports:
	input:
		expand('Renumbering_QC_Reports/{h}.{a}.{y}.report',a=ALGOS,y=YEARS,h=HIERS)
	output:
		'Renumbering_QC_Reports/all_reports.txt'
	shell:
		"""
		cat {input} > {output}
		"""

rule split_coms:
	input:
		'RelabeledSubComs/{y}/{h}.{a}.{y}.coms.txt'
	output:
		indicatorfile='Split100SubComs/{y}/.{h}.{a}.{y}.completed'
	resources:
		threads=20
	shell:
		"""
		mkdir -p Split100SubComs
		mkdir -p Split100SubComs/{wildcards.y}
		split -l 100 {input} Split100SubComs/{wildcards.y}/{wildcards.h}.{wildcards.a}.{wildcards.y}. --additional-suffix=.coms.txt
		touch {output.indicatorfile}
		"""




