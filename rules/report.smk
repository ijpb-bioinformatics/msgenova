rule copy_config:
	"""
	Copy config advanced in 00_logs
	"""
	output:
		"results/00_logs/config_advanced"
	threads: get_thread
	resources:
		mem_mb=get_mem
	params:
		wd=WORKDIR,
		config=config["repo_script"]
	run:
		with open(output[0],'w') as yaml_file:
			yaml.dump(config, yaml_file,default_flow_style=False)


#	shell:
#		"""
#		cp {params.config}/config_advanced {params.wd}/{output}
#		"""
#

def get_regions_file(wildcards):
	return config["regions"]

rule copy_region:
	"""
	If a region file is present, copy region file
	"""
	input:
		get_regions_file
	output:
		"results/genome/regions.bed"
	threads: get_thread
	resources:
		mem_mb=get_mem
	params:
		wd=WORKDIR
	shell:
		"""
		cp {input} {params.wd}/{output}
		"""

rule copy_sample_file:
	"""
	Copy sample file for report creation
	"""
	input:
		config["sample"]
	output:
		"results/genome/sample_sheet"
	threads: get_thread
	resources:
		mem_mb=get_mem
	params:
		wd=WORKDIR
	shell:
		"""
		cp {input} {params.wd}/{output}
		"""

rule extract_flagstat:
	"""
	Concatenate flagstat output to simplify report
	"""
	input:
		expand("results/02_mapping/flagstat/{s.sample}.flagstat",s=SAMPLE.itertuples()),
	output:
		"results/02_mapping/flagstat/concatenate_flagstat.txt"
	threads: get_thread
	resources:
		mem_mb=get_mem
	params:
		wd=WORKDIR
	shell:
		"""
			for file in {input}
			do
				name=`basename {params.wd}/$file | sed 's/\.flagstat//g' `
				total_reads=`cat {params.wd}/$file |grep 'in\ total' | head -n1 | cut -d" " -f 1`
				total_mapped=`cat {params.wd}/$file | grep 'mapped' | head -n1 | cut -d" " -f 1`
				supp=`cat {params.wd}/$file | grep 'supplementary' | head -n1 | cut -d" " -f 1`
				echo -e $name'\t'$total_reads'\t'$total_mapped'\t'$supp >> {params.wd}/{output}
			done
		"""

def get_input_report(wildcards):
	list=[]
	list.append("results/01_sequence_qc/"+NAME_PROJECT+"_multiqc_trim.html")
	list.extend(expand("results/02_mapping/depth/{s.sample}.depth",s=SAMPLE.itertuples()),)
	list.extend(expand("results/02_mapping/coverage/{s.sample}.coverage",s=SAMPLE.itertuples()),)
	list.extend(expand("results/02_mapping/flagstat/{s.sample}.flagstat",s=SAMPLE.itertuples()),)
	list.append("results/02_mapping/flagstat/concatenate_flagstat.txt")
	list.append("results/genome/sample_sheet")
	list.append("results/01_sequence_qc/log/trimmomatic.log")
	list.append("results/00_logs/config_advanced")
	list.append("results/03_snv_indels_calling/"+NAME_PROJECT+"_snv_indel.vcf")
	list.append("results/03_snv_indels_calling/"+NAME_PROJECT+"snpeff.html")
	if get_vector(config,"regions")[0] == "TRUE":
		list.append("results/genome/regions.bed")
	if get_vector(config,"vector")[0] == "TRUE":
		checkpoint_output=checkpoints.cut_vector_file.get(**wildcards).output[0]
		list.extend(expand("results/05_tdnascan/{s.sample}/{v}/5.{v}_insertion.annotated.bed", s=SAMPLE.itertuples(),v=glob_wildcards(os.path.join(checkpoint_output, "{v}.fa")).v))
		list.append("results/04_pindel/"+NAME_PROJECT+"_sv.vcf.gz")
		#list.extend(expand("results/04_pindel/{s.sample}_{ext_pindel}",s=SAMPLE.itertuples(),ext_pindel=["BP","CloseEndMapped","D","INV","LI","RP","SI","TD"]),)
		list.append("results/genome/TDNA_sequence.fasta")
	if ( (get_vector(config,"vector")[0] == "FALSE")  & (config["SV"] == "TRUE" )) :
		list.append("results/04_pindel/"+NAME_PROJECT+"_sv.vcf.gz")
		#list.extend(expand("results/04_pindel/{s.sample}_{ext_pindel}",s=SAMPLE.itertuples(),ext_pindel=["BP","CloseEndMapped","D","INV","LI","RP","SI","TD"]),) 
	return list

rule report:
	"""
	Create final report comporting all analysis results
	"""
	input:
		#"results/04_pindel/pindel_results.vcf.gz"
		get_input_report
		#multiqc="results/01_sequence_qc/"+NAME_PROJECT+"_multiqc_trim.html",
		#file_depth=expand("results/02_mapping/depth/{s.sample}.depth",s=SAMPLE.itertuples()),
		#flag="results/03_mapping/flagstat/Extract_data.flagstat",
		#flag=expand("results/02_mapping/flagstat/{s.sample}.flagstat",s=SAMPLE.itertuples()),
		#sample_file=config["sample"],
		#file_coverage=expand("results/02_mapping/coverage/{s.sample}.coverage",s=SAMPLE.itertuples()),
		#expand("results/02_mapping/coverage/{s.sample}.coverage",s=SAMPLE.itertuples()),
		#reference=config["reference"],
		#qc_trimmo="results/01_sequence_qc/log/trimmomatic.log",
		#haplotype_caller_html="results/03_snv_indels_calling/"+NAME_PROJECT+"snpeff.html",
		#haplotype_caller_vcf="results/03_snv_indels_calling/"+NAME_PROJECT+"_snv_indel.vcf",
		#config="results/00_logs/config_advanced",
		#insertion=aggregate_vector,
		#pindel=expand("results/04_pindel/{s.sample}_{ext_pindel}",s=SAMPLE.itertuples(),ext_pindel=["BP","CloseEndMapped","D","INV","LI","RP","SI","TD"]),
	output:
		"results/06_report/msgenova_report.html"
	params:
		workdir=config["repo_script"],
		wd=WORKDIR,
		DP=get_DP,
		AR=get_AR,
	conda:
		"../envs/R.yaml"
	threads: get_thread
	resources:
		mem_mb=get_mem
	shell:
		"Rscript -e \"rmarkdown::render('{params.workdir}/script/ms_report.Rmd', output_file = '{params.wd}/{output}', params = list(result_dir='{params.wd}/results/', DP.min='{params.DP}', AR.min='{params.AR}'))\""
		# "Rscript -e \"rmarkdown::render('/save/project/ijpb/bioinfo-code/src/essai_report.Rmd', output_file = '{params.wd}/{output}', params = list(result_dir='/work/gadam/msgenova_reduce/results/', DP.min='5', AR.min='0.2'),intermediates_dir='{params.wd}/results/')\""
