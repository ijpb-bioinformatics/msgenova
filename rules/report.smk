rule copy_config:
	"""
	Copy config advanced in 00_logs
	"""
	output:
		"results/00_logs/config_advanced"
	threads: get_thread("copy_config")
	resources:
		mem_mb=get_mem("copy_config")
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
		"results/00_logs/regions.bed"
	threads: get_thread("copy_region")
	resources:
		mem_mb=get_mem("copy_region")
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
		"results/00_logs/sample_sheet"
	threads: get_thread("copy_sample_file")
	resources:
		mem_mb=get_mem("copy_sample_file")
	params:
		wd=WORKDIR
	shell:
		"""
		cp {input} {params.wd}/{output}
		"""

rule copy_environments:
	output:
		directory("results/00_logs/envs/")
	threads: get_thread("copy_environments")
	resources: 
		mem_mb=get_mem("copy_environments")
	params:
		install_dir=config["repo_script"]
	shell:
		"""
		mkdir results/00_logs/envs
		cp {params.install_dir}/envs/env* results/00_logs/envs/
		"""

rule extract_flagstat:
	"""
	Concatenate flagstat output to simplify report
	"""
	input:
		expand("results/02_mapping/flagstat/{s.sample}.flagstat",s=SAMPLE.itertuples()),
	output:
		"results/02_mapping/flagstat/concatenate_flagstat.txt"
	threads: get_thread("extract_flagstat")
	resources:
		mem_mb=get_mem("extract_flagstat")
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

#		


# list.append("results/04_delly/"+NAME_PROJECT+"_sv_delly.vcf.gz")

def get_input_report(wildcards):
	list=[]
	list.append("results/01_sequence_qc/"+NAME_PROJECT+"_multiqc_trim.html")
	list.extend(expand("results/02_mapping/depth/{s.sample}.depth",s=SAMPLE.itertuples()),)
	list.extend(expand("results/02_mapping/coverage/{s.sample}.coverage",s=SAMPLE.itertuples()),)
	list.extend(expand("results/02_mapping/flagstat/{s.sample}.flagstat",s=SAMPLE.itertuples()),)
	list.extend(expand("results/04_delly/{s.sample}.vcf.gz",s=SAMPLE.itertuples()),)
	list.append("results/02_mapping/flagstat/concatenate_flagstat.txt")
	list.append("results/00_logs/sample_sheet")
	list.append("results/01_sequence_qc/log/trimmomatic.log")
	list.append("results/00_logs/config_advanced")
	list.append("results/03_snv_indels_calling/"+NAME_PROJECT+"_snv_indel.vcf")
	list.append("results/03_snv_indels_calling/"+NAME_PROJECT+"snpeff.html")
	list.append(directory("results/00_logs/envs/"))
	if get_vector(config,"regions")[0] == "TRUE":
		list.append("results/00_logs/regions.bed")
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
		get_input_report
	output:
		"results/msgenova_report.html"
	params:
		workdir=config["repo_script"],
		wd=WORKDIR,
		DP_min=get_DP_min,
		AR_min=get_AR_min,
		AR=get_AR
	conda:
		"../envs/env_R.yaml"
	threads: get_thread("report")
	resources:
		mem_mb=get_mem("report")
	shell:
		"Rscript -e \"rmarkdown::render('{params.workdir}/script/report.Rmd', output_file = '{params.wd}/{output}', intermediates_dir='{params.wd}/results',  params = list(result_dir='{params.wd}/results/', DP.min='{params.DP_min}', AR.min='{params.AR_min}',qual.min='30', projectName='NAME_PROJECT', AR='{params.AR}'))\""
		# "Rscript -e \"rmarkdown::render('/save/project/ijpb/bioinfo-code/src/essai_report.Rmd', output_file = '{params.wd}/{output}', params = list(result_dir='/work/gadam/msgenova_reduce/results/', DP.min='5', AR.min='0.2'),intermediates_dir='{params.wd}/results/')\""
