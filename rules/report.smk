rule copy_config:
	"""
	Copy config advanced in 00_logs
	"""
	output:
		"results/00_logs/config_advanced"
	params:
		config=config["repo_script"]
	threads: get_thread
	resources:
		mem_mb=get_mem
	params:
		wd=WORKDIR
	shell:
		"""
		cp {params.config}/config_advanced {params.wd}/{output}
		"""

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

def get_input_report(wildcards):
	list=[]
	list.append("results/01_sequence_qc/"+NAME_PROJECT+"_multiqc_trim.html")
	list.extend(expand("results/02_mapping/depth/{s.sample}.depth",s=SAMPLE.itertuples()),)
	list.extend(expand("results/02_mapping/coverage/{s.sample}.coverage",s=SAMPLE.itertuples()),)
	list.extend(expand("results/02_mapping/flagstat/{s.sample}.flagstat",s=SAMPLE.itertuples()),)
	list.append(config["sample"])
	list.append("results/01_sequence_qc/log/trimmomatic.log")
	list.append("results/00_logs/config_advanced")
	list.append("results/03_snv_indels_calling/"+NAME_PROJECT+"_snv_indel.vcf")
	list.append("results/03_snv_indels_calling/"+NAME_PROJECT+"snpeff.html")
	if get_vector(config,"regions")[0] == "TRUE":
		list.append("results/genome/regions.bed")
	if get_vector(config,"vector")[0] == "TRUE":
		checkpoint_output=checkpoints.cut_vector_file.get(**wildcards).output[0]
		list.extend(expand("results/05_tdnascan/{s.sample}/{v}/5.{v}_insertion.annotated.new.bed", s=SAMPLE.itertuples(),v=glob_wildcards(os.path.join(checkpoint_output, "{v}.fa")).v))
		list.extend(expand("results/04_pindel/{s.sample}_{ext_pindel}",s=SAMPLE.itertuples(),ext_pindel=["BP","CloseEndMapped","D","INV","LI","RP","SI","TD"]),)
		list.append("results/genome/TDNA_sequence.fasta")
	if ( (get_vector(config,"vector")[0] == "FALSE")  & (config["SV"] == "TRUE" )) :
		list.extend(expand("results/04_pindel/{s.sample}_{ext_pindel}",s=SAMPLE.itertuples(),ext_pindel=["BP","CloseEndMapped","D","INV","LI","RP","SI","TD"]),) 
	return list

rule report:
	"""
	Create final report comporting all analysis results
	"""
	input:
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
		DP=get_DP,
		AR=get_AR,
	conda:
		"../envs/R.yaml"
	threads: get_thread
	resources:
		mem_mb=get_mem
	shell:
		"Rscript -e \"rmarkdown::render('{params.workdir}/script/ms_report.Rmd', output_file = '{params.workdir}/{output}', params = list(result_dir='{params.workdir}/results/', DP.min='{params.DP}', AR.min='{params.AR}'))\""

